from __future__ import annotations
import csv
from pathlib import Path
from typing import NewType
from collections import namedtuple
from collections.abc import Iterable
from logging import getLogger, Logger

import numpy as np
import numpy.typing as npt
import numpy.lib.recfunctions as rcf

from .data import Data


# A haplotype block consists of
# 1) pop   - A population label (str), like 'YRI'
# 2) chrom - A chromosome name (str), like 'chr19' or simply '19'
# 3) bp    - The end position of the block in bp (int), like 1001038
# 4) cm    - The end position of the block in cM (float), like 43.078
HapBlock = [("pop", "U6"), ("chrom", "U10"), ("bp", np.uint32), ("cm", np.float64)]

# This tuple lists the haplotype blocks in a sample, one set for each chromosome
# Let's define a type alias, "SampleBlocks", for future use...
SampleBlocks = NewType(
    "SampleBlocks", "list[npt.NDArray[HapBlock], npt.NDArray[HapBlock]]]"  # type: ignore
)


class Breakpoints(Data):
    """
    A class for processing breakpoints from a file

    Attributes
    ----------
    data : dict[str, SampleBlocks]
        The haplotype blocks for each chromosome in each sample
        This dict maps samples (as strings) to their haplotype blocks (as SampleBlocks)
    fname : Path | str
        The path to the file containing the data
    labels : dict | None
        A dictionary containing population labels. It maps each label to the unique
        integers in the "pop" field of :py:attr:`~.Breakpoints.data`
    log: Logger
        A logging instance for recording debug statements.

    Examples
    --------
    >>> breakpoints = Breakpoints.load('tests/data/test.bp')
    """

    def __init__(self, fname: Path | str, log: Logger = None):
        super().__init__(fname, log)
        self._ext = "bp"
        self.labels = None

    @classmethod
    def load(
        cls: Breakpoints, fname: Path | str, samples: set[str] = None
    ) -> Breakpoints:
        """
        Load breakpoints from a TSV file

        Read the file contents and standardize the Breakpoints

        Parameters
        ----------
        fname
            See documentation for :py:attr:`~.Data.fname`
        samples : set[str], optional
            See documentation for :py:meth:`~.Breakpoints.read`

        Returns
        -------
        Breakpoints
            A Breakpoints object with the data loaded into its properties
        """
        breakpoints = cls(fname)
        breakpoints.read(samples)
        return breakpoints

    def read(self, samples: set[str] = None):
        """
        Read breakpoints from a TSV file into a data structure stored in :py:attr:`~.Breakpoints.data`

        Parameters
        ----------
        samples : set[str], optional
            A subset of the samples for which to extract breakpoints

            Defaults to loading breakpoints for all samples

        Raises
        ------
        AssertionError
            If the provided file doesn't follow the expected format
        """
        super().read()
        # call self.__iter__() to obtain the data
        self.data = dict(self.__iter__(samples))
        self.log.info(f"Loaded {len(self.data)} samples from .{self._ext} file")

    def __iter__(self, samples: set[str] = None) -> Iterable[str, SampleBlocks]:
        """
        Read breakpoints from a TSV line by line without storing more than a single
        sample at a time

        Parameters
        ----------
        samples : set[str], optional
            See documentation for :py:meth:`~.Breakpoints.read`

        Returns
        ------
        Iterable[str, SampleBlocks]
            An iterator over each sample in the file, where the sample if specified
            first as a string, and then followed by its SampleBlocks
        """
        # TODO: add a region parameter
        bps = self.hook_compressed(self.fname, mode="r")
        bp_text = csv.reader(bps, delimiter="\t")
        samp = None
        blocks = {}
        for line in bp_text:
            # ignore all of the comment lines
            if line[0].startswith("#"):
                continue
            # check: is this a sample name or a hap block line?
            if len(line) == 1:
                line = line[0]
                # get the end of the sample name (ex: "_1" or "_2")
                strand_num = line.rsplit("_", 1)[-1]
                if not (strand_num in ("1", "2")):
                    self.log.warning(
                        f"Ignoring improperly formatted line in bp file: '{line}'"
                    )
                    continue
                strand_num = int(strand_num) - 1
                if strand_num:
                    if samp is None:
                        self.log.error(f"bp file does not start with first strand")
                    elif samp != line[:-2]:
                        self.log.error(
                            f"'{samp}_1' was not followed by '{samp}_2' in the bp file"
                        )
                else:
                    if samp is not None and (samples is None or samp in samples):
                        # output the previous sample
                        yield samp, [np.array(b, dtype=HapBlock) for b in blocks]
                    samp = line[:-2]
                    blocks = [[], []]
            elif len(line) == 4:
                blocks[strand_num].append(tuple(line))
            else:
                self.log.warning(
                    f"Ignoring improperly formatted line in bp file: '{line}'"
                )
        if samp is not None and (samples is None or samp in samples):
            # output the previous sample
            yield samp, [np.array(b, dtype=HapBlock) for b in blocks]
        bps.close()

    def encode(self, labels: tuple[str] = None):
        """
        Replace each ancestral label in :py:attr:`~.Breakpoints.data` with an
        equivalent integer. Store a dictionary mapping these integers back to their
        respective labels in :py:attr:`~.Breakpoints.labels`.

        This method modifies :py:attr:`~.Breakpoints.data` in place.

        Parameters
        ----------
        labels: tuple[str], optional
            A list of population labels. The order of the labels in this list will be
            kept in the respective labels.
        """
        if not (self.labels is None):
            raise ValueError("The data has already been encoded.")
        # save the order of the fields for later reordering
        names = [f[0] for f in HapBlock]
        # initialize labels dict and label counter
        if labels is None:
            labels = {}
        else:
            labels = {pop: i for i, pop in enumerate(labels)}
        pop_count = len(labels)
        seen = set()
        for sample, blocks in self.data.items():
            for strand_num in range(len(blocks)):
                # initialize and fill the array of integers
                ints = np.zeros(len(blocks[strand_num]), dtype=[("pop", np.uint8)])
                for i, pop in enumerate(blocks[strand_num]["pop"]):
                    if pop not in labels:
                        labels[pop] = pop_count
                        pop_count += 1
                    ints[i] = labels[pop]
                    seen.add(pop)
                # replace the "pop" labels
                arr = rcf.drop_fields(blocks[strand_num], ["pop"])
                blocks[strand_num] = rcf.merge_arrays((arr, ints), flatten=True)[names]
        self.labels = {k: v for k, v in labels.items() if k in seen}

    def recode(self):
        """
        Replace each integer in :py:attr:`~.Breakpoints.data` with an
        equivalent ancestral label. Use the dictionary mapping these integers back to
        their respective ancestral labels stored in :py:attr:`~.Breakpoints.labels`.

        This method modifies :py:attr:`~.Breakpoints.data` in place.
        """
        if self.labels is None:
            raise ValueError("The data has already been recoded.")
        dtype = dict(HapBlock)
        names = list(dtype.keys())
        dtype = dtype["pop"]
        map_func = np.vectorize({v: k for k, v in self.labels.items()}.get)
        for sample, blocks in self.data.items():
            for strand_num in range(len(blocks)):
                # initialize and fill the array of pop labels
                pops = map_func(blocks[strand_num]["pop"]).astype([("pop", dtype)])
                # replace the "pop" labels
                arr = rcf.drop_fields(blocks[strand_num], ["pop"])
                blocks[strand_num] = rcf.merge_arrays((arr, pops), flatten=True)[names]
        self.labels = None

    @staticmethod
    def _find_blocks(
        blocks: npt.NDArray[np.uint32], positions: npt.NDArray[np.uint32]
    ) -> npt.NDArray[np.uint32]:
        """
        For each position in the list of positions on a chromosome, locate the index of
        its ancestral block within the numpy array of blocks

        Uses binary search to make things fast

        Parameters
        ----------
        blocks: npt.NDArray[np.uint32]
            The end positions of each ancestral block, sorted in ascending order
        positions: npt.NDArray[np.uint32]
            The position of each variant on a chromosome. These are the positions at
            which we'd like to query a sample's ancestry.

        Returns
        -------
        npt.NDArray[np.uint16]
            For each position in positions, return the index of its block in blocks
        """
        # use binary search to get the positions
        indices = np.searchsorted(blocks, positions, side="left")
        # if any indices exceed the length of the blocks, raise an error
        if np.any(indices >= len(blocks)):
            problem_position = positions[indices >= len(blocks)][0]
            raise ValueError(
                f"Position {problem_position} exceeds the range of the provided "
                "breakpoint blocks."
            )
        return indices

    def population_array(
        self,
        variants: np.array,
        samples: tuple[str] = None,
    ) -> npt.NDArray:
        """
        Output an array denoting the population labels of each variant for each sample

        Parameters
        ----------
        variants : np.array
            Variant-level meta information in a mixed np array of dtypes:
            CHROM (str) and POS (int)
        samples : tuple[str], optional
            A subset of samples to include in the output, ordered by their given order

        Returns
        -------
        npt.NDArray
            An array of shape: samples x variants x 2

            The array will have the same dtype as the population labels in the "pop"
            field of :py:attr:`~.Breakpoints.data`. Use :py:meth:`~.Breakpoints.encode`
            or :py:meth:`~.Breakpoints.recode` to change this.
        """
        if samples is None:
            data = self.data
        else:
            data = {samp: self.data[samp] for samp in samples}
        # initialize the return matrix
        dtype = HapBlock[0][1] if self.labels is None else np.uint8
        arr = np.empty((len(data), len(variants), 2), dtype=dtype)
        # iterate through the variants belonging to each chromosome
        gts_chroms = set(variants["chrom"])
        self.log.info(
            f"Obtaining ancestry for {len(data)} samples and {len(variants)} "
            f"variants in {len(gts_chroms)} chromosomes"
        )
        for chrom in gts_chroms:
            var_idxs = variants["chrom"] == chrom
            positions = variants["pos"][var_idxs]
            # obtain the population labels of each sample
            for samp_idx, samp_blocks in enumerate(data.values()):
                for strand_num in range(len(samp_blocks)):
                    blocks = samp_blocks[strand_num]
                    chrom_block = blocks[blocks["chrom"] == chrom]
                    if not len(chrom_block):
                        samp_id = tuple(data.keys())[samp_idx]
                        raise ValueError(
                            f"Chromosome {chrom} in the genotypes is absent in the "
                            f"breakpoints for sample {samp_id}_{strand_num+1}. Check "
                            "that your 'chr' prefixes match!"
                        )
                    # TODO: raise an exception if the end positions in chrom_block
                    # aren't sorted
                    # Now try to figure out the right population labels using binary
                    # search and then store them in the result matrix
                    try:
                        arr[samp_idx, var_idxs, strand_num] = chrom_block["pop"][
                            self._find_blocks(chrom_block["bp"], positions)
                        ]
                    except ValueError as e:
                        samp_id = tuple(data.keys())[samp_idx]
                        if str(e).startswith("Position "):
                            raise ValueError(
                                f"The breakpoints for chromosome {chrom} in sample"
                                f" {samp_id}_{strand_num+1} do not specify an ancestry"
                                " for one of the requested variants."
                            ) from e
                        else:
                            raise e
        return arr

    def write(self):
        """
        Write the breakpoints in this class to a file at :py:attr:`~.Breakpoints.fname`

        Examples
        --------
        To write to a file, you must first initialize a Breakpoints object and then
        fill out the names, data, and samples properties:

        >>> from haptools.data import Breakpoints, HapBlock
        >>> breakpoints = Breakpoints('simple.bp')
        >>> breakpoints.data = {
        >>>     'HG00096': [
        >>>         np.array([('YRI','chr1',10114,4.3),('CEU','chr1',10116,5.2)], dtype=HapBlock)
        >>>         np.array([('CEU','chr1',10114,4.3),('YRI','chr1',10116,5.2)], dtype=HapBlock)
        >>>     ], 'HG00097': [
        >>>         np.array([('YRI','chr1',10114,4.3),('CEU','chr2',10116,5.2)], dtype=HapBlock)
        >>>         np.array([('CEU','chr1',10114,4.3),('YRI','chr2',10116,5.2)], dtype=HapBlock)
        >>>     ]
        >>> }
        >>> breakpoints.write()
        """
        with self.hook_compressed(self.fname, mode="w") as bkpts:
            csv_writer = csv.writer(
                bkpts, delimiter="\t", dialect="unix", quoting=csv.QUOTE_NONE
            )
            for samp, blocks in self.data.items():
                for strand_num in range(len(blocks)):
                    bkpts.write(f"{samp}_{strand_num+1}\n")
                    for block in blocks[strand_num]:
                        csv_writer.writerow(block)
