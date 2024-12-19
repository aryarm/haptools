#!/usr/bin/env python

# To test: ./clumpSTR.py --summstats-snps tests/eur_gwas_pvalue_chr19.LDL.glm.linear --clump-snp-field ID --clump-field p-value --clump-chrom-field CHROM --clump-pos-field position --clump-p1 0.2 --out test.clump
import math
from logging import Logger, getLogger

import numpy as np

from .data import Genotypes, GenotypesVCF, GenotypesPLINK, GenotypesTR, GenotypesPLINKTR


class Variant:
    def __init__(self, varid, chrom, pos, pval, vartype):
        self.varid = varid
        self.chrom = chrom
        self.pos = pos
        self.pval = pval
        self.vartype = vartype

    def __str__(self):
        return "%s %s %s %s %s" % (
            self.varid,
            self.chrom,
            self.pos,
            self.pval,
            self.vartype,
        )


class SummaryStats:
    """
    Load and process summary statistics

    Attributes
    ----------
    summstats: list[Variant]
        list of Variant objects that represent all summary
        statistics in the file.
    log: Logger
        A logging instance for recording debug statements.

    Examples
    --------
    Loading a summary stats file, grabbing an index variant, and its
    candidate variants to calculate LD.

    >>> summstats = SummaryStats(log)
    >>> summstats.Load(summstats_strs, vartype="SNP", pthresh=clump_p2,
            id_field=clump_id_field, p_field=clump_field,
            chrom_field=clump_chrom_field, pos_field=clump_pos_field)
    >>> indexvar = summstats.GetNextIndexVariant(clump_p1)
    >>> candidates = summstats.QueryWindow(indexvar, clump_kb)
    """

    def __init__(self, log: Logger = None):
        self.summstats = []
        self.log = log or getLogger(self.__class__.__name__)

    def Load(
        self,
        statsfile,
        vartype="SNP",
        pthresh=1.0,
        id_field="SNP",
        p_field="P",
        chrom_field="CHR",
        pos_field="POS",
    ):
        """
        Load summary statistics
        Ignore variants with pval < pthresh
        Not yet implemented
        """
        summstats = []  # List of Variants

        # First, parse header line to get col. numbers
        f = open(statsfile, "r")
        header = f.readline()
        if header.startswith("#"):
            header = header[1:]
        header_items = [item.strip() for item in header.split()]
        try:
            snp_col = header_items.index(id_field)
        except ValueError:
            raise ValueError("Could not find %s in header" % id_field)
        try:
            p_col = header_items.index(p_field)
        except ValueError:
            raise ValueError("Could not find %s in header" % p_field)
        try:
            chrom_col = header_items.index(chrom_field)
        except ValueError:
            raise ValueError("Could not find %s in header" % chrom_field)
        try:
            pos_col = header_items.index(pos_field)
        except ValueError:
            raise ValueError("Could not find %s in header" % pos_field)

        # Now, load in stats. Skip things with pval>pthresh
        line = f.readline()
        while line.strip() != "":
            items = [item.strip() for item in line.strip().split()]
            if float(items[p_col]) > pthresh:
                line = f.readline()
                continue
            summstats.append(
                Variant(
                    items[snp_col],
                    items[chrom_col],
                    int(items[pos_col]),
                    float(items[p_col]),
                    vartype,
                )
            )
            line = f.readline()
        f.close()
        self.summstats.extend(summstats)

    def GetNextIndexVariant(self, index_pval_thresh):
        """
        Get the next index variant, which is the
        variant with the best p-value
        If no more variants below the clump-p1 threshold,
        return None
        Not yet implemented
        """
        best_var = None
        best_var_p = 1.0
        for variant in self.summstats:
            if variant.pval < best_var_p and variant.pval < index_pval_thresh:
                best_var = variant
                best_var_p = variant.pval
        return best_var

    def QueryWindow(self, indexvar, window_kb):
        """
        Find all candidate variants in the specified
        window around the index variant
        Not yet implemented
        """
        # First get stats on the indexvariant
        chrom = indexvar.chrom
        pos = indexvar.pos

        # Find candidates in the window
        candidates = []
        for variant in self.summstats:
            if variant.chrom == chrom and abs(variant.pos - pos) / 1000 < window_kb:
                candidates.append(variant)

        return candidates

    def RemoveClump(self, clumpvars):
        """
        Remove the variants from a clump
        from further consideration
        """
        keepvars = []
        for variant in self.summstats:
            if variant not in clumpvars:
                keepvars.append(variant)
        self.summstats = keepvars


def _SortSamples(samples):
    """
    Sort samples along with their indices.
    """
    # Create indices to track for each sample
    inds = np.arange(len(samples))

    # Sort indices and samples together
    sorted_data = [[sample, ind] for sample, ind in sorted(zip(samples, inds))]

    # Grab samples and inds separately
    sorted_samples = [sample for sample, _ in sorted_data]
    sorted_inds = [ind for _, ind in sorted_data]

    return sorted_samples, sorted_inds


def GetOverlappingSamples(snpgts, strgts):
    """
    Get indices of overlapping samples for snps and strs

    Parameters
    ----------
    snpgts: GenotypesVCF
        SNP Genotypes object
    strgts: GenotypesTR
        STR Genotypes object

    Returns
    -------
    snp_samples: list(int)
        Indices of overlapping samples for snps
    str_samples: list(int)
        Indices of overlapping samples for strs
    """
    # Sort samples and respective indices in sample array together
    snp_match_inds = []
    str_match_inds = []
    snp_samples, snp_inds = _SortSamples(snpgts.samples)
    str_samples, str_inds = _SortSamples(strgts.samples)

    # Since both lists are sorted we iterate over each list using counters
    # to determine where potential matching samples are
    snp_counter = 0
    str_counter = 0
    while snp_counter < len(snp_inds) and str_counter < len(str_inds):
        if str_samples[str_counter] < snp_samples[snp_counter]:
            str_counter += 1
        elif str_samples[str_counter] == snp_samples[snp_counter]:
            snp_match_inds.append(snp_inds[snp_counter])
            str_match_inds.append(str_inds[str_counter])
            snp_counter += 1
            str_counter += 1
        else:
            snp_counter += 1

    return snp_match_inds, str_match_inds


def LoadVariant(var, gts, log):
    """
    Extract vector of genotypes for this variant
    """
    # Grab variant from snps or strs depending on variant type
    var_ind = (gts.variants["pos"] == int(var.pos)) & (
        gts.variants["chrom"] == var.chrom
    )
    data = gts.data[:, var_ind, :]
    data = data.reshape(data.shape[0], data.shape[2])
    return data


def _CalcChiSQ(f00, f01, f10, f11, gt_counts, n):
    """
    Calculate Chi-squared test stat for given freqs.
    """
    chisq_exp = np.zeros((3, 3))
    root_exp = np.zeros((3, 3))

    # calculate expected values for a given root
    root_exp[0, 0] = n * f00**2
    root_exp[0, 1] = 2 * n * f00 * f01
    root_exp[0, 2] = n * f01**2
    root_exp[1, 0] = 2 * n * f00 * f10
    root_exp[1, 1] = 2 * n * f01 * f10 + 2 * n * f00 * f11
    root_exp[1, 2] = 2 * n * f01 * f11
    root_exp[2, 0] = n * f10**2
    root_exp[2, 1] = 2 * n * f10 * f11
    root_exp[2, 2] = n * f11**2
    for i in range(3):
        for j in range(3):
            if root_exp[i, j] > 0.0:
                chisq_exp = (gt_counts[i, j] - root_exp[i, j]) ** 2 / root_exp[i, j]

    return np.sum(chisq_exp)


def _CalcLDStats(f00, p, q, gt_counts, n):
    """
    Given frequency of gt 0|0 (f11) and major and minor allele freqs p and q calculate stats.
    """
    f01 = p - f00
    f10 = q - f00
    f11 = 1 - (f00 + f01 + f10)
    D = (f00 * f11) - (f01 * f10)
    if D >= 0.0:
        Dmax = min(p * (1.0 - q), q * (1.0 - p))
    else:
        Dmax = min(p * q, (1 - p) * (1 - q))
    Dprime = D / Dmax
    r_squared = (D**2) / (p * (1 - p) * q * (1 - q))

    return (
        round(Dprime, 6),
        round(r_squared, 6),
        _CalcChiSQ(f00, f01, f10, f11, gt_counts, n),
    )


def _CalcBestRoot(real_roots, minhap, maxhap, p, q, gt_counts, n):
    """
    Given a list of real roots (max possible 3) calculate the best root.
    The best root is the one with the lowest chisq test statistic value.
    """
    # determine the best root by grabbing the one with the lowest chisq
    best_Dprime = 0
    best_rsquared = 0
    best_chisq = np.inf
    for root_freq00 in real_roots:
        # calculate LD stats given root freq is within bounds
        if root_freq00 >= minhap - 0.00001 and root_freq00 <= maxhap + 0.00001:
            Dprime, r_squared, chisq = _CalcLDStats(root_freq00, p, q, gt_counts, n)
            if chisq < best_chisq:
                best_Dprime = Dprime
                best_rsquared = r_squared
                best_chisq = chisq

    return best_Dprime, best_rsquared


def ComputeExactLD(candidate_gt, index_gt, log):
    """
    Compute exact solution of haplotype frequencies to calculate r squared value.
    NOTE currently this approach only works for biallelic variants since having more variants
    causes the equation we're solving for to be a cubic but instead to the degree of n where
    n is the total number of alelles which also invalidates STRs.

    Parameters
    ----------
    candidate_gt: np.array
        array of size (genotypes,) where genotypes is the number of samples
    index_gt: np.array
        array of size (genotypes,) where genotypes is the number of samples
    log: Logger
        A logging instance for recording debug statements.

    Returns
    -------
    r_squared: float
        R squared value inferred from ML solution.
    """
    # load in 3x3 array where axes are genotypes (0,1,2) for each variant
    # y-axis = candidate gt, x-axis = index_gt
    gt_counts = np.zeros((3, 3))
    for gt1 in range(3):
        # subset candidate gts to genotype gt
        subset_gt1 = candidate_gt == gt1
        for gt2 in range(3):
            subset_index_gt = index_gt[subset_gt1]
            gt_counts[gt1, gt2] = np.sum(subset_index_gt == gt2)

    n = np.sum(gt_counts)
    p = (2.0 * np.sum(gt_counts[0, :]) + np.sum(gt_counts[1, :])) / (2.0 * n)
    q = (2.0 * np.sum(gt_counts[:, 0]) + np.sum(gt_counts[:, 1])) / (2.0 * n)

    num_alt = 2.0 * gt_counts[0, 0] + gt_counts[0, 1] + gt_counts[1, 0]
    a = 4.0 * n
    b = 2.0 * n * (1.0 - 2.0 * p - 2.0 * q) - 2.0 * num_alt - gt_counts[1, 1]
    c = (
        -num_alt * (1.0 - 2.0 * p - 2.0 * q)
        - gt_counts[1, 1] * (1.0 - p - q)
        + 2.0 * n * p * q
    )
    d = -num_alt * p * q

    minhap = num_alt / (2.0 * float(n))
    maxhap = (num_alt + gt_counts[1, 1]) / (2.0 * float(n))

    xN = -b / (3.0 * a)
    d2 = (math.pow(b, 2) - 3.0 * a * c) / (9 * math.pow(a, 2))
    yN = a * math.pow(xN, 3) + b * math.pow(xN, 2) + c * xN + d
    yN2 = math.pow(yN, 2)
    h2 = 4 * math.pow(a, 2) * math.pow(d2, 3)

    # store all real roots to cubic to iterate over and determine which is best
    real_roots = []

    # three possible scenarios of solutions
    if yN2 > h2:
        # calculate real root alpha
        number1 = 0.0
        number2 = 0.0
        if (1.0 / (2.0 * a) * (-yN + math.pow((yN2 - h2), 0.5))) < 0:
            number1 = -math.pow(
                -(1.0 / (2.0 * a) * (-yN + math.pow((yN2 - h2), 0.5))), 1.0 / 3.0
            )
        else:
            number1 = math.pow(
                (1.0 / (2.0 * a) * (-yN + math.pow((yN2 - h2), 0.5))), 1.0 / 3.0
            )

        if (1.0 / (2.0 * a) * (-yN - math.pow((yN2 - h2), 0.5))) < 0:
            number2 = -math.pow(
                -(1.0 / (2.0 * a) * (-yN - math.pow((yN2 - h2), 0.5))), 1.0 / 3.0
            )
        else:
            number2 = math.pow(
                (1.0 / (2.0 * a) * (-yN - math.pow((yN2 - h2), 0.5))), 1.0 / 3.0
            )

        # singular real root
        alpha = xN + number1 + number2

        # store real root alpha
        real_roots = [alpha]

    elif yN2 == h2:
        # Calculate three real roots alpha beta and gamma
        delta = math.pow((yN / 2.0 * a), (1.0 / 3.0))
        alpha = xN + delta
        beta = xN + delta
        gamma = xN - 2.0 * delta

        # store all real roots
        real_roots = [alpha, beta, gamma]

    elif yN2 < h2:
        # calculate 3 real roots alpha beta and gamma
        h = math.pow(h2, 0.5)
        theta = (math.acos(-yN / h)) / 3.0
        delta = math.pow(d2, 0.5)
        alpha = xN + 2.0 * delta * math.cos(theta)
        beta = xN + 2.0 * delta * math.cos(2.0 * math.pi / 3.0 + theta)
        gamma = xN + 2.0 * delta * math.cos(4.0 * math.pi / 3.0 + theta)

        # store all real roots
        real_roots = [alpha, beta, gamma]

    else:
        raise Exception(f"Can't calculate r squared from given values {yN2} and {h2}")

    # Solve for best roots
    best_Dprime, best_rsquared = _CalcBestRoot(
        real_roots, minhap, maxhap, p, q, gt_counts, n
    )
    return best_Dprime, best_rsquared


def _FilterGts(candidate_gt, index_gt, log):
    """
    Filter invalid values from gts, 254 and 255 since -2 and -1 encode for these
    once converted to uint8 and sum alleles together.
    """
    # Check Alleles for invalid values and remove samples from both sets of alleles
    miss = np.iinfo(np.uint8).max - 1
    valid_gts = np.all(candidate_gt < miss, axis=1) & np.all(index_gt < miss, axis=1)
    candidate_gt = candidate_gt[valid_gts, :]
    index_gt = index_gt[valid_gts, :]
    log.debug(f"Valid Genotype Indices: {valid_gts}")
    log.debug(f"Candidate GTs: {candidate_gt}")
    log.debug(f"Index GTs: {index_gt}")

    # Sum alleles
    candidate_gt = np.sum(candidate_gt, axis=1).flatten()
    index_gt = np.sum(index_gt, axis=1).flatten()
    return candidate_gt, index_gt


def ComputeLD(candidate_gt, index_gt, LD_type, log):
    """
    Compute the LD between two variants
    """
    # Filter invalid gt values
    candidate_gt, index_gt = _FilterGts(candidate_gt, index_gt, log)
    if not (np.size(candidate_gt) and np.size(index_gt)):
        return None, 0

    # Check if all values in either array are the same
    if np.unique(candidate_gt).shape[0] == 1 or np.unique(index_gt).shape[0] == 1:
        log.debug("GTs between one of the variants are constant across all samples.")
        return None, np.nan

    # Compute and Maximum likelihood solution or Pearson r2
    if LD_type == "Exact":
        return ComputeExactLD(candidate_gt, index_gt, log)
    elif LD_type == "Pearson":
        return None, np.corrcoef(index_gt, candidate_gt)[0, 1] ** 2


def WriteClump(indexvar, clumped_vars, outf):
    """
    Write a clump to the output file
    Not yet implemented
    """
    outf.write(
        "\t".join(
            [
                indexvar.varid,
                indexvar.chrom,
                str(indexvar.pos),
                str(indexvar.pval),
                indexvar.vartype,
                ",".join([str(item) for item in clumped_vars]),
            ]
        )
        + "\n"
    )


def clumpstr(
    summstats_snps,
    summstats_strs,
    gts_snps,
    gts_strs,
    clump_p1,
    clump_p2,
    clump_id_field,
    clump_field,
    clump_chrom_field,
    clump_pos_field,
    clump_kb,
    clump_r2,
    LD_type,
    out,
    log,
):
    ###### User checks ##########
    # if summstats_snps, also need gts_snps
    log.debug(f"Validating SNP files {summstats_snps} or {gts_snps}")
    if summstats_snps or gts_snps:
        try:
            assert gts_snps
            assert summstats_snps
        except:
            raise Exception(
                f"One of summstats-snps {summstats_snps} and gts-snps {gts_snps} is "
                "not present. Please ensure both have been inputted correctly."
            )

    # if summstats_strs, also need gts_strs
    log.debug(f"Validating STR files {summstats_strs} or {gts_strs}")
    if summstats_strs or gts_strs:
        try:
            assert gts_strs
            assert summstats_strs
        except:
            raise Exception(
                f"One of summstats-strs {summstats_strs} and gts-strs {gts_strs} is "
                "not present. Please ensure both have been inputted correctly."
            )

    if summstats_strs and LD_type == "Exact":
        raise Exception(
            "The exact method of computing LD can only be used with biallelic loci. "
            + "STRs are not compatible with the exact LD compute method. "
        )

    ###### Load summary stats ##########
    summstats = SummaryStats(log)
    if summstats_snps is not None:
        summstats.Load(
            summstats_snps,
            vartype="SNP",
            pthresh=clump_p2,
            id_field=clump_id_field,
            p_field=clump_field,
            chrom_field=clump_chrom_field,
            pos_field=clump_pos_field,
        )
    if summstats_strs is not None:
        summstats.Load(
            summstats_strs,
            vartype="STR",
            pthresh=clump_p2,
            id_field=clump_id_field,
            p_field=clump_field,
            chrom_field=clump_chrom_field,
            pos_field=clump_pos_field,
        )

    ###### Set up genotypes ##########
    snpgts = None
    strgts = None
    gts = None
    if gts_snps:
        log.debug("Loading SNP Genotypes.")
        if str(gts_snps).endswith("pgen"):
            snpgts = GenotypesPLINK.load(gts_snps)
        else:
            snpgts = GenotypesVCF.load(gts_snps)
    if gts_strs:
        log.debug("Loading STR Genotypes.")
        if str(gts_strs).endswith("pgen"):
            strgts = GenotypesPLINKTR.load(gts_strs)
        else:
            strgts = GenotypesTR.load(gts_strs)

    if gts_snps and gts_strs:
        log.debug("Calculating set of overlapping samples between STRs and SNPs.")
        # Grab all shared samples between snp list and str list
        # NOTE samples are returned such that resulting data is matching
        snp_samples, str_samples = GetOverlappingSamples(snpgts, strgts)
        log.debug(f"Shared samples: {snp_samples}")

        # NOTE snpgts has data, variants, and samples where data is alleles (samples x variants x alleles)
        #      variants has id, pos, chrom, ref, alt for snps
        #      samples is list of samples corresponding to x-axis of samples
        snpgts.data = snpgts.data[snp_samples, :, :]
        snpgts.samples = tuple(np.array(snpgts.samples)[snp_samples])
        strgts.data = strgts.data[str_samples, :, :]
        strgts.samples = tuple(np.array(strgts.samples)[str_samples])

        # Merge STR and SNP GTs
        # NOTE if Genotypes is not used and GenotypesVCF is instead it will error because
        #      GenotypesVCF requires alleles to be present and Genotypes does not
        gts = Genotypes.merge_variants((snpgts, strgts), fname=None)
    elif gts_snps:
        gts = snpgts
    elif gts_strs:
        gts = strgts
    else:
        raise Exception("Unable to load valid genotype data.")

    ###### Setup output file ##########
    outf = open(out, "w")
    outf.write("\t".join(["ID", "CHROM", "POS", "P", "VARTYPE", "CLUMPVARS"]) + "\n")

    ###### Perform clumping ##########
    indexvar = summstats.GetNextIndexVariant(clump_p1)
    while indexvar is not None:
        # Load indexvar gts
        indexvar_gt = LoadVariant(indexvar, gts, log)
        # Collect candidate variants within range of index variant
        candidates = summstats.QueryWindow(indexvar, clump_kb)

        log.debug(f"Current index variant: {indexvar}")

        # calculate LD between candidate vars and index var
        clumpvars = []
        for c in candidates:
            # load candidate variant c genotypes
            candidate_gt = LoadVariant(c, gts, log)
            Dprime, r2 = ComputeLD(candidate_gt, indexvar_gt, LD_type, log)
            # If using pearson Dprime is not calculated
            log.debug(
                f"D' and r2 between {indexvar} with {c}\nD' = {Dprime}, r^2 = {r2}"
            )
            if r2 > clump_r2:
                clumpvars.append(c)
        WriteClump(indexvar, clumpvars, outf)
        summstats.RemoveClump(clumpvars + [indexvar])
        indexvar = summstats.GetNextIndexVariant(clump_p1)
    outf.close()
