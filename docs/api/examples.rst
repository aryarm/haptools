.. _api-examples:


examples
========

Converting a ``.blocks.det`` file into a ``.hap`` file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can use the :ref:`data API <api-data>` to convert `a PLINK 1.9 .blocks.det file <https://www.cog-genomics.org/plink/1.9/formats#blocks>`_ into a ``.hap`` file.

As an example, let's say we would like to convert `the following simple.blocks.det file <https://github.com/cast-genomics/haptools/blob/main/tests/data/simple.blocks.det>`_.

.. include:: ../../tests/data/simple.blocks.det
  :literal:

.. code-block:: python

    from haptools import data

    # load the genotypes file
    # you can use either a VCF or PGEN file
    gt = data.GenotypesVCF.load("tests/data/simple.vcf.gz")
    gt = data.GenotypesPLINK.load("tests/data/simple.pgen")

    # load the haplotypes
    hp = data.Haplotypes("output.hap")
    hp.data = {}

    # iterate through lines of the .blocks.det file
    with open("tests/data/simple.blocks.det") as blocks_file:
        for idx, line in enumerate(blocks_file.read().splitlines()[1:]):
            # initialize variables and parse line from the blocks file
            hap_id = f"H{idx}"
            chrom, bp1, bp2, kb, nsnps, snps = line.strip().split()

            # create a haplotype line in the .hap file
            hp.data[hap_id] = data.Haplotype(
                chrom=chrom, start=int(bp1), end=int(bp2), id=hap_id
            )

            # fetch alleles from the genotypes file
            snp_gts = gt.subset(variants=tuple(snps.split("|")))

            # create variant lines for each haplotype
            # Note that the .blocks.det file doesn't specify an allele, so
            # we simply choose the first allele (ie the REF allele) for this example
            hp.data[hap_id].variants = tuple(
                data.Variant(
                    start=v["pos"],
                    end=v["pos"] + len(v["alleles"][0]),
                    id=v["id"],
                    allele=v["alleles"][0],
                )
                for v in snp_gts.variants
            )

    hp.write()

.. _api-examples-snps2hap:

Converting a ``.snplist`` file into a ``.hap`` file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
How would you convert a :doc:`.snplist file </formats/snplist>` into a ``.hap`` file suitable for use by ``simphenotype``?

The basic idea is to encode each SNP as a haplotype containing only a single allele. For example, let's say your ``.snplist`` file has two SNPs like this.

.. include:: ../../tests/data/apoe.snplist
  :literal:

Then your ``.hap`` file might look something like this.

.. include:: ../../tests/data/apoe.hap
  :literal:

You can use the :ref:`data API <api-data>` and the :ref:`simphenotype API <api-haptools-sim_phenotype>` to create such a file.

.. code-block:: python

    from haptools import data
    from haptools.sim_phenotype import Haplotype

    variants = {}
    # load variants from the snplist file
    with open("tests/data/apoe.snplist") as snplist_file:
        for line in snplist_file.readlines():
            # parse variant ID and beta from file
            ID, beta = line.split("\t")
            variants[ID] = float(beta)

    # load the genotypes file
    gt = data.GenotypesVCF("tests/data/apoe.vcf.gz")
    gt.read(variants=variants.keys())

    # initialize an empty haplotype file
    hp = data.Haplotypes("output.hap", haplotype=Haplotype)
    hp.data = {}

    for variant in gt.variants:
        ID, chrom, pos, alleles = variant[["id", "chrom", "pos", "alleles"]]        
        # we arbitrarily choose to use the ALT allele but alleles[0] will give you REF
        allele = alleles[1]
        end = pos + len(allele)

        # create a haplotype line in the .hap file
        hp.data[ID] = Haplotype(chrom=chrom, start=pos, end=end, id=ID, beta=variants[ID])

        # create a variant line for each haplotype
        hp.data[ID].variants = (data.Variant(start=pos, end=end, id=ID, allele=allele),)

    hp.write()

.. _api-examples-bp2anc:

Converting a ``.bp`` file into a ``.hanc`` per-site ancestry file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can obtain the ancestry of a list of variants directly from a ``.bp`` file using the :ref:`data API <api-data-bp2anc>`.

**Input:**

* Breakpoints in a :ref:`.bp file <formats-breakpoints>`
* A list of variants in a :ref:`PLINK2 PVAR file <formats-genotypesplink>`

**Output:**

* An ``.hanc`` per-site ancestry file as described in `the admix-simu documentation <https://github.com/williamslab/admix-simu/tree/master?tab=readme-ov-file#per-site-ancestry-values>`_:

.. include:: ../../tests/data/simple.hanc
  :literal:

.. code-block:: python

    import numpy as np
    from pathlib import Path
    from haptools import data

    output = Path("output.hanc")

    # load breakpoints from the bp file and encode each population label as an int
    breakpoints = data.Breakpoints.load("tests/data/simple.bp")
    breakpoints.encode()
    print(breakpoints.labels)

    # load the SNPs array from a PVAR file
    snps = data.GenotypesPLINK("tests/data/simple.pgen")
    snps.read_variants()
    snps = snps.variants[["chrom", "pos"]]

    # create array of per-site ancestry values
    arr = breakpoints.population_array(variants=snps)
    # reshape from n x p x 2 to n*2 x p
    # so rows are haplotypes and columns are variants
    arr = arr.transpose((0, 2, 1)).reshape(-1, arr.shape[1])

    # write to haplotype ancestry file
    np.savetxt(output, arr, fmt="%i", delimiter="")
