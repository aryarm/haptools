.. _formats-genotypes:


Genotypes
=========

.. figure:: https://github.com/CAST-genomics/haptools/assets/23412689/6da88941-7520-4c19-beaa-27f540f6b047

	The time required to load various genotype file formats.

VCF/BCF
~~~~~~~

Genotype files must be specified as VCF or BCF files. They can be bgzip-compressed.

To be loaded properly, VCFs must follow the VCF specification. VCFs with duplicate variant IDs do not follow the specification; the IDs must be unique. Please validate your VCF using a tool like `gatk ValidateVariants <https://gatk.broadinstitute.org/hc/en-us/articles/360037057272-ValidateVariants>`_ before using haptools.

.. _formats-genotypesplink:

PLINK2 PGEN
~~~~~~~~~~~

There is also experimental support for `PLINK2 PGEN <https://github.com/chrchang/plink-ng/blob/master/pgen_spec/pgen_spec.pdf>`_ files (accomponied by PVAR and PSAM files) in some commands. These files can be loaded and created much more quickly than VCFs, so we highly recommend using them if you're working with large datasets. See the documentation for the :class:`GenotypesPLINK` class in :ref:`the API docs <api-data-genotypesplink>` for more information.

If you run out memory when using PGEN files, consider reading/writing variants from the file in chunks via the ``--chunk-size`` parameter.

Converting from VCF to PGEN
---------------------------
To convert a VCF containing only SNPs to PGEN, use the following command.

.. code-block:: bash

	plink2 --snps-only 'just-acgt' --vcf input.vcf --make-pgen --out output

To convert a VCF containing tandem repeats to PGEN, use this command, instead.

.. code-block:: bash

	plink2 --vcf-half-call m --make-pgen 'pvar-cols=vcfheader,qual,filter,info' --vcf input.vcf --make-pgen --out output

If you are seeing cryptic errors with haptools and your PGEN file, please validate it first:

.. code-block:: bash

	plink2 --pfile output --validate

Tandem repeats
~~~~~~~~~~~~~~
VCFs containing tandem repeats usually have a *type* indicating the tool from which they originated. We support whichever types are supported by `TRTools <https://trtools.readthedocs.io/en/stable/CALLERS.html>`_.

We do our best to infer the *type* of a TR VCF automatically. However, there will be some instances when it cannot be inferred.
Users of TRTools know to specify :code:`--vcftype` in that situation. However, most haptools commands do not yet support the :code:`--vcftype` parameter. Instead, you can modify the header of your VCF to trick haptools into being able to infer the *type*.

For example, if you're using HipSTR, you can add :code:`##command=hipstr...`. Refer to `this code in TRTools <https://trtools.readthedocs.io/en/stable/trtools.utils.tr_harmonizer.html#trtools.utils.tr_harmonizer.InferVCFType>`_ for more details.

Please note that all of this also applies to PVAR files created from TR VCFs.
