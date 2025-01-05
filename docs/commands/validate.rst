.. _commands-validate:


validate
========

Validate the formatting of a :doc:`.hap file </formats/haplotypes>`. Output warnings/errors explaining how the formatting of your ``.hap`` file may be improved.

If a :ref:`.pvar file <formats-genotypesplink>` file is provided, the SNPs and TRs present in the ``.hap`` file will be checked for existence in the ``.pvar`` file.

.. note::

  This command will not check that your ``.hap`` file is properly sorted. It only checks formatting.

Usage
~~~~~
.. code-block:: bash

  haptools validate \
  --sort \
  --genotypes PATH \
  --verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
  HAPFILE

Examples
~~~~~~~~
.. code-block:: bash

  haptools validate tests/data/validate/basic.hap

Outputs a message specifying the amount of errors and warnings.

.. code-block::

  [    INFO] Completed .hap file validation with 0 errors and 0 warnings.

All warnings and errors will be logged if there are any.

.. code-block:: bash

  haptools validate tests/data/validate/no_version.hap

.. code-block::

  [ WARNING] No version declaration found. Assuming to use the latest version.
  [    INFO] Completed .hap file validation with 0 errors and 1 warnings.
  Error: Found several warnings and / or errors in the .hap file

All ``.hap`` files must be sorted before they can be validated, so we try our best to sort your ``.hap`` file internally before performing any validation checks.
If your ``.hap`` file is already sorted, you should use the ``--sorted`` parameter. It will speed things up a bit by skipping the sorting step. If your ``.hap`` file is indexed, it will be assumed to be sorted regardless.

.. code-block:: bash

  haptools validate --sorted tests/data/simple.hap

As mentioned before, one can use the ``--genotypes`` flag to provide a ``.pvar`` file with which to compare the existence of variant IDs.
The following will check if all of the variant IDs in the ``.hap`` file appear in the ``.pvar`` file.

.. code-block:: bash

  haptools validate --genotypes tests/data/simple.pvar tests/data/simple.hap

.. note::

  We accept a PVAR file instead of a VCF in order to avoid reading lots of information
  which is not relevant to the validation process. However, any VCF subsetted to just
  its first 8 fields is a valid PVAR file. So you can easily create a PVAR file from a
  VCF using ``cut -f -8`` or ``plink2 --make-just-pvar``.

Detailed Usage
~~~~~~~~~~~~~~

.. click:: haptools.__main__:main
  :prog: haptools
  :show-nested:
  :commands: validate
