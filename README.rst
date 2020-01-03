=====================
VCF Consensus Builder
=====================


.. image:: https://img.shields.io/pypi/v/vcf_consensus_builder.svg
        :target: https://pypi.python.org/pypi/vcf_consensus_builder

.. image:: https://img.shields.io/travis/peterk87/vcf_consensus_builder.svg
        :target: https://travis-ci.org/peterk87/vcf_consensus_builder

.. image:: https://readthedocs.org/projects/vcf-consensus-builder/badge/?version=latest
        :target: https://vcf-consensus-builder.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Build a consensus sequence from a VCF and reference sequence masking low and no coverage positions.

You *could* use ``bcftools consensus`` but then you would need to apply the low and no coverage position masking *after* ``bcftools`` has generated the consensus, which may be tricky.


* Free software: MIT license
* Documentation: https://vcf-consensus-builder.readthedocs.io.


Features
--------

* Masks low and no coverage positions in the reference (default: 0X and <5X) with ``N`` and ``-`` by default
* No need to ``bgzip`` the VCF file or index it like ``bcftools consensus`` requires.
*

Usage
-----


Install
~~~~~~~

Install with ``pip`` from PyPI with

.. code-block::

    pip install vcf_consensus_builder


Show Help
~~~~~~~~~

Help message:

.. code-block::

    $ vcf_consensus_builder --help
    Usage: vcf_consensus_builder [OPTIONS]

      Build a consensus sequence from a VCF and ref sequence masking low and no
      coverage positions.

    Options:
      -v, --vcf-file PATH      VCF file path (v4)  [required]
      -d, --depths-file PATH   samtools depth output file (no headers)  [required]
      -r, --ref-fasta PATH     Reference sequence FASTA file (single sequence
                               entry only!)  [required]
      -o, --output-fasta TEXT  Output consensus sequence FASTA file path (default
                               write to stdout)
      --low-coverage INTEGER   Low coverage threshold; replace positions with less
                               than this depth with "N" by default
      --no-coverage INTEGER    No coverage threshold; replace positions with less
                               than or equal this depth with "-" by default
      --low-cov-char TEXT      Low coverage character ("N" by default)
      --no-cov-char TEXT       No coverage character ("-" by default)
      --sample-name TEXT       Optional sample name for output fasta header ID
      --help                   Show this message and exit.


Basic usage
~~~~~~~~~~~

Run on the test data including in the repo

.. code-block::

    # Clone this repo and enter it
    $ git clone https://github.com/peterk87/vcf_consensus_builder.git --depth=1
    $ cd vcf_consensus_builder/
    # run vcf_consensus_builder on test data
    $ vcf_consensus_builder -v tests/data/test.vcf \
        -d tests/data/test-depths.tsv \
        -r tests/data/ref.fa
    # produces the following to stdout
    >sample1 ref="ref ref"
    NACCGTANACAATAN--


Masking of no and low coverage positions in reference sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``vcf_consensus_builder`` first masks no and low coverage positions in the reference sequence file and then applies the ``ALT`` variants in the VCF.

**NOTE:** ``vcf_consensus_builder`` does not perform any VCF variant filtering. It is assumed that the VCF input file contains only variants you wish to see in your consensus sequence. Please use ``bcftools filter`` with appropriate filtering/exclusion expressions to get the variants you wish to see represented in your consensus (see https://samtools.github.io/bcftools/howtos/filtering.html for more info about how to filter your VCF file)

Given this reference sequence

.. code-block::

    >ref
    NGCCAAGTCTNCGACATN-

And this  ``samtools depth`` output

.. code-block::

    sample1	ref	1	4
    sample1	ref	2	9
    sample1	ref	3	9
    sample1	ref	4	9
    sample1	ref	5	9
    sample1	ref	6	9
    sample1	ref	7	10
    sample1	ref	8	10
    sample1	ref	9	10
    sample1	ref	10	10
    sample1	ref	11	3
    sample1	ref	12	9
    sample1	ref	13	9
    sample1	ref	14	9
    sample1	ref	15	9
    sample1	ref	16	9
    sample1	ref	17	5
    sample1	ref	18	4
    sample1	ref	19	0
    sample1	ref	20	0

The low (below 5X) and no (0X) coverage positions in the reference sequence will be replaced with ``N`` and ``-``, respectively.

The masked reference sequence will be:

.. code-block::

    >ref
    NGCCAAGTCTNCGACATN-

This masked sequence will be used for generating the consensus sequence.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
