#!/usr/bin/env python

"""Tests for `vcf_consensus_builder` package."""

import re
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import SeqIO
from click.testing import CliRunner

from vcf_consensus_builder import cli

VCF_DEL = 'tests/data/test.vcf'
VCF_INS = 'tests/data/test2.vcf'
REF_FASTA = 'tests/data/ref.fa'
DEPTHS = 'tests/data/test-depths.tsv'
OUTPUT_FASTA = 'out.fa'


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()

    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert re.search(r'--help\s+Show this message and exit.', help_result.output, flags=re.DOTALL) is not None

    # Test replacing multi-char deletion and SNPs and setting sample name via command-line
    with TemporaryDirectory(prefix='vcf_consensus_builder', dir='/tmp') as tempdir:
        temppath = Path(tempdir)
        full_fasta_output = temppath / OUTPUT_FASTA
        sample_name = 'SAMPLE1'
        result = runner.invoke(cli.main,
                               ['-v', VCF_DEL,
                                '-d', DEPTHS,
                                '-r', REF_FASTA,
                                '-o', full_fasta_output,
                                '--sample-name', sample_name])
        assert result.exit_code == 0
        rec = SeqIO.read(full_fasta_output, 'fasta')
        assert str(rec.seq) == 'NACCGTANACAATAN--', 'There must be a deletion of 3 characters in the middle of the seq'
        assert rec.id == sample_name

    # Test changing no coverage threshold
    with TemporaryDirectory(prefix='vcf_consensus_builder', dir='/tmp') as tempdir:
        temppath = Path(tempdir)
        full_fasta_output = temppath / OUTPUT_FASTA
        sample_name = 'SAMPLE1'
        result = runner.invoke(cli.main,
                               ['-v', VCF_DEL,
                                '-d', DEPTHS,
                                '-r', REF_FASTA,
                                '-o', full_fasta_output,
                                '--sample-name', sample_name,
                                '--no-coverage', 4])
        assert result.exit_code == 0
        rec = SeqIO.read(full_fasta_output, 'fasta')
        assert str(rec.seq) == '-ACCGTA-ACAATA---', 'Positions below 5X coverage must be replaced with "-"'
        assert rec.id == sample_name

    # Test replacing low and no coverage characters with other characters than default N and - respectively
    with TemporaryDirectory(prefix='vcf_consensus_builder', dir='/tmp') as tempdir:
        temppath = Path(tempdir)
        full_fasta_output = temppath / OUTPUT_FASTA
        sample_name = 'SAMPLE1'
        result = runner.invoke(cli.main,
                               ['-v', VCF_DEL,
                                '-d', DEPTHS,
                                '-r', REF_FASTA,
                                '-o', full_fasta_output,
                                '--sample-name', sample_name,
                                '--no-cov-char', '=',
                                '--low-cov-char', '@'])
        assert result.exit_code == 0
        rec = SeqIO.read(full_fasta_output, 'fasta')
        assert str(rec.seq) == '@ACCGTA@ACAATA@==', \
            'No coverage positions must be replaced with "=". Low coverage (<5X) positions must be replaced with "@".'
        assert rec.id == sample_name

    # Test replacing multi-char insertion and SNPs
    with TemporaryDirectory(prefix='vcf_consensus_builder', dir='/tmp') as tempdir:
        temppath = Path(tempdir)
        full_fasta_output = temppath / OUTPUT_FASTA
        result = runner.invoke(cli.main,
                               ['-v', VCF_INS,
                                '-d', DEPTHS,
                                '-r', REF_FASTA,
                                '-o', full_fasta_output])
        assert result.exit_code == 0
        rec = SeqIO.read(full_fasta_output, 'fasta')
        assert str(rec.seq) == 'NACCGTATTTGTCNACAATAN--', 'There must be an insertion of "TTT" in the middle of the seq'
        assert rec.id == 'sample1', \
            'FASTA ID must be the first sample name in the VCF if not explicitly specified as a command-line arg'
