"""Console script for vcf_consensus_builder."""
import sys
import click

from vcf_consensus_builder.log import init_console_logger
from .vcf_consensus_builder import consensus


@click.command()
@click.option('-v', '--vcf-file',
              type=click.Path(exists=True),
              required=True,
              help='VCF file path (v4)')
@click.option('-d', '--depths-file',
              type=click.Path(exists=True),
              required=True,
              help='samtools depth output file (no headers)')
@click.option('-r', '--ref-fasta',
              type=click.Path(exists=True),
              required=True,
              help='Reference sequence FASTA file (single sequence entry only!)')
@click.option('-o', '--output-fasta',
              default=sys.stdout,
              help='Output consensus sequence FASTA file path (default write to stdout)')
@click.option('--low-coverage', type=int, default=5,
              help='Low coverage threshold; replace positions with less than this depth with "N" by default')
@click.option('--no-coverage', type=int, default=0,
              help='No coverage threshold; replace positions with less than or equal this depth with "-" by default')
@click.option('--low-cov-char', type=str, default='N',
              help='Low coverage character ("N" by default)')
@click.option('--no-cov-char', type=str, default='-',
              help='No coverage character ("-" by default)')
@click.option('--sample-name', type=str, required=False,
              help='Optional sample name for output fasta header ID')
@click.option('-V', '--verbose', default=0, count=True, help='Verbosity of logging')
def main(vcf_file,
         depths_file,
         ref_fasta,
         output_fasta,
         low_coverage,
         no_coverage,
         low_cov_char,
         no_cov_char,
         sample_name='SAMPLE',
         verbose=0):
    """Build a consensus sequence from a VCF and ref sequence masking low and no coverage positions."""
    init_console_logger(verbose)
    consensus(vcf_file=vcf_file,
              depths_file=depths_file,
              ref_fasta=ref_fasta,
              output_fasta=output_fasta,
              low_coverage=low_coverage,
              no_coverage=no_coverage,
              low_cov_char=low_cov_char,
              no_cov_char=no_cov_char,
              sample_name=sample_name)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
