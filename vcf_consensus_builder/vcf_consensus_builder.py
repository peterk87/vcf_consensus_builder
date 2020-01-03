"""Main module."""
import logging
from io import TextIOWrapper

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from vcf_consensus_builder.io import read_vcf

logger = logging.getLogger(__name__)


def replace_low_depth_positions(ref_fasta: str,
                                depths_file: str,
                                no_coverage: int = 0,
                                low_coverage: int = 5,
                                no_cov_char: str = '-',
                                low_cov_char: str = 'N',
                                *args,
                                **kwargs) -> SeqRecord:
    """Replace no and low coverage depth positions with specified characters.

    Args:
        ref_fasta: FASTA file path
        depths_file: Samtools depth tab-delimited file path
        no_coverage: No coverage depth threshold. At and below this coverage positions are changed to the `no_cov_char`
        low_coverage: Low coverage depth threshold. Below this coverage positions are changed to the `low_cov_char`
        no_cov_char: Character to substitute at no coverage positions (default: "-")
        low_cov_char: Character to substitute at low coverage positions (default: "N")

    Returns:
        SeqRecord with low coverage depth positions substituted with `unmapped_char`
    """
    assert len(no_cov_char) == 1, '"no_cov_char" must be a str of length of 1, e.g. "-"/single dash character.'
    assert len(low_cov_char) == 1, '"low_cov_char" must be a str of length of 1, e.g. "-"/single dash character.'
    assert no_coverage <= low_coverage, '"low_coverage" must be greater than or equal to "no_coverage"'
    seq_rec = SeqIO.read(ref_fasta,
                         format='fasta')
    logger.debug(f'seq_rec length: {len(seq_rec)}')
    df: pd.DataFrame = pd.read_csv(depths_file,
                                   header=None,
                                   names=['genome', 'position', 'coverage'], sep='\t')
    no_coverage_positions: pd.Series = df[df.coverage <= no_coverage].position - 1
    logger.info(f'No ({no_coverage}X) coverage positions: {no_coverage_positions.size}')
    low_coverage_positions: pd.Series = df[(df.coverage > no_coverage) & (df.coverage < low_coverage)].position - 1
    logger.info(f'Low (<{low_coverage}X) coverage positions: {low_coverage_positions.size}')
    if low_coverage_positions.size == 0 and no_coverage_positions.size == 0:
        logger.info(f'No positions with low (<{low_coverage}X) or no ({no_coverage}X) coverage. '
                    f'No need to mask any positions in the reference sequence')
        return seq_rec
    mutable_seq = seq_rec.seq.tomutable()
    for position in low_coverage_positions:
        mutable_seq[position] = low_cov_char
    for position in no_coverage_positions:
        mutable_seq[position] = no_cov_char
    seq_rec.seq = mutable_seq.toseq()
    return seq_rec


def consensus_segment(seq: str,
                      curr_position: int,
                      ref_variant: str,
                      alt_variant: str,
                      prev_position: int = 0) -> (str, int):
    segment_seq = seq[prev_position:(curr_position - 1)] + alt_variant
    next_position = curr_position + len(ref_variant) - 1
    if len(alt_variant) != len(ref_variant):
        logger.info(f'At position {curr_position}, ALT="{alt_variant}" '
                    f'(n={len(alt_variant)} VS REF="{ref_variant}" (n={len(ref_variant)})')
        logger.info(f'Previous position={prev_position} | curr_position={curr_position} | next={next_position}')
    return segment_seq, next_position


def create_cons_seq(seq: str, df_vcf: pd.DataFrame) -> str:
    """Create consensus sequence given a reference sequence and a table of a Snippy vcf_to_tab

    Args:
        seq: Reference sequence
        df_vcf: VCF file DataFrame

    Returns:
        Consensus sequence
    """
    segments = []
    prev_position = 0
    for _, curr_var in df_vcf.iterrows():
        if prev_position > curr_var.POS - 1:
            logger.warning(
                f'Skipping variant (ALT={curr_var.ALT}) at {curr_var.POS} (previous position ({prev_position}) >= POS)')
            continue
        segment, prev_position = consensus_segment(seq=seq,
                                                   curr_position=curr_var.POS,
                                                   ref_variant=curr_var.REF,
                                                   alt_variant=curr_var.ALT,
                                                   prev_position=prev_position)
        segments.append(segment)
    # append the rest of the reference sequence
    segments.append(seq[prev_position:])
    return ''.join(segments)


def consensus(ref_fasta,
              vcf_file,
              depths_file,
              sample_name,
              output_fasta,
              no_coverage: int = 0,
              low_coverage: int = 5,
              no_cov_char: str = '-',
              low_cov_char: str = 'N',
              *args,
              **kwargs):
    ref_seq_record: SeqRecord = replace_low_depth_positions(**locals())

    if logger.isEnabledFor(logging.DEBUG):
        from collections import Counter
        logger.debug(f'ref_seq_record: {Counter(str(ref_seq_record.seq)).most_common()}')
    df_vcf_tsv: pd.DataFrame = read_vcf(vcf_file)
    logger.debug(f'df_vcf_tsv shape: {df_vcf_tsv.shape}')
    if sample_name is None:
        from .io import VCF_COL_DTYPES
        sample_name = list(set(df_vcf_tsv.columns) - set(VCF_COL_DTYPES.keys()))[0]
    consensus_seq: str = create_cons_seq(str(ref_seq_record.seq), df_vcf_tsv)
    logger.debug(f'consensus_seq length: {len(consensus_seq)}')
    with open(output_fasta, 'w') if not isinstance(output_fasta, TextIOWrapper) else output_fasta as f:
        f.write(f'>{sample_name} ref="{ref_seq_record.id} {ref_seq_record.description}"\n')
        # 70 characters per line to keep in line with NCBI FASTA output
        chars_per_line = 70
        for i in range(0, len(consensus_seq), chars_per_line):
            f.write(f'{consensus_seq[i:i + chars_per_line]}\n')
