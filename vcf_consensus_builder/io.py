"""IO functions"""
from os import PathLike
from typing import Dict

import pandas as pd

VCF_COL_DTYPES: Dict = dict(CHROM='category',
                            POS='uint32',
                            ID='category',
                            REF='category',
                            ALT='category',
                            QUAL=float,
                            FILTER='category',
                            INFO=str,
                            FORMAT=str)


def read_vcf(vcf_file: PathLike) -> pd.DataFrame:
    """Read VCF file into a DataFrame"""
    with open(vcf_file) as fh:
        vcf_cols = []
        for line in fh:
            if line.startswith('#CHROM'):
                vcf_cols = line[1:].strip().split('\t')
                break
        df = pd.read_table(fh,
                           comment='#',
                           header=None,
                           names=vcf_cols,
                           dtype=VCF_COL_DTYPES)
    return df
