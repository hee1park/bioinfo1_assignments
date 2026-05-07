#!/usr/bin/env python3
"""Plot ribosome footprint density around the start codon from fivepcounts-filtered-RPF-siLuc.txt."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def load_counts(path: Path) -> pd.DataFrame:
    cols = [
        'chrom',
        'read_start',
        'read_end',
        'count',
        'tx_chrom',
        'tx_start',
        'tx_end',
        'tx_id',
        'start_codon',
        'strand',
    ]
    return pd.read_csv(
        path,
        sep='\t',
        header=None,
        names=cols,
        dtype={
            'chrom': str,
            'read_start': int,
            'read_end': int,
            'count': int,
            'tx_chrom': str,
            'tx_start': int,
            'tx_end': int,
            'tx_id': str,
            'start_codon': int,
            'strand': str,
        },
    )


def summarize_position_counts(df: pd.DataFrame, window: int) -> pd.Series:
    df = df.copy()
    df['rel_pos'] = df['read_start'] - df['start_codon']
    df = df[df['rel_pos'].between(-window, window)]
    summary = df.groupby('rel_pos')['count'].sum().reindex(range(-window, window + 1), fill_value=0)
    return summary


def plot_start_codon_density(counts: pd.Series, output_path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(counts.index, counts.values/1000, width=0.8, color="#4f97ca", edgecolor='black')

    ax.axvline(0, color='red', linewidth=1.5, label='start codon')
    #ax.axvspan(-12, 12, color='orange', alpha=0.15)

    ax.set_xlim(counts.index.min() - 1, counts.index.max() + 1)
    ax.set_xlabel('Relative position to start codon of 5\'-end of reads')
    ax.set_ylabel('Raw count (x1000)')
    ax.set_title(title)
    ax.grid(axis='y', linestyle=':', alpha=0.6)
    ax.legend(loc='upper right')

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    print(f'Wrote plot to {output_path}')


def main() -> None:
    parser = argparse.ArgumentParser(description='Plot RPF density around the start codon for siLuc.')
    parser.add_argument('--input', default='data-work/fivepcounts-filtered-RPF-siLuc.txt', help='Input fivepcounts file')
    parser.add_argument('--output', default='results/start_codon_rpf_siLuc.png', help='Output plot PNG file')
    parser.add_argument('--window', type=int, default=50, help='Window around the start codon in nucleotides')
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = load_counts(input_path)
    counts = summarize_position_counts(df, args.window)
    plot_start_codon_density(counts, output_path, 'RPF density around start codon (siLuc)')


if __name__ == '__main__':
    main()
