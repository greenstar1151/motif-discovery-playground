#!/usr/bin/env python3
"""
Morbius Benchmark Accuracy Evaluation Script

Compare predicted motif positions with ground truth (GT) file.

Usage:
    python scripts/evaluate_accuracy.py <gt_file.tsv> <result_file.tsv> [options]

Options:
    --tolerance N       Position tolerance for fuzzy matching (default: 0)
    --motif-length N    Motif length for overlap calculation (default: infer from GT)
    --output FILE       Save detailed results to CSV file
    --verbose           Print per-sequence comparison details

Output:
    - Exact match accuracy
    - Overlap match accuracy (predicted within GT motif region)
    - Tolerance match accuracy
    - Position error statistics
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np


def load_tsv(filepath: str) -> pd.DataFrame:
    """Load TSV file with seq_id, position, motif columns."""
    df = pd.read_csv(filepath, sep='\t', dtype={'seq_id': str})
    
    # Validate required columns
    required_cols = ['seq_id', 'position', 'motif']
    missing = set(required_cols) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    return df


def evaluate_accuracy(
    gt_df: pd.DataFrame,
    result_df: pd.DataFrame,
    tolerance: int = 0,
    motif_length: int | None = None
) -> dict:
    """
    Compare predictions with ground truth.
    
    Args:
        gt_df: Ground truth DataFrame
        result_df: Prediction result DataFrame
        tolerance: Position tolerance for fuzzy matching
        motif_length: Motif length for overlap calculation (inferred from GT if None)
    
    Returns:
        Dictionary with accuracy metrics
    """
    # Merge on seq_id
    merged = gt_df.merge(
        result_df,
        on='seq_id',
        how='left',
        suffixes=('_gt', '_pred')
    )
    
    # Count missing predictions
    missing_preds = merged['position_pred'].isna().sum()
    
    # Filter to sequences with predictions
    has_pred = merged['position_pred'].notna()
    matched = merged[has_pred].copy()
    
    if len(matched) == 0:
        return {
            'total_sequences': len(gt_df),
            'predicted_sequences': 0,
            'missing_predictions': len(gt_df),
            'exact_match': 0,
            'exact_accuracy': 0.0,
            'overlap_match': 0,
            'overlap_accuracy': 0.0,
            'tolerance_match': 0,
            'tolerance_accuracy': 0.0,
        }
    
    # Calculate position difference
    matched['pos_diff'] = (matched['position_pred'] - matched['position_gt']).abs()
    
    # Infer GT motif length (IMPLANTLENGTH) from GT if not provided
    if motif_length is None:
        motif_length = int(matched['motif_gt'].str.len().median())
    
    # Exact match
    exact_match = (matched['pos_diff'] == 0).sum()
    
    # Overlap match: predicted position falls within the implanted motif region
    # i.e., gt_pos <= pred_pos < gt_pos + IMPLANTLENGTH
    # This means: pred_pos - gt_pos >= 0 AND pred_pos - gt_pos < IMPLANTLENGTH
    pos_offset = matched['position_pred'] - matched['position_gt']
    overlap_match = ((pos_offset >= 0) & (pos_offset < motif_length)).sum()
    
    # Tolerance match
    tolerance_match = (matched['pos_diff'] <= tolerance).sum()
    
    # Position error statistics
    pos_errors = matched['pos_diff'].values
    
    total = len(gt_df)
    predicted = len(matched)
    
    return {
        'total_sequences': total,
        'predicted_sequences': predicted,
        'missing_predictions': int(missing_preds),
        'motif_length_gt': motif_length,
        
        'exact_match': int(exact_match),
        'exact_accuracy': 100.0 * exact_match / total,
        
        'overlap_match': int(overlap_match),
        'overlap_accuracy': 100.0 * overlap_match / total,
        
        'tolerance': tolerance,
        'tolerance_match': int(tolerance_match),
        'tolerance_accuracy': 100.0 * tolerance_match / total,
        
        'mean_pos_error': float(np.mean(pos_errors)),
        'median_pos_error': float(np.median(pos_errors)),
        'max_pos_error': int(np.max(pos_errors)),
        'min_pos_error': int(np.min(pos_errors)),
        'std_pos_error': float(np.std(pos_errors)),
    }


def print_results(metrics: dict, verbose: bool = False):
    """Print evaluation results in a formatted way."""
    print("\n" + "=" * 60)
    print("  Morbius Benchmark Accuracy Evaluation")
    print("=" * 60)
    
    print(f"\n[Dataset Statistics]")
    print(f"  Total sequences      : {metrics['total_sequences']:,}")
    print(f"  Predicted sequences  : {metrics['predicted_sequences']:,}")
    print(f"  Missing predictions  : {metrics['missing_predictions']:,}")
    print(f"  GT motif length      : {metrics['motif_length_gt']} bp")
    
    print(f"\n[Accuracy Metrics]")
    print(f"  Exact match          : {metrics['exact_match']:,} / {metrics['total_sequences']:,}")
    print(f"                       : {metrics['exact_accuracy']:.2f}%")
    
    print(f"\n  Overlap match        : {metrics['overlap_match']:,} / {metrics['total_sequences']:,}")
    print(f"  (within GT region)   : {metrics['overlap_accuracy']:.2f}%")
    
    if metrics['tolerance'] > 0:
        print(f"\n  Tolerance match      : {metrics['tolerance_match']:,} / {metrics['total_sequences']:,}")
        print(f"  (±{metrics['tolerance']} bp)        : {metrics['tolerance_accuracy']:.2f}%")
    
    print(f"\n[Position Error Statistics]")
    print(f"  Mean error           : {metrics['mean_pos_error']:.2f} bp")
    print(f"  Median error         : {metrics['median_pos_error']:.2f} bp")
    print(f"  Std dev              : {metrics['std_pos_error']:.2f} bp")
    print(f"  Min / Max            : {metrics['min_pos_error']} / {metrics['max_pos_error']} bp")
    
    print("\n" + "=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate motif discovery accuracy against ground truth',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('gt_file', help='Ground truth TSV file')
    parser.add_argument('result_file', help='Prediction result TSV file')
    parser.add_argument('--tolerance', type=int, default=0,
                        help='Position tolerance for fuzzy matching (default: 0)')
    parser.add_argument('--motif-length', type=int, default=None,
                        help='Motif length for overlap calculation (default: infer from GT)')
    parser.add_argument('--output', type=str, default=None,
                        help='Save detailed merged results to CSV file')
    parser.add_argument('--verbose', action='store_true',
                        help='Print per-sequence comparison details')
    parser.add_argument('--json', action='store_true',
                        help='Output metrics as JSON')
    
    args = parser.parse_args()
    
    # Check files exist
    if not Path(args.gt_file).exists():
        print(f"Error: Ground truth file not found: {args.gt_file}", file=sys.stderr)
        return 1
    
    if not Path(args.result_file).exists():
        print(f"Error: Result file not found: {args.result_file}", file=sys.stderr)
        return 1
    
    try:
        # Load data
        gt_df = load_tsv(args.gt_file)
        result_df = load_tsv(args.result_file)
        
        print(f"[Info] Loaded GT: {len(gt_df):,} sequences")
        print(f"[Info] Loaded Result: {len(result_df):,} sequences")
        
        # Evaluate
        metrics = evaluate_accuracy(
            gt_df, result_df,
            tolerance=args.tolerance,
            motif_length=args.motif_length
        )
        
        if args.json:
            import json
            print(json.dumps(metrics, indent=2))
        else:
            print_results(metrics, verbose=args.verbose)
        
        # Save detailed results if requested
        if args.output:
            merged = gt_df.merge(
                result_df,
                on='seq_id',
                how='left',
                suffixes=('_gt', '_pred')
            )
            merged['pos_diff'] = (merged['position_pred'] - merged['position_gt']).abs()
            merged.to_csv(args.output, index=False)
            print(f"\n[Info] Detailed results saved to: {args.output}")
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
