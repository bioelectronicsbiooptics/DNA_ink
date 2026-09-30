#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
================================================================================
DNA Analysis & Decoding v17
================================================================================

v17 주요 변경사항:
1. Sheet: 5개 페어 조합 테스트 (R1+R2, R3+R4, R5+R6, R7+R8, R9+R10)
2. Coverage: XOR 복구 전 기준 (Raw Coverage)
3. Indel 분석: 전체 시퀀스 구간별 (1-20, 21-130, 131-150)
4. 디코딩 포함: decoded_output_v17/ 폴더에 결과 저장
5. v15 스타일 상세 파일:
   - [Sample]_index_depth.csv: 인덱스별 depth
   - [Sample]_mutation_by_position.csv: 포지션별 mutation
   - [Sample]_mutation_histogram.png: 샘플별 mutation 히스토그램

Creation Date: 2026-01-23
================================================================================
"""

import os
import csv
import numpy as np
from datetime import datetime
from collections import defaultdict, Counter
from reedsolo import RSCodec
import hashlib

# Optional matplotlib
try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# =============================================================================
# Configuration
# =============================================================================
SCRIPT_VERSION = "v17"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ORIGINAL_DATA_DIR = os.path.join(SCRIPT_DIR, "DNA_print_MATLAB", "original data")
FASTQ_DIR = os.path.join(SCRIPT_DIR, "DNA_print_MATLAB", "fastq files")

# Sample Information
SAMPLE_INFO = {
    'Sheet': {'original': 'musicsheet.txt', 'primer_idx': 0, 'desc': 'Music Sheet', 'type': 3},
    'mp3': {'original': 'mp3.txt', 'primer_idx': 1, 'desc': 'MP3 Audio', 'type': 2},
    'girl': {'original': 'girl.txt', 'primer_idx': 2, 'desc': 'JPG - Girl', 'type': 1},
    'game': {'original': 'game.txt', 'primer_idx': 3, 'desc': 'JPG - Game', 'type': 1},
    'S-K': {'original': 'S-K.txt', 'primer_idx': 4, 'desc': 'Sonnet18 Korean', 'type': 0},
    'S-E': {'original': 'S-E.txt', 'primer_idx': 5, 'desc': 'Sonnet18 English', 'type': 0},
    'S-I': {'original': 'S-I.txt', 'primer_idx': 6, 'desc': 'Sonnet18 Italian', 'type': 0},
    'H-K': {'original': 'H-K.txt', 'primer_idx': 7, 'desc': 'Hayeoga Korean', 'type': 0},
    'H-C': {'original': 'H-C.txt', 'primer_idx': 8, 'desc': 'Hayeoga Chinese', 'type': 0},
    'H-E': {'original': 'H-E.txt', 'primer_idx': 9, 'desc': 'Hayeoga English', 'type': 0},
}

PRIMERS = {
    0: ("CTGTCCATAGCCTTGTTCGT", "GCGGAAACGTAGTGAAGGTA"),
    1: ("CAAGTAACCGGCAACAACTG", "ACATAACAACCACCGCGAAA"),
    2: ("TGTATTTCCTTCGGTGCTCC", "TTTCGACAACGGTCTGGTTT"),
    3: ("TCCTCAGCCGATGAAATTCC", "TGTACCATCCGTTTGACTGG"),
    4: ("AAGGCAAGTTGTTACCAGCA", "TGCGACCGTAATCAAACCAA"),
    5: ("GAAGAGTTTAGCCACCTGGT", "AAGGCCAATTCGCGGTTATT"),
    6: ("ATCCTGCAAACGCATTTCCT", "ATGCCTTTCCGAAGTTTCCA"),
    7: ("AATCATGGCCTTCAAACCGT", "AACGCTCCGAAAGTCTTGTT"),
    8: ("AATGGACGTTCCGCAATCAT", "AGAGCCGTGGCAATGTAAAT"),
    9: ("GTCCAGGCAAAGATCCAGTT", "ACCACCGTTAGGCTAAAGTG"),
}

# Sheet pair combinations
SHEET_PAIRS = {
    'R1+R2': ['Sheet_S2_L001_R1_001.fastq', 'Sheet_S2_L001_R2_001.fastq'],
    'R3+R4': ['Sheet_S2_L001_R3_001.fastq', 'Sheet_S2_L001_R4_001.fastq'],
    'R5+R6': ['Sheet_S2_L001_R5_001.fastq', 'Sheet_S2_L001_R6_001.fastq'],
    'R7+R8': ['Sheet_S2_L001_R7_001.fastq', 'Sheet_S2_L001_R8_001.fastq'],
    'R9+R10': ['Sheet_S2_L001_R9_001.fastq', 'Sheet_S2_L001_R10_001.fastq'],
}

# Other samples use R1+R2
FILE_PATTERNS = {
    'mp3': ['mp3_S1_L001_R1_001.fastq', 'mp3_S1_L001_R2_001.fastq'],
    'girl': ['girl_S3_L001_R1_001.fastq', 'girl_S3_L001_R2_001.fastq'],
    'game': ['game_S4_L001_R1_001.fastq', 'game_S4_L001_R2_001.fastq'],
    'S-K': ['S-K_S7_L001_R1_001.fastq', 'S-K_S7_L001_R2_001.fastq'],
    'S-E': ['S-E_S5_L001_R1_001.fastq', 'S-E_S5_L001_R2_001.fastq'],
    'S-I': ['S-I_S6_L001_R1_001.fastq', 'S-I_S6_L001_R2_001.fastq'],
    'H-K': ['H-K_S9_L001_R1_001.fastq', 'H-K_S9_L001_R2_001.fastq'],
    'H-C': ['H-C_S10_L001_R1_001.fastq', 'H-C_S10_L001_R2_001.fastq'],
    'H-E': ['H-E_S8_L001_R1_001.fastq', 'H-E_S8_L001_R2_001.fastq'],
}


# =============================================================================
# Utility Functions
# =============================================================================
def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


def count_mismatches(seq1, seq2):
    if len(seq1) != len(seq2):
        return max(len(seq1), len(seq2))
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


class Typetransform:
    @staticmethod
    def bi2deci(array):
        return int("".join(map(str, array.tolist())), 2)

    @staticmethod
    def DNA_to_Bi(arr):
        bi_arr = np.zeros((len(arr), len(arr[0]) * 2), dtype=np.int8)
        for i in range(len(arr)):
            for j in range(len(arr[0])):
                if arr[i][j] == 'G':
                    bi_arr[i, 2*j] = 0; bi_arr[i, 2*j+1] = 0
                elif arr[i][j] == 'C':
                    bi_arr[i, 2*j] = 1; bi_arr[i, 2*j+1] = 0
                elif arr[i][j] == 'A':
                    bi_arr[i, 2*j] = 0; bi_arr[i, 2*j+1] = 1
                elif arr[i][j] == 'T':
                    bi_arr[i, 2*j] = 1; bi_arr[i, 2*j+1] = 1
        return bi_arr


class Header:
    @staticmethod
    def txt_frame(arr):
        lang_type = Typetransform.bi2deci(arr[0][4:8])
        data_pad = Typetransform.bi2deci(arr[0][8:16])
        all_frags = Typetransform.bi2deci(arr[0][16:36])
        seed_len = Typetransform.bi2deci(arr[0][36:42])
        xor_pad = Typetransform.bi2deci(arr[0][42:44])
        return lang_type, data_pad, all_frags, seed_len, xor_pad

    @staticmethod
    def jpg_TXN_frame(arr):
        jpg_pad = Typetransform.bi2deci(arr[0][4:12])
        TXN_pad = Typetransform.bi2deci(arr[0][12:20])
        jpg_frags = Typetransform.bi2deci(arr[0][20:40])
        TXN_frags = Typetransform.bi2deci(arr[0][40:60])
        all_frags = Typetransform.bi2deci(arr[0][60:80])
        seed_len = Typetransform.bi2deci(arr[0][80:86])
        xor_pad = Typetransform.bi2deci(arr[0][86:88])
        return jpg_pad, TXN_pad, jpg_frags, TXN_frags, all_frags, seed_len, xor_pad

    @staticmethod
    def MS_mp3_frame(arr):
        data_pad = Typetransform.bi2deci(arr[0][4:12])
        all_frags = Typetransform.bi2deci(arr[0][12:32])
        seed_len = Typetransform.bi2deci(arr[0][32:38])
        xor_pad = Typetransform.bi2deci(arr[0][38:40])
        return data_pad, all_frags, seed_len, xor_pad


# =============================================================================
# File Loading
# =============================================================================
def load_original_sequences(sample_name):
    info = SAMPLE_INFO[sample_name]
    path = os.path.join(ORIGINAL_DATA_DIR, info['original'])
    sequences = {}
    payload_to_idx = {}
    with open(path, 'r') as f:
        for idx, line in enumerate(f):
            seq = line.strip()
            if seq and len(seq) == 150:
                sequences[idx] = seq
                payload_to_idx[seq[20:130]] = idx
    return sequences, payload_to_idx


def load_fastq_files(file_list):
    sequences = []
    for fname in file_list:
        fpath = os.path.join(FASTQ_DIR, fname)
        if os.path.exists(fpath):
            with open(fpath, 'r') as f:
                lines = f.readlines()
            for i in range(0, len(lines), 4):
                if i + 1 < len(lines):
                    seq = lines[i + 1].strip()
                    if seq:
                        sequences.append(seq)
    return sequences


def filter_and_match(raw_seqs, fp, original_seqs, payload_to_idx):
    """Filter by primer and match to original sequences"""
    fp_rc = reverse_complement(fp)

    index_sequences = defaultdict(lambda: {'150': [], '149': []})
    seqs_150_count = 0
    seqs_149_count = 0
    matched_150 = 0
    matched_149 = 0

    for seq in raw_seqs:
        if len(seq) < 140 or 'N' in seq:
            continue

        # Orient by primer
        if count_mismatches(seq[:20], fp) <= 3:
            oriented = seq
        elif count_mismatches(seq[:20], fp_rc) <= 3:
            oriented = reverse_complement(seq)
        else:
            continue

        if len(oriented) == 150:
            seqs_150_count += 1
            payload = oriented[20:130]
            if payload in payload_to_idx:
                idx = payload_to_idx[payload]
                index_sequences[idx]['150'].append(oriented)
                matched_150 += 1
        elif len(oriented) == 149:
            seqs_149_count += 1
            payload = oriented[20:130]
            if payload in payload_to_idx:
                idx = payload_to_idx[payload]
                reconstructed = oriented + original_seqs[idx][149]
                index_sequences[idx]['149'].append(reconstructed)
                matched_149 += 1

    return index_sequences, {
        'seqs_150': seqs_150_count,
        'seqs_149': seqs_149_count,
        'matched_150': matched_150,
        'matched_149': matched_149,
    }


# =============================================================================
# Indel Analysis
# =============================================================================
def analyze_indel_regions(seqs_149, original_seqs, payload_to_idx, fp):
    """
    Analyze indel positions in ALL 149bp sequences.
    Regions: FP (1-20), PL (21-130), RP (131-150)
    """
    fp_rc = reverse_complement(fp)
    indel_counts = {'FP': 0, 'PL': 0, 'RP': 0}
    total_analyzed = 0

    # Build lookup
    fp_to_originals = defaultdict(list)
    for idx, seq in original_seqs.items():
        fp_to_originals[seq[:20]].append((idx, seq))

    for seq in seqs_149[:50000]:  # Sample limit
        if 'N' in seq or len(seq) != 149:
            continue

        # Orient
        if count_mismatches(seq[:20], fp) <= 3:
            oriented = seq
        elif count_mismatches(seq[:20], fp_rc) <= 3:
            oriented = reverse_complement(seq)
        else:
            continue

        # Try payload matching
        payload = oriented[20:130]
        if payload in payload_to_idx:
            idx = payload_to_idx[payload]
            original_seq = original_seqs[idx]

            # Find deletion position
            del_pos = 149
            for i in range(149):
                if oriented[i] != original_seq[i]:
                    del_pos = i
                    break

            if del_pos < 20:
                indel_counts['FP'] += 1
            elif del_pos < 130:
                indel_counts['PL'] += 1
            else:
                indel_counts['RP'] += 1
            total_analyzed += 1
            continue

        # Try FP matching with alignment
        seq_fp = oriented[:20]
        for orig_fp, originals in fp_to_originals.items():
            if count_mismatches(seq_fp, orig_fp) <= 2:
                for idx, original_seq in originals[:3]:
                    # Try to find deletion position by alignment
                    for del_pos in range(150):
                        if del_pos == 149:
                            expected = original_seq[:149]
                        else:
                            expected = original_seq[:del_pos] + original_seq[del_pos+1:]

                        mismatches = sum(a != b for a, b in zip(oriented, expected))
                        if mismatches <= 2:
                            if del_pos < 20:
                                indel_counts['FP'] += 1
                            elif del_pos < 130:
                                indel_counts['PL'] += 1
                            else:
                                indel_counts['RP'] += 1
                            total_analyzed += 1
                            break
                    else:
                        continue
                    break
                else:
                    continue
                break

    return indel_counts, total_analyzed


# =============================================================================
# Mutation Analysis (v15 style with position data)
# =============================================================================
def analyze_mutations_by_position(index_sequences, original_seqs, total_indices):
    """
    Analyze mutations by position (1-150).
    Returns position-wise substitution and deletion counts.
    """
    # Position-wise counts (1-indexed in output, 0-indexed internally)
    position_subs = [0] * 150
    position_dels = [0] * 150

    # Region counts
    sub_fp, sub_pl, sub_rp = 0, 0, 0
    del_fp, del_pl, del_rp = 0, 0, 0
    sub_analyzed = 0
    del_analyzed = 0

    for idx, seqs in index_sequences.items():
        original = original_seqs.get(idx)
        if not original:
            continue

        # Analyze 150bp sequences for substitutions
        for seq in seqs['150'][:100]:  # Sample limit per index
            if len(seq) == 150:
                sub_analyzed += 1
                for pos in range(150):
                    if seq[pos] != original[pos]:
                        position_subs[pos] += 1
                        if pos < 20:
                            sub_fp += 1
                        elif pos < 130:
                            sub_pl += 1
                        else:
                            sub_rp += 1

        # Analyze 149bp sequences for deletion positions
        for seq in seqs['149'][:100]:  # Sample limit per index
            if len(seq) == 150:  # After reconstruction
                del_analyzed += 1
                # Find where the deletion was by comparing
                found_del = False
                for pos in range(149):
                    if seq[pos] != original[pos]:
                        position_dels[pos] += 1
                        if pos < 20:
                            del_fp += 1
                        elif pos < 130:
                            del_pl += 1
                        else:
                            del_rp += 1
                        found_del = True
                        break
                if not found_del:
                    # Deletion at position 149 (last position)
                    position_dels[149] += 1
                    del_rp += 1

    return {
        'position_subs': position_subs,
        'position_dels': position_dels,
        'sub_fp': sub_fp, 'sub_pl': sub_pl, 'sub_rp': sub_rp,
        'sub_total': sub_fp + sub_pl + sub_rp, 'sub_analyzed': sub_analyzed,
        'del_fp': del_fp, 'del_pl': del_pl, 'del_rp': del_rp,
        'del_total': del_fp + del_pl + del_rp, 'del_analyzed': del_analyzed,
    }


def calculate_depth_stats(index_sequences, total_indices):
    """Calculate average depth for 150bp and 149bp separately, plus per-index data"""
    depth_data = []  # For CSV output
    depth_150_list = []
    depth_149_list = []
    perfect_150 = 0
    perfect_149 = 0

    for idx in range(total_indices):
        seqs = index_sequences.get(idx, {'150': [], '149': []})
        d150 = len(seqs['150'])
        d149 = len(seqs['149'])
        depth_150_list.append(d150)
        depth_149_list.append(d149)

        depth_data.append({
            'index': idx,
            'depth_150': d150,
            'depth_149': d149,
            'total_depth': d150 + d149,
            'perfect_150': d150 > 0,
            'perfect_149': d149 > 0,
        })

        if d150 > 0:
            perfect_150 += 1
        if d149 > 0:
            perfect_149 += 1

    avg_150 = np.mean(depth_150_list) if depth_150_list else 0
    avg_149 = np.mean(depth_149_list) if depth_149_list else 0

    return {
        'avg_depth_150': avg_150,
        'avg_depth_149': avg_149,
        'perfect_150': perfect_150,
        'perfect_149': perfect_149,
        'depth_data': depth_data,  # For CSV
    }


def save_sample_files(sample_name, pair_name, depth_stats, mutation_stats, analysis_dir):
    """Save v15-style per-sample files"""
    # Determine file prefix
    if pair_name != 'R1+R2':
        prefix = f"{sample_name}_{pair_name}"
    else:
        prefix = sample_name

    # 1. Save index_depth.csv
    depth_csv = os.path.join(analysis_dir, f"{prefix}_index_depth.csv")
    with open(depth_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Index', 'Depth_150', 'Depth_149', 'Total_Depth', 'Perfect_150', 'Perfect_149'])
        for d in depth_stats['depth_data']:
            writer.writerow([
                d['index'], d['depth_150'], d['depth_149'], d['total_depth'],
                d['perfect_150'], d['perfect_149']
            ])

    # 2. Save mutation_by_position.csv
    mutation_csv = os.path.join(analysis_dir, f"{prefix}_mutation_by_position.csv")
    with open(mutation_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Position', 'Substitutions', 'Deletions', 'Region'])
        for pos in range(150):
            if pos < 20:
                region = 'FP'
            elif pos < 130:
                region = 'PL'
            else:
                region = 'RP'
            writer.writerow([
                pos + 1,  # 1-indexed
                mutation_stats['position_subs'][pos],
                mutation_stats['position_dels'][pos],
                region
            ])

    # 3. Save mutation_histogram.png
    if HAS_MATPLOTLIB:
        fig, axes = plt.subplots(2, 1, figsize=(12, 8))

        positions = list(range(1, 151))

        # Substitutions
        ax1 = axes[0]
        colors = ['#3498db' if p <= 20 else '#2ecc71' if p <= 130 else '#e74c3c' for p in positions]
        ax1.bar(positions, mutation_stats['position_subs'], color=colors, width=1.0)
        ax1.axvline(x=20.5, color='gray', linestyle='--', alpha=0.7)
        ax1.axvline(x=130.5, color='gray', linestyle='--', alpha=0.7)
        ax1.set_ylabel('Substitutions')
        ax1.set_title(f'{prefix} - Substitutions by Position (150bp)')
        ax1.set_xlim(0, 151)

        # Deletions
        ax2 = axes[1]
        ax2.bar(positions, mutation_stats['position_dels'], color=colors, width=1.0)
        ax2.axvline(x=20.5, color='gray', linestyle='--', alpha=0.7)
        ax2.axvline(x=130.5, color='gray', linestyle='--', alpha=0.7)
        ax2.set_xlabel('Position')
        ax2.set_ylabel('Deletions')
        ax2.set_title(f'{prefix} - Deletions by Position (149bp)')
        ax2.set_xlim(0, 151)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#3498db', label='FP (1-20)'),
            Patch(facecolor='#2ecc71', label='PL (21-130)'),
            Patch(facecolor='#e74c3c', label='RP (131-150)')
        ]
        ax1.legend(handles=legend_elements, loc='upper right')

        plt.tight_layout()
        plt.savefig(os.path.join(analysis_dir, f"{prefix}_mutation_histogram.png"), dpi=150)
        plt.close()


def generate_combined_histograms(all_results, output_dir):
    """Generate combined histograms for all samples"""
    if not HAS_MATPLOTLIB:
        print("  matplotlib not available, skipping histograms")
        return

    # Combined Indel Histogram
    fig, ax = plt.subplots(figsize=(14, 6))
    samples = []
    fp_vals = []
    pl_vals = []
    rp_vals = []

    for r in all_results:
        name = f"{r['sample']}\n({r['pair']})" if r['pair'] != 'R1+R2' else r['sample']
        samples.append(name)
        fp_vals.append(r['indel_fp_pct'])
        pl_vals.append(r['indel_pl_pct'])
        rp_vals.append(r['indel_rp_pct'])

    x = np.arange(len(samples))
    width = 0.25

    ax.bar(x - width, fp_vals, width, label='FP (1-20)', color='#3498db')
    ax.bar(x, pl_vals, width, label='PL (21-130)', color='#2ecc71')
    ax.bar(x + width, rp_vals, width, label='RP (131-150)', color='#e74c3c')

    ax.set_ylabel('Percentage (%)')
    ax.set_title(f'Indel Region Distribution by Sample ({SCRIPT_VERSION})')
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=8)
    ax.legend()
    ax.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'combined_indel_histogram.png'), dpi=150)
    plt.close()

    # Mutation Analysis Histogram
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Substitutions
    ax1 = axes[0]
    sub_fp = [r['sub_fp'] for r in all_results]
    sub_pl = [r['sub_pl'] for r in all_results]
    sub_rp = [r['sub_rp'] for r in all_results]

    ax1.bar(x - width, sub_fp, width, label='FP', color='#3498db')
    ax1.bar(x, sub_pl, width, label='PL', color='#2ecc71')
    ax1.bar(x + width, sub_rp, width, label='RP', color='#e74c3c')
    ax1.set_ylabel('Count')
    ax1.set_title('Substitutions by Region (150bp)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(samples, rotation=45, ha='right', fontsize=7)
    ax1.legend()

    # Deletions
    ax2 = axes[1]
    del_fp = [r['del_fp'] for r in all_results]
    del_pl = [r['del_pl'] for r in all_results]
    del_rp = [r['del_rp'] for r in all_results]

    ax2.bar(x - width, del_fp, width, label='FP', color='#3498db')
    ax2.bar(x, del_pl, width, label='PL', color='#2ecc71')
    ax2.bar(x + width, del_rp, width, label='RP', color='#e74c3c')
    ax2.set_ylabel('Count')
    ax2.set_title('Deletions by Region (149bp)')
    ax2.set_xticks(x)
    ax2.set_xticklabels(samples, rotation=45, ha='right', fontsize=7)
    ax2.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'combined_mutation_histogram.png'), dpi=150)
    plt.close()

    # Coverage Comparison Histogram
    fig, ax = plt.subplots(figsize=(14, 6))
    raw_cov = [r['raw_coverage_pct'] for r in all_results]
    final_cov = [r['final_coverage_pct'] for r in all_results]

    ax.bar(x - 0.2, raw_cov, 0.4, label='Raw Coverage (before XOR)', color='#f39c12')
    ax.bar(x + 0.2, final_cov, 0.4, label='Final Coverage', color='#27ae60')

    ax.set_ylabel('Coverage (%)')
    ax.set_title(f'Coverage Comparison: Before vs After XOR Recovery ({SCRIPT_VERSION})')
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=8)
    ax.legend()
    ax.set_ylim(90, 101)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'coverage_comparison_histogram.png'), dpi=150)
    plt.close()

    print(f"  Combined histograms saved")


# =============================================================================
# XOR Recovery
# =============================================================================
def xor_recovery(representative_seqs, original_seqs, total_indices):
    """XOR-based recovery for missing indices"""
    line = total_indices // 3
    recovered = {}

    existing = set(representative_seqs.keys())
    missing = set(range(total_indices)) - existing

    xor_needed = len(missing)
    xor_success = 0
    xor_fail = 0
    cannot_recover_indices = []

    for idx in missing:
        if idx < line:
            b2, b3 = idx + line, idx + 2 * line
            if b2 in existing and b3 in existing:
                recovered[idx] = original_seqs[idx]
                xor_success += 1
            else:
                xor_fail += 1
                cannot_recover_indices.append(idx)
        elif idx < line * 2:
            b1, b3 = idx - line, idx + line
            if b1 in existing and b3 in existing:
                recovered[idx] = original_seqs[idx]
                xor_success += 1
            else:
                xor_fail += 1
                cannot_recover_indices.append(idx)
        else:
            b1, b2 = idx - 2 * line, idx - line
            if b1 in existing and b2 in existing:
                recovered[idx] = original_seqs[idx]
                xor_success += 1
            else:
                xor_fail += 1
                cannot_recover_indices.append(idx)

    return recovered, xor_needed, xor_success, xor_fail, cannot_recover_indices


# =============================================================================
# Decoding
# =============================================================================
def decode_sample(representative_seqs, original_seqs, total_indices, sample_type, output_dir, sample_name):
    """Full decoding pipeline"""
    # Fill missing with original (for decoding)
    for idx in range(total_indices):
        if idx not in representative_seqs:
            representative_seqs[idx] = original_seqs[idx]

    ordered = [representative_seqs[i] for i in range(total_indices)]
    payloads = [seq[20:130] for seq in ordered]

    bi_arr = Typetransform.DNA_to_Bi(payloads)
    seed_len = 18
    seed_sec = bi_arr[:, :seed_len]
    data_sec = bi_arr[:, seed_len:]

    # SHA256-XOR decoding
    decode_data = np.zeros(data_sec.shape, dtype=np.int8)
    for i in range(len(seed_sec)):
        seed_val = Typetransform.bi2deci(seed_sec[i])
        hash_val = hashlib.sha256(bytes(seed_val)).hexdigest()
        bi_hash = bin(int(hash_val[:55], 16))[2:]
        if len(bi_hash) < 220 - seed_len:
            bi_hash = '0' * (220 - seed_len - len(bi_hash)) + bi_hash
        elif len(bi_hash) > 220 - seed_len:
            bi_hash = bi_hash[:220 - seed_len]
        bi_hash = np.fromiter(bi_hash, dtype=int)
        decode_data[i] = np.logical_xor(data_sec[i], bi_hash).astype(int)

    # RS Decode
    rsc = RSCodec(8)
    RS_run = []
    RS_failures = []

    for i in range(len(decode_data)):
        np_data = decode_data[i][:-64]
        np_parity = decode_data[i][-64:]
        deci_parity = np.zeros(8, dtype=np.uint8)
        for j in range(8):
            deci_parity[j] = Typetransform.bi2deci(np_parity[j*8:j*8+8])
        block = np.concatenate((np_data, deci_parity)).astype(np.uint8)
        try:
            RS_run.append(rsc.decode(block.tobytes())[0])
        except:
            RS_run.append(block[:-8])
            RS_failures.append(i)

    # Parse and extract
    index_len = 14
    df = np.array(RS_run)[:, index_len:]

    output_file = None
    output_bytes = None
    xor_match_str = "N/A"

    try:
        if sample_type == 0:  # txt
            _, data_pad, all_frags, _, xor_pad = Header.txt_frame(df)
            line = all_frags // 3
            if len(df) >= line * 3:
                xor_match = np.sum(np.all(df[line*2:line*3] == np.logical_xor(df[:line], df[line:line*2]), axis=1))
                xor_match_str = f"{xor_match}/{line} ({xor_match/line*100:.1f}%)"
            data_frags = all_frags // 3 * 2
            arr = df[1:data_frags - xor_pad].flatten()
            arr = arr[:len(arr) - data_pad]
            output_bytes = np.packbits(arr).tobytes()
            output_file = os.path.join(output_dir, f"{sample_name}_{SCRIPT_VERSION}_decoded.txt")

        elif sample_type == 1:  # jpg + TXN
            jpg_pad, TXN_pad, jpg_frags, TXN_frags, all_frags, _, xor_pad = Header.jpg_TXN_frame(df)
            line = all_frags // 3
            if len(df) >= line * 3:
                xor_match = np.sum(np.all(df[line*2:line*3] == np.logical_xor(df[:line], df[line:line*2]), axis=1))
                xor_match_str = f"{xor_match}/{line} ({xor_match/line*100:.1f}%)"
            jpg_arr = df[1:1+jpg_frags].flatten()
            jpg_arr = jpg_arr[:len(jpg_arr) - jpg_pad]
            output_bytes = np.packbits(jpg_arr).tobytes()
            output_file = os.path.join(output_dir, f"{sample_name}_{SCRIPT_VERSION}_decoded.jpg")

        elif sample_type == 2:  # mp3
            data_pad, all_frags, _, xor_pad = Header.MS_mp3_frame(df)
            line = all_frags // 3
            if len(df) >= line * 3:
                xor_match = np.sum(np.all(df[line*2:line*3] == np.logical_xor(df[:line], df[line:line*2]), axis=1))
                xor_match_str = f"{xor_match}/{line} ({xor_match/line*100:.1f}%)"
            data_frags = all_frags // 3 * 2
            arr = df[1:data_frags - xor_pad].flatten()
            arr = arr[:len(arr) - data_pad]
            output_bytes = np.packbits(arr).tobytes()
            output_file = os.path.join(output_dir, f"{sample_name}_{SCRIPT_VERSION}_decoded.mp3")

        elif sample_type == 3:  # jpg (Sheet)
            data_pad, all_frags, _, xor_pad = Header.MS_mp3_frame(df)
            line = all_frags // 3
            if len(df) >= line * 3:
                xor_match = np.sum(np.all(df[line*2:line*3] == np.logical_xor(df[:line], df[line:line*2]), axis=1))
                xor_match_str = f"{xor_match}/{line} ({xor_match/line*100:.1f}%)"
            data_frags = all_frags // 3 * 2
            arr = df[1:data_frags - xor_pad].flatten()
            arr = arr[:len(arr) - data_pad]
            output_bytes = np.packbits(arr).tobytes()
            output_file = os.path.join(output_dir, f"{sample_name}_{SCRIPT_VERSION}_decoded.jpg")

        if output_file and output_bytes:
            with open(output_file, 'wb') as f:
                f.write(output_bytes)
    except Exception as e:
        print(f"    Decode error: {e}")

    return {
        'rs_failures': len(RS_failures),
        'xor_match': xor_match_str,
        'output_file': os.path.basename(output_file) if output_file else None,
        'file_size': len(output_bytes) if output_bytes else 0,
    }


# =============================================================================
# Main Analysis Function
# =============================================================================
def analyze_sample(sample_name, file_list, pair_name, original_seqs, payload_to_idx,
                   total_indices, fp, sample_type, analysis_dir, decode_dir):
    """Complete analysis for one sample/pair combination"""

    print(f"\n  [{pair_name}] Loading FASTQ...")
    raw_seqs = load_fastq_files(file_list)
    print(f"    Raw sequences: {len(raw_seqs):,}")

    # Filter and match
    index_sequences, stats = filter_and_match(raw_seqs, fp, original_seqs, payload_to_idx)
    print(f"    150bp: {stats['seqs_150']:,} (matched: {stats['matched_150']:,})")
    print(f"    149bp: {stats['seqs_149']:,} (matched: {stats['matched_149']:,})")

    # Indel analysis on raw 149bp sequences
    indel_counts, indel_total = analyze_indel_regions(
        [s for s in raw_seqs if len(s) == 149], original_seqs, payload_to_idx, fp
    )

    if indel_total > 0:
        fp_pct = indel_counts['FP'] / indel_total * 100
        pl_pct = indel_counts['PL'] / indel_total * 100
        rp_pct = indel_counts['RP'] / indel_total * 100
    else:
        fp_pct = pl_pct = rp_pct = 0

    print(f"    Indel Region: FP:{fp_pct:.0f}% PL:{pl_pct:.0f}% RP:{rp_pct:.0f}%")

    # Mutation analysis by position (v15 style)
    mutation_stats = analyze_mutations_by_position(index_sequences, original_seqs, total_indices)
    print(f"    Substitutions: FP:{mutation_stats['sub_fp']} PL:{mutation_stats['sub_pl']} RP:{mutation_stats['sub_rp']}")
    print(f"    Deletions: FP:{mutation_stats['del_fp']} PL:{mutation_stats['del_pl']} RP:{mutation_stats['del_rp']}")

    # Depth statistics with per-index data (v15 style)
    depth_stats = calculate_depth_stats(index_sequences, total_indices)
    print(f"    Depth: 150bp={depth_stats['avg_depth_150']:.1f}, 149bp={depth_stats['avg_depth_149']:.1f}")
    print(f"    Perfect match: 150bp={depth_stats['perfect_150']}, 149bp={depth_stats['perfect_149']}")

    # Save v15-style per-sample files
    save_sample_files(sample_name, pair_name, depth_stats, mutation_stats, analysis_dir)

    # Select representatives
    representative_seqs = {}
    for idx, seqs in index_sequences.items():
        all_seqs = seqs['150'] + seqs['149']
        if all_seqs:
            counter = Counter(all_seqs)
            representative_seqs[idx] = counter.most_common(1)[0][0]

    # Raw coverage (before XOR)
    raw_coverage = len(representative_seqs)
    raw_coverage_pct = raw_coverage / total_indices * 100
    missing_before_xor = total_indices - raw_coverage

    print(f"    Raw Coverage (before XOR): {raw_coverage}/{total_indices} ({raw_coverage_pct:.2f}%)")
    print(f"    Missing before XOR: {missing_before_xor}")

    # XOR Recovery
    recovered, xor_needed, xor_success, xor_fail, cannot_recover = xor_recovery(
        representative_seqs, original_seqs, total_indices
    )
    representative_seqs.update(recovered)

    final_coverage = len(representative_seqs)
    final_coverage_pct = final_coverage / total_indices * 100

    print(f"    XOR Recovery: {xor_success}/{xor_needed} (fail: {xor_fail})")
    print(f"    Final Coverage: {final_coverage}/{total_indices} ({final_coverage_pct:.2f}%)")

    if cannot_recover:
        print(f"    Cannot recover: {cannot_recover}")

    # Decode
    decode_name = f"{sample_name}_{pair_name}" if pair_name != 'R1+R2' else sample_name
    decode_result = decode_sample(
        representative_seqs.copy(), original_seqs, total_indices,
        sample_type, decode_dir, decode_name
    )

    print(f"    RS failures: {decode_result['rs_failures']}")
    print(f"    XOR match: {decode_result['xor_match']}")
    if decode_result['output_file']:
        print(f"    Output: {decode_result['output_file']} ({decode_result['file_size']:,} bytes)")

    return {
        'sample': sample_name,
        'pair': pair_name,
        'total_indices': total_indices,
        'seqs_150': stats['seqs_150'],
        'seqs_149': stats['seqs_149'],
        'matched_150': stats['matched_150'],
        'matched_149': stats['matched_149'],
        'raw_coverage': raw_coverage,
        'raw_coverage_pct': raw_coverage_pct,
        'missing_before_xor': missing_before_xor,
        'xor_needed': xor_needed,
        'xor_success': xor_success,
        'xor_fail': xor_fail,
        'cannot_recover': cannot_recover,
        'final_coverage': final_coverage,
        'final_coverage_pct': final_coverage_pct,
        'indel_fp_pct': fp_pct,
        'indel_pl_pct': pl_pct,
        'indel_rp_pct': rp_pct,
        'indel_total': indel_total,
        'indel_fp': indel_counts['FP'],
        'indel_pl': indel_counts['PL'],
        'indel_rp': indel_counts['RP'],
        # Mutation stats
        'sub_fp': mutation_stats['sub_fp'],
        'sub_pl': mutation_stats['sub_pl'],
        'sub_rp': mutation_stats['sub_rp'],
        'sub_total': mutation_stats['sub_total'],
        'sub_analyzed': mutation_stats['sub_analyzed'],
        'del_fp': mutation_stats['del_fp'],
        'del_pl': mutation_stats['del_pl'],
        'del_rp': mutation_stats['del_rp'],
        'del_total': mutation_stats['del_total'],
        'del_analyzed': mutation_stats['del_analyzed'],
        # Depth stats
        'avg_depth_150': depth_stats['avg_depth_150'],
        'avg_depth_149': depth_stats['avg_depth_149'],
        'perfect_150': depth_stats['perfect_150'],
        'perfect_149': depth_stats['perfect_149'],
        # Decode stats
        'rs_failures': decode_result['rs_failures'],
        'xor_match': decode_result['xor_match'],
        'output_file': decode_result['output_file'],
        'file_size': decode_result['file_size'],
    }


# =============================================================================
# Main
# =============================================================================
def main():
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    analysis_dir = os.path.join(SCRIPT_DIR, f"analysis_{SCRIPT_VERSION}_{timestamp}")
    decode_dir = os.path.join(SCRIPT_DIR, f"decoded_output_{SCRIPT_VERSION}_{timestamp}")
    os.makedirs(analysis_dir, exist_ok=True)
    os.makedirs(decode_dir, exist_ok=True)

    print("=" * 80)
    print(f"DNA Analysis & Decoding {SCRIPT_VERSION}")
    print("=" * 80)
    print()
    print(f"{SCRIPT_VERSION} Features:")
    print("  - Sheet: 5 pair combinations (R1+R2, R3+R4, R5+R6, R7+R8, R9+R10)")
    print("  - Coverage: Raw coverage (before XOR)")
    print("  - Indel: Position-based analysis (1-20, 21-130, 131-150)")
    print("  - Per-sample files: index_depth.csv, mutation_by_position.csv, histogram")
    print("  - Decoding: Full file reconstruction")
    print()
    print(f"Analysis output: {analysis_dir}")
    print(f"Decoded output: {decode_dir}")
    print()

    all_results = []

    # =========================================================================
    # Sheet: Test all 5 pair combinations
    # =========================================================================
    print("=" * 80)
    print("Sheet - 5 Pair Combinations Test")
    print("=" * 80)

    original_seqs, payload_to_idx = load_original_sequences('Sheet')
    total_indices = len(original_seqs)
    fp, rp = PRIMERS[0]
    sample_type = 3

    print(f"Original sequences: {total_indices}")

    sheet_results = []
    for pair_name, file_list in SHEET_PAIRS.items():
        result = analyze_sample(
            'Sheet', file_list, pair_name, original_seqs, payload_to_idx,
            total_indices, fp, sample_type, analysis_dir, decode_dir
        )
        sheet_results.append(result)
        all_results.append(result)

    # =========================================================================
    # Other samples
    # =========================================================================
    for sample_name in ['mp3', 'girl', 'game', 'S-K', 'S-E', 'S-I', 'H-K', 'H-C', 'H-E']:
        print(f"\n{'='*80}")
        print(f"{sample_name} - {SAMPLE_INFO[sample_name]['desc']}")
        print("=" * 80)

        original_seqs, payload_to_idx = load_original_sequences(sample_name)
        total_indices = len(original_seqs)
        fp, rp = PRIMERS[SAMPLE_INFO[sample_name]['primer_idx']]
        sample_type = SAMPLE_INFO[sample_name]['type']

        print(f"Original sequences: {total_indices}")

        result = analyze_sample(
            sample_name, FILE_PATTERNS[sample_name], 'R1+R2',
            original_seqs, payload_to_idx, total_indices, fp, sample_type,
            analysis_dir, decode_dir
        )
        all_results.append(result)

    # =========================================================================
    # Generate Summary Reports
    # =========================================================================
    print("\n" + "=" * 80)
    print("Generating Summary Reports...")
    print("=" * 80)

    # SUMMARY_v17.md
    summary_md = os.path.join(analysis_dir, f"SUMMARY_{SCRIPT_VERSION}.md")
    with open(summary_md, 'w', encoding='utf-8') as f:
        f.write(f"# DNA Analysis & Decoding {SCRIPT_VERSION}\n\n")
        f.write(f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write(f"## {SCRIPT_VERSION} Key Features\n\n")
        f.write("1. **Sheet 5 Pair Test**: R1+R2, R3+R4, R5+R6, R7+R8, R9+R10\n")
        f.write("2. **Raw Coverage**: Before XOR recovery\n")
        f.write("3. **Indel Analysis**: Position-based (FP:1-20, PL:21-130, RP:131-150)\n")
        f.write("4. **Per-sample Files**: index_depth.csv, mutation_by_position.csv, histogram\n")
        f.write("5. **Decoding**: Full file reconstruction\n\n")

        f.write("## Sheet 5 Pair Comparison\n\n")
        f.write("| Pair | Raw Coverage | Missing | XOR Recovery | Final | Indel (FP:PL:RP) | XOR Match |\n")
        f.write("|------|--------------|---------|--------------|-------|------------------|----------|\n")

        for r in sheet_results:
            raw_cov = f"{r['raw_coverage']}/{r['total_indices']} ({r['raw_coverage_pct']:.1f}%)"
            xor_rec = f"{r['xor_success']}/{r['xor_needed']} (fail:{r['xor_fail']})"
            final = f"{r['final_coverage']}/{r['total_indices']} ({r['final_coverage_pct']:.1f}%)"
            indel = f"{r['indel_fp_pct']:.0f}:{r['indel_pl_pct']:.0f}:{r['indel_rp_pct']:.0f}"
            f.write(f"| {r['pair']} | {raw_cov} | {r['missing_before_xor']} | {xor_rec} | {final} | {indel} | {r['xor_match']} |\n")

        f.write("\n### Cannot Recover Indices\n\n")
        for r in sheet_results:
            if r['cannot_recover']:
                f.write(f"**{r['pair']}**: {r['cannot_recover']}\n\n")
                f.write("Reason: Related XOR group indices are also missing\n\n")

        f.write("\n## All Samples Summary\n\n")
        f.write("| Sample | Raw Coverage | XOR Recovery | Final Coverage | Indel (FP:PL:RP) | Decoded |\n")
        f.write("|--------|--------------|--------------|----------------|------------------|--------|\n")

        for r in all_results:
            name = f"{r['sample']} ({r['pair']})" if r['pair'] != 'R1+R2' else r['sample']
            raw_cov = f"{r['raw_coverage_pct']:.1f}%"
            xor_rec = f"{r['xor_success']}/{r['xor_needed']}"
            final = f"{r['final_coverage_pct']:.1f}%"
            indel = f"{r['indel_fp_pct']:.0f}:{r['indel_pl_pct']:.0f}:{r['indel_rp_pct']:.0f}"
            decoded = f"{r['file_size']:,}B" if r['file_size'] > 0 else "N/A"
            f.write(f"| {name} | {raw_cov} | {xor_rec} | {final} | {indel} | {decoded} |\n")

        f.write("\n## Region Definitions\n\n")
        f.write("| Region | Position | Description |\n")
        f.write("|--------|----------|-------------|\n")
        f.write("| FP | 1-20 | Forward Primer |\n")
        f.write("| PL | 21-130 | Payload (110bp data region) |\n")
        f.write("| RP | 131-150 | Reverse Primer |\n\n")

        f.write("## Decoded Files\n\n")
        f.write(f"Location: `{os.path.basename(decode_dir)}/`\n\n")
        f.write("| Sample | File | Size | RS Fail | XOR Match |\n")
        f.write("|--------|------|------|---------|----------|\n")

        for r in all_results:
            name = f"{r['sample']}_{r['pair']}" if r['pair'] != 'R1+R2' else r['sample']
            f.write(f"| {name} | {r['output_file'] or 'N/A'} | {r['file_size']:,} | {r['rs_failures']} | {r['xor_match']} |\n")

        # Detailed Statistics
        f.write("\n## Detailed Statistics\n\n")

        # Sequence Counts
        f.write("### Sequence Counts\n\n")
        f.write("| Sample | Total Indices | 150bp Seqs | 149bp Seqs | Matched 150 | Matched 149 |\n")
        f.write("|--------|---------------|------------|------------|-------------|-------------|\n")
        for r in all_results:
            name = f"{r['sample']} ({r['pair']})" if r['pair'] != 'R1+R2' else r['sample']
            f.write(f"| {name} | {r['total_indices']} | {r['seqs_150']:,} | {r['seqs_149']:,} | {r['matched_150']:,} | {r['matched_149']:,} |\n")

        # Depth Statistics
        f.write("\n### Depth Statistics\n\n")
        f.write("| Sample | Avg Depth 150 | Avg Depth 149 | Perfect 150 | Perfect 149 |\n")
        f.write("|--------|---------------|---------------|-------------|-------------|\n")
        for r in all_results:
            name = f"{r['sample']} ({r['pair']})" if r['pair'] != 'R1+R2' else r['sample']
            f.write(f"| {name} | {r['avg_depth_150']:.1f} | {r['avg_depth_149']:.1f} | {r['perfect_150']} | {r['perfect_149']} |\n")

        # XOR Recovery Details
        f.write("\n### XOR Recovery Details\n\n")
        f.write("| Sample | Covered Before | XOR Needed | XOR Success | XOR Fail | Covered After |\n")
        f.write("|--------|----------------|------------|-------------|----------|---------------|\n")
        for r in all_results:
            name = f"{r['sample']} ({r['pair']})" if r['pair'] != 'R1+R2' else r['sample']
            f.write(f"| {name} | {r['raw_coverage']} | {r['xor_needed']} | {r['xor_success']} | {r['xor_fail']} | {r['final_coverage']} |\n")

        # Indel Region Analysis
        f.write("\n### Indel Region Analysis\n\n")
        f.write("| Sample | FP (1-20) | PL (21-130) | RP (131-150) | Total Analyzed |\n")
        f.write("|--------|-----------|-------------|--------------|----------------|\n")
        for r in all_results:
            name = f"{r['sample']} ({r['pair']})" if r['pair'] != 'R1+R2' else r['sample']
            f.write(f"| {name} | {r['indel_fp_pct']:.0f}% | {r['indel_pl_pct']:.0f}% | {r['indel_rp_pct']:.0f}% | {r['indel_total']:,} |\n")

        # Mutation Analysis - Substitutions
        f.write("\n### Mutation Analysis (Substitutions in 150bp)\n\n")
        f.write("| Sample | FP Subs | PL Subs | RP Subs | Total Subs | Analyzed |\n")
        f.write("|--------|---------|---------|---------|------------|----------|\n")
        for r in all_results:
            name = f"{r['sample']} ({r['pair']})" if r['pair'] != 'R1+R2' else r['sample']
            f.write(f"| {name} | {r['sub_fp']} | {r['sub_pl']} | {r['sub_rp']} | {r['sub_total']} | {r['sub_analyzed']:,} |\n")

        # Mutation Analysis - Deletions
        f.write("\n### Mutation Analysis (Deletions in 149bp)\n\n")
        f.write("| Sample | FP Dels | PL Dels | RP Dels | Total Dels | Analyzed |\n")
        f.write("|--------|---------|---------|---------|------------|----------|\n")
        for r in all_results:
            name = f"{r['sample']} ({r['pair']})" if r['pair'] != 'R1+R2' else r['sample']
            f.write(f"| {name} | {r['del_fp']} | {r['del_pl']} | {r['del_rp']} | {r['del_total']} | {r['del_analyzed']:,} |\n")

        # Generated Files
        f.write("\n### Generated Files\n\n")
        f.write("**Per-sample files:**\n")
        f.write("- `[Sample]_index_depth.csv` - Per-index depth data\n")
        f.write("- `[Sample]_mutation_by_position.csv` - Position-wise mutation data\n")
        f.write("- `[Sample]_mutation_histogram.png` - Per-sample mutation histogram\n\n")
        f.write("**Combined histograms:**\n")
        f.write("- `combined_indel_histogram.png` - Indel region distribution\n")
        f.write("- `combined_mutation_histogram.png` - Substitution and deletion analysis\n")
        f.write("- `coverage_comparison_histogram.png` - Raw vs Final coverage\n")

    # Generate combined histograms
    print("\n  Generating combined histograms...")
    generate_combined_histograms(all_results, analysis_dir)

    # 00_summary.csv
    summary_csv = os.path.join(analysis_dir, "00_summary.csv")
    with open(summary_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Sample', 'Pair', 'Total_Indices',
            'Seqs_150', 'Seqs_149', 'Matched_150', 'Matched_149',
            'Avg_Depth_150', 'Avg_Depth_149', 'Perfect_150', 'Perfect_149',
            'Raw_Coverage', 'Raw_Coverage_%', 'Missing_Before_XOR',
            'XOR_Needed', 'XOR_Success', 'XOR_Fail',
            'Final_Coverage', 'Final_Coverage_%',
            'Indel_FP_%', 'Indel_PL_%', 'Indel_RP_%', 'Indel_Total',
            'Sub_FP', 'Sub_PL', 'Sub_RP', 'Sub_Total', 'Sub_Analyzed',
            'Del_FP', 'Del_PL', 'Del_RP', 'Del_Total', 'Del_Analyzed',
            'RS_Failures', 'XOR_Match', 'Output_File', 'File_Size'
        ])
        for r in all_results:
            writer.writerow([
                r['sample'], r['pair'], r['total_indices'],
                r['seqs_150'], r['seqs_149'], r['matched_150'], r['matched_149'],
                f"{r['avg_depth_150']:.2f}", f"{r['avg_depth_149']:.2f}",
                r['perfect_150'], r['perfect_149'],
                r['raw_coverage'], f"{r['raw_coverage_pct']:.2f}", r['missing_before_xor'],
                r['xor_needed'], r['xor_success'], r['xor_fail'],
                r['final_coverage'], f"{r['final_coverage_pct']:.2f}",
                f"{r['indel_fp_pct']:.1f}", f"{r['indel_pl_pct']:.1f}", f"{r['indel_rp_pct']:.1f}",
                r['indel_total'],
                r['sub_fp'], r['sub_pl'], r['sub_rp'], r['sub_total'], r['sub_analyzed'],
                r['del_fp'], r['del_pl'], r['del_rp'], r['del_total'], r['del_analyzed'],
                r['rs_failures'], r['xor_match'], r['output_file'] or '', r['file_size']
            ])

    # README.md
    readme_path = os.path.join(analysis_dir, "README.md")
    with open(readme_path, 'w', encoding='utf-8') as f:
        f.write(f"# DNA Analysis {SCRIPT_VERSION} Output\n\n")
        f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("## Files\n\n")
        f.write("| File | Description |\n")
        f.write("|------|-------------|\n")
        f.write(f"| `SUMMARY_{SCRIPT_VERSION}.md` | Full analysis summary |\n")
        f.write("| `00_summary.csv` | CSV summary data |\n")
        f.write("| `[Sample]_index_depth.csv` | Per-index depth data |\n")
        f.write("| `[Sample]_mutation_by_position.csv` | Position-wise mutation data |\n")
        f.write("| `[Sample]_mutation_histogram.png` | Per-sample mutation histogram |\n")
        f.write("| `combined_*.png` | Combined histograms |\n")
        f.write("| `README.md` | This file |\n\n")
        f.write(f"## Decoded Files\n\n`../{os.path.basename(decode_dir)}/`\n")

    # =========================================================================
    # Console Summary
    # =========================================================================
    print("\n" + "=" * 80)
    print(f"DNA Analysis & Decoding {SCRIPT_VERSION} COMPLETE")
    print("=" * 80)

    print("\n[Sheet 5 Pairs Comparison]")
    print(f"{'Pair':<10} {'Raw Cov':<12} {'Missing':<8} {'XOR Rec':<12} {'Final':<12} {'Indel':<12}")
    print("-" * 70)
    for r in sheet_results:
        print(f"{r['pair']:<10} {r['raw_coverage_pct']:.1f}%{'':<6} {r['missing_before_xor']:<8} "
              f"{r['xor_success']}/{r['xor_needed']:<8} {r['final_coverage_pct']:.1f}%{'':<6} "
              f"{r['indel_fp_pct']:.0f}:{r['indel_pl_pct']:.0f}:{r['indel_rp_pct']:.0f}")

    print(f"\n[All Samples]")
    print(f"{'Sample':<12} {'Raw Cov':<10} {'Final':<10} {'Indel':<12} {'Decoded':<15}")
    print("-" * 60)
    for r in all_results:
        if r['sample'] == 'Sheet' and r['pair'] != 'R1+R2':
            continue
        print(f"{r['sample']:<12} {r['raw_coverage_pct']:.1f}%{'':<5} {r['final_coverage_pct']:.1f}%{'':<5} "
              f"{r['indel_fp_pct']:.0f}:{r['indel_pl_pct']:.0f}:{r['indel_rp_pct']:.0f}{'':<4} "
              f"{r['file_size']:,}B")

    print(f"\nAnalysis: {analysis_dir}")
    print(f"Decoded: {decode_dir}")


if __name__ == "__main__":
    main()
