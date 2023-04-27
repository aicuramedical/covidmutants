#!/usr/bin/env python3
# coding=utf-8
# title           :sortseq.py
# description     :Counting of mutants
# date            :20230427
# version         :1.3.2
# copyright       :Elena Graf
# notes           :Needs numpy and pandas.
# ==============================================================================

import argparse
import os
import re

import numpy as np
import pandas as pd

release = '1.3.2'


def mutsec(filepath: str, output: str):
    # allowed length of FP and RP respectively
    fp_length = 32
    rp_length = 32

    # How many variations are allowed?
    allowed_characters = 3

    # Import tsv Data
    df = pd.read_csv(filepath, sep='\t')

    # Calculations
    # Only FAIL
    df_only_fail = df.drop(df[df['PASS/FAIL'] == 'PASS'].index)
    df_only_fail.index = range(len(df_only_fail))

    # Exclude Data where FORWARD_DIFF or REVERSE_DIFF are too short
    # drop rows where FP not complete
    incomplete_fp = []
    for i in range(len(df_only_fail)):
        if len(df_only_fail.FORWARD_DIFF[i]) < fp_length:
            incomplete_fp.append(i)

    df_full_fp = df_only_fail.drop(incomplete_fp)
    df_full_fp.index = range(len(df_full_fp))

    # drop rows where RP not complete
    incomplete_rp = []
    for i in range(len(df_full_fp)):
        if len(df_full_fp.REVERSE_DIFF[i]) < rp_length:
            incomplete_rp.append(i)

    df_full_primers = df_full_fp.drop(incomplete_rp)
    df_full_primers.index = range(len(df_full_primers))

    # Exclude Data with too many changes
    # Filter FORWARD_DIFF
    indices_too_many_changes_fp = []  # indices where REVERSE_DIFF has too many changes
    for i in range(len(df_full_primers)):
        forward_diff = df_full_primers.FORWARD_DIFF[i]
        if forward_diff.count('.') < len(forward_diff) - allowed_characters:  # check number of .
            indices_too_many_changes_fp.append(i)

    df_allowed_fp_changes = df_full_primers.drop(indices_too_many_changes_fp)
    df_allowed_fp_changes.index = range(len(df_allowed_fp_changes))

    # Filter REVERSE_DIFF
    indices_too_many_changes_rp = []  # indices where REVERSE_DIFF has too many changes
    for i in range(len(df_allowed_fp_changes)):
        reverse_diff = df_allowed_fp_changes.REVERSE_DIFF[i]
        if reverse_diff.count('.') < len(reverse_diff) - allowed_characters:  # check number of .
            indices_too_many_changes_rp.append(i)

    df_allowed_primer_changes = df_allowed_fp_changes.drop(indices_too_many_changes_rp)
    df_allowed_primer_changes.index = range(len(df_allowed_primer_changes))

    # Check for IUD codes or -
    not_allowed = "[-, N, n, V, v, B, b, H, h, D, d, K, k, S, s, W, w, M, m, Y, y, R, r]"
    # Filter through FORWARD_DIFF
    not_allowed_fp = []
    for i in range(len(df_allowed_primer_changes)):
        forward_diff = df_allowed_primer_changes.FORWARD_DIFF[i]
        if len(re.findall(not_allowed, forward_diff)) > 0:
            not_allowed_fp.append(i)

    df_allowed_fp = df_allowed_primer_changes.drop(not_allowed_fp)
    df_allowed_fp.index = range(len(df_allowed_fp))

    # Filter through REVERSE_DIFF
    not_allowed_rp = []
    for i in range(len(df_allowed_fp)):
        reverse_diff = df_allowed_fp.REVERSE_DIFF[i]
        if len(re.findall(not_allowed, reverse_diff)) > 0:
            not_allowed_rp.append(i)

    df_allowed_primer = df_allowed_fp.drop(not_allowed_rp)
    df_allowed_primer.index = range(len(df_allowed_primer))

    # Find Pairs and count occurrences
    # Find Pairs -> create arrays with them -> delete them from df -> append pairs dataframe to other one
    forward_diff_pairs = []
    reverse_diff_pairs = []
    indices_pairs = []
    for i in range(len(df_allowed_primer)):
        forward_diff = df_allowed_primer.FORWARD_DIFF[i]
        reverse_diff = df_allowed_primer.REVERSE_DIFF[i]
        if forward_diff.count('.') < len(forward_diff) and reverse_diff.count('.') < len(reverse_diff):
            forward_diff_pairs.append(forward_diff)
            reverse_diff_pairs.append(reverse_diff)
            indices_pairs.append(i)

    data_pairs = {
        'FP variation': forward_diff_pairs,
        'RP variation': reverse_diff_pairs
    }

    df_pairs = pd.DataFrame(data=data_pairs)
    df_pairs_counter = df_pairs.value_counts().reset_index(name='variation_counts')[
        ['variation_counts', 'FP variation', 'RP variation']].sort_values(by=['FP variation', 'RP variation'],
                                                                          ascending=False)
    df_pairs_counter.index = range(len(df_pairs_counter))

    # Delete Rows with Pairs from df_allowed_primer
    df_allowed_primer = df_allowed_primer.drop(indices_pairs)
    df_allowed_primer.index = range(len(df_allowed_primer))

    # Find and count occurrences in either FORWARD_DIFF or REVERSE_DIFF

    for i in range(len(df_allowed_primer)):
        if len(df_allowed_primer.REVERSE_DIFF[i]) != 32:
            print(len(df_allowed_primer.REVERSE_DIFF[i]))

    df_allowed_primer_sorted = df_allowed_primer.sort_values(by=['FORWARD_DIFF', 'REVERSE_DIFF'], ascending=False)
    df_allowed_primer_sorted.index = range(len(df_allowed_primer_sorted))

    count_forward_primers = np.array(np.unique(df_allowed_primer_sorted.FORWARD_DIFF, return_counts=True))
    count_reverse_primers = np.unique(df_allowed_primer_sorted.REVERSE_DIFF, return_counts=True)

    # Filter out variations with no change
    # Forward
    not_empty_forward_primer = []
    not_empty_forward_primer_counter = []
    counter_forward = 0
    for i in range(len(count_forward_primers[0])):  # avoid error
        if count_forward_primers[0][i].count('.') == len(count_forward_primers[0][i]):
            not_empty_forward_primer = np.delete(count_forward_primers[0], i)
            not_empty_forward_primer_counter = np.delete(count_forward_primers[1], i)
            counter_forward += 1

    if counter_forward == 0:
        not_empty_forward_primer = count_forward_primers[0]
        not_empty_forward_primer_counter = count_forward_primers[1]

    # Reverse
    not_empty_reverse_primer = []
    not_empty_reverse_primer_counter = []
    counter_reverse = 0
    for i in range(len(count_reverse_primers[0])):  # avoid error
        if count_reverse_primers[0][i].count('.') == len(count_reverse_primers[0][i]):
            not_empty_reverse_primer = np.delete(count_reverse_primers[0], i)
            not_empty_reverse_primer_counter = np.delete(count_reverse_primers[1], i)
            counter_reverse += 1

    if counter_reverse == 0:
        not_empty_reverse_primer = count_reverse_primers[0]
        not_empty_reverse_primer_counter = count_reverse_primers[1]

    not_empty_forward_primer = np.flip(not_empty_forward_primer)
    not_empty_forward_primer_counter = np.flip(not_empty_forward_primer_counter)

    not_empty_reverse_primer = np.flip(not_empty_reverse_primer)
    not_empty_reverse_primer_counter = np.flip(not_empty_reverse_primer_counter)

    # Create DataFrame with Data

    type_of_change = []
    forward_primer_change_total = []
    reverse_primer_change_total = []
    change_counts = []

    for i in range(len(not_empty_forward_primer) + len(not_empty_reverse_primer)):
        if i in range(len(not_empty_forward_primer)):
            type_of_change.append('FORWARD DIFF')
            forward_primer_change_total.append(not_empty_forward_primer[i])
            reverse_primer_change_total.append("")
            change_counts.append(not_empty_forward_primer_counter[i])
        else:
            type_of_change.append('REVERSE DIFF')
            forward_primer_change_total.append("")
            reverse_primer_change_total.append(not_empty_reverse_primer[i - len(not_empty_forward_primer)])
            change_counts.append(not_empty_reverse_primer_counter[i - len(not_empty_forward_primer)])

    d = {
        'variation_counts': change_counts,
        'FP variation': forward_primer_change_total,
        'RP variation': reverse_primer_change_total,
    }

    df_final = pd.DataFrame(data=d)
    df_final_with_pairs = df_pairs_counter._append(df_final, ignore_index=True)

    # Export to Excel Table
    df_final.to_csv(output, sep='\t')

    # General Calculations
    number_of_pass = len(df[df['PASS/FAIL'] == 'PASS'])
    number_of_accepted_variations = sum(np.array(df_final_with_pairs['variation_counts']))

    print('Total amount of data: ' + str(len(df)))
    print('Number of PASS: ' + str(number_of_pass))
    print('Number of accepted variations: ' + str(number_of_accepted_variations))
    print('Number of not accepted variations: ' + str(len(df) - number_of_pass - number_of_accepted_variations))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarise results.')
    parser.add_argument('-i', '--in', dest='infile', type=str, required=True,
                        help='Input file (oligomutk.py output) ')
    parser.add_argument('-o', '--out', dest='outfile', type=str, default='output.tsv',
                        help='Result/summary file (TSV formatted)')
    args = parser.parse_args()

    if os.path.exists(args.infile):
        mutsec(os.path.abspath(args.infile), os.path.abspath(args.outfile))

# Run oligomutk
# Run the script
#    sortseq.py --in oligomutk_output.tsv --out results.tsv
