#!/usr/bin/env python3
# coding=utf-8
# title           :oligomutk.py
# description     :Covid primer monitor / Mutant Screening
# date            :20220317
# version         :1.2.0
# copyright       :Michael Bekaert
# notes           :Needs KAT, emboss (water) and Biopython
# ==============================================================================

import argparse
import gzip
import json
import os
import subprocess
import sys
import tempfile
from functools import partial

from Bio import SeqIO
from joblib import Parallel, delayed

release = '1.2.0'


def seq_diff(seq1, seq2):
    ret = ''
    diff = 0
    if len(seq1) == len(seq2):
        for i in range(0, len(seq1)):
            if seq1[i] == seq2[i]:
                ret += '.'
            else:
                diff += 1
                ret += seq2[i].lower() if seq1[i] != '-' else seq2[i].upper()

    return [ret, diff]


def collect_new(file, forward, reverse, path='/tmp/seq', trust=False, verbose=False):
    print('Pre-processing the sequences', file=sys.stdout)
    forward_len = len(forward)
    reverse_len = len(reverse)
    with open(path + '_forward.fa', 'w') as outfile:
        outfile.write('>forward\n' + str(forward).upper() + '\n')
    with open(path + '_reverse.fa', 'w') as outfile:
        outfile.write('>reverse\n' + str(reverse).upper() + '\n')

    if trust:
        if file.endswith('.gz'):
            subprocess.call('gunzip -c ' + file + ' >' + path + '_genomes.fa', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        else:
            os.symlink(file, path + '_genomes.fa')
    else:
        # filter the input fasta file for wrong, duplicated or misspell records
        _open = partial(gzip.open, mode='rt') if file.endswith('.gz') else open

        # check duplicated names and fix names
        names = []
        with open(path + '_genomes.fa', 'w') as genome, _open(file) as f:
            for record in SeqIO.parse(f, 'fasta'):
                name = record.description.replace(' ', '_').replace('\t', '_')
                if name not in names:
                    names.append(name)
                    genome.write('>' + name + '\n' + str(record.seq).upper().replace('U', 'T') + '\n')
                else:
                    print('"' + name + '" already exists!', file=sys.stderr)

    if os.path.getsize(path + '_genomes.fa') > 32 and forward_len > 0 and reverse_len > 0:
        return [path + '_genomes.fa', path + '_forward.fa', forward_len, path + '_reverse.fa', reverse_len]

    return False


def run_step1(file, forward, forward_len, reverse, reverse_len, path='/tmp/seq', cpu=1, verbose=False):
    print('Processing the sequences with KAT', file=sys.stdout)
    genomes = {}

    cmd = 'kat sect -t ' + str(cpu) + ' -m ' + str(forward_len) + ' -o ' + path + '_forward ' + str(file) + ' ' + str(
        forward)
    if verbose:
        print('... ' + cmd)
    subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if os.path.exists(path + '_forward-stats.tsv'):
        with open(path + '_forward-stats.tsv', 'r') as fwin:
            for line in fwin:
                if not line.startswith('seq_name'):
                    tab = line.split('\t')
                    if len(tab) >= 9 and len(tab[0]) > 0 and len(tab[8]) > 0:
                        if int(tab[8]) > 0:
                            genomes[tab[0]] = {'forward': 1, 'forward_diff': '.' * forward_len}  # , 'forward_var': 0, 'forward_len': forward_len, 'forward_line': line}
                        else:
                            genomes[tab[0]] = {}  # 'forward_line': line}
        os.remove(path + '_forward-stats.tsv')
        os.remove(path + '_forward-counts.cvg')

    cmd = 'kat sect -t ' + str(cpu) + ' -m ' + str(reverse_len) + ' -o ' + path + '_reverse ' + str(file) + ' ' + str(reverse)
    if verbose:
        print('... ' + cmd)
    subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if os.path.exists(path + '_reverse-stats.tsv'):
        with open(path + '_reverse-stats.tsv', 'r') as rvin:
            for line in rvin:
                if not line.startswith('seq_name'):
                    tab = line.split('\t')
                    if len(tab) >= 9 and len(tab[0]) > 0 and len(tab[8]) > 0:
                        if int(tab[8]) > 0:
                            genomes[tab[0]].update({'reverse': 1, 'reverse_diff': '.' * reverse_len})  # , 'reverse_var': 0, 'reverse_len': reverse_len, 'reverse_line': line})
                            # else:
                            #    genomes[tab[0]].update({'reverse_line': line})
        os.remove(path + '_reverse-stats.tsv')
        os.remove(path + '_reverse-counts.cvg')

    return genomes


def run_step2(genomes, file, forward, reverse, probe=None, amplicon=None, all=False, output=None, path='/tmp/seq', max_diff=5, cpu=1, verbose=False):
    def get_alignment(record):
        if verbose:
            print('  - Aligning: ' + str(record.id), file=sys.stdout)
        log_issue = 0

        try:
            tmp_name = next(tempfile._get_candidate_names())
            with open(path + tmp_name + '_local.fa', 'w') as outaln:
                outaln.write('>local\n' + str(record.seq) + '\n')

            cmd = 'water -gapopen 10 -gapextend 0.5 -nobrief --auto -aformat3 fasta -outfile ' + path + tmp_name + '_local.aln.fa ' + path + tmp_name + '_local.fa ' + path + '_both.fa'
            if verbose:
                print('... ' + cmd)

            subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)

            local = None
            if os.path.exists(path + tmp_name + '_local.aln.fa'):
                for record2 in SeqIO.parse(path + tmp_name + '_local.aln.fa', 'fasta'):
                    if record2.id == 'local':
                        local = str(record2.seq)
                    elif record2.id != 'local' and local is not None:
                        (seq, variations) = seq_diff(str(record2.seq), local)
                        if len(seq) > 1:
                            if variations > max_diff or abs(len(str(record2.seq)) - len(local)) > max_diff:
                                log_issue += 1
                                genomes[record.id]['issue'] = log_issue
                            genomes[record.id][record2.id.lower() + '_diff'] = seq
                            # genomes[record.id][record2.id.lower() + '_var'] = variations
                            # genomes[record.id][record2.id.lower() + '_len'] = len(str(record2.seq))
                        local = None
                    else:
                        print(
                            'Warning: "water" has generated a critical error for ' + record.id + ' sequence will be skipped.',
                            file=sys.stderr)
            if output is not None and record.id in check_sequence:
                with open(output + ('.issue' if log_issue > 0 else '.fail') + '.fasta', 'a') as outfile:
                    outfile.write('>' + str(record.id) + '\n' + str(record.seq) + '\n')

            if os.path.exists(path + tmp_name + '_local.aln.fa'):
                os.remove(path + tmp_name + '_local.aln.fa')
            if os.path.exists(path + tmp_name + '_local.fa'):
                os.remove(path + tmp_name + '_local.fa')
        except:
            print('Error: "water" has generated a critical error for ' + record.id + ' sequence will be skipped.',
                  file=sys.stderr)

    print('Checking the intriguing sequences with local alignments', file=sys.stdout)

    check_sequence = set()

    if len(genomes) > 0:
        for genome in genomes.keys():
            flag = set()
            if not ('forward' in genomes[genome] and 'reverse' in genomes[genome]):
                check_sequence.add(str(genome))

    if len(check_sequence) > 0:
        subprocess.call('cat ' + forward + ' ' + reverse + '> ' + path + '_both.fa', stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)

        if (isinstance(probe, list) and len(probe) > 0) or amplicon is not None:
            with open(path + '_both.fa', 'a') as outfile:
                if isinstance(probe, list) and len(probe) > 0:
                    probes = 0
                    for item in probe:
                        probes += 1
                        outfile.write('>probe' + str(probes) + '\n' + str(item).upper() + '\n')
                if amplicon is not None:
                    outfile.write('>amplicon\n' + str(amplicon).upper() + '\n')

        Parallel(n_jobs=cpu, require='sharedmem')(delayed(get_alignment)(record) for record in SeqIO.parse(file, 'fasta') if all or record.id in check_sequence)

        if os.path.exists(path + '_both.fa'):
            os.remove(path + '_both.fa')

    return genomes


def report(genomes, probe=None, amplicon=None):
    print('Reporting', file=sys.stdout)

    log_pass = 0
    log_fail = 0
    log_issue = 0
    log = []
    if len(genomes) > 0:
        for genome in genomes.keys():
            flag = set()
            if 'issue' in genomes[genome]:
                log_issue += 1
                flag.add('Sequencing accuracy to be checked manually')
            elif 'forward' in genomes[genome] and 'reverse' in genomes[genome]:
                log_pass += 1
            else:
                log_fail += 1
                if 'forward' not in genomes[genome]:
                    if 'forward_diff' in genomes[genome]:
                        flag.add('forward primer region has a sequence variation')
                    else:
                        flag.add('forward primer was not found')
                elif 'reverse' not in genomes[genome]:
                    if 'reverse_diff' in genomes[genome]:
                        flag.add('reverse primer region has a sequence variation')
                    else:
                        flag.add('reverse primer was not found')
            tmp = ''
            if isinstance(probe, list) and len(probe) > 0:
                probes = 0
                for item in probe:
                    probes += 1
                    if 'probe' + str(probes) + '_diff' in genomes[genome]:
                        tmp += '\t' + genomes[genome]['probe' + str(probes) + '_diff']
                    else:
                        tmp += '\t'

            tmp += (('\t' + (genomes[genome]['amplicon_diff'] if 'amplicon_diff' in genomes[genome] else '')) if amplicon is not None else '')

            log.append(str(genome) + '\t' + ('PASS' if len(flag) == 0 else ('ISSUE' if log_issue > 0 else 'FAIL')) + '\t' + (genomes[genome]['forward_diff'] if 'forward_diff' in genomes[genome] else '') + '\t' + (genomes[genome]['reverse_diff'] if 'reverse_diff' in genomes[genome] else '') + tmp + '\t' + '|'.join(flag))
    return [log_pass, log_fail, log_issue, log]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fasta', '-i', '-in', dest='genomesfile', type=str, required=True, help='Specify covid fasta genomes file')
    parser.add_argument('-forward', '-1', dest='forward', type=str, required=True, help='Specify forward primer')
    parser.add_argument('-reverse', '-2', dest='reverse', type=str, required=True, help='Specify reverse primer')
    parser.add_argument('-probe', dest='probe', nargs='+', type=str, help='Specify probe sequence(s)')
    parser.add_argument('-amplicon', dest='amplicon', type=str, help='Specify amplicon sequence')
    parser.add_argument('-tmp', dest='tmp', type=str, default=tempfile._get_default_tempdir(), help='temporary folder location')
    parser.add_argument('-max', '-m', dest='max_diff', type=int, default=5, help='Maximum sequence variations before querying sequencing accuracy')
    parser.add_argument('-output', '-out', '-o', dest='outfile', type=str, default='output', help='Specify output file prefix')
    parser.add_argument('-threads', '-t', dest='threads', type=int, default=1, help='Maximum threads')
    parser.add_argument('-all', '-a', dest='all', action='store_true', default=False, help='Test all sequences (Slow)')
    parser.add_argument('-sort', '-s', dest='run_sort', action='store_true', default=False, help='Run "sortseq.py" script')
    parser.add_argument('-json', dest='json', action='store_true', default=False, help='Export results matrix as JSON file (Experimental)')
    parser.add_argument('-trust', dest='trust', action='store_true', default=False, help='Trust input file and don\'t run sanity checks')
    parser.add_argument('-verbose', '-v', dest='verbose', action='store_true', default=False, help='Become very chatty (Development)')

    args = parser.parse_args()

    if os.path.exists(args.genomesfile):
        path = os.path.abspath(os.path.join(args.tmp, next(tempfile._get_candidate_names())))
        ret = collect_new(os.path.abspath(args.genomesfile), args.forward, args.reverse, path=path, trust=args.trust, verbose=args.verbose)
        if ret:
            (file, forward, forward_len, reverse, reverse_len) = ret
            genomes = run_step1(file, forward, forward_len, reverse, reverse_len, cpu=args.threads, path=path, verbose=args.verbose)
            genomes = run_step2(genomes, file, forward, reverse, probe=args.probe, amplicon=args.amplicon, all=args.all, output=args.outfile, max_diff=args.max_diff, cpu=args.threads, path=path, verbose=args.verbose)

            (log_pass, log_fail, log_issue, log) = report(genomes, probe=args.probe, amplicon=args.amplicon)
            print('\nANALYSED: ' + str(log_pass + log_fail) + '\n PASS:    ' + str(log_pass) + '\n FAIL:    ' + str(log_fail) + ('\n QUERY:   ' + str(log_issue) if log_issue > 0 else ''), file=sys.stdout)

            tmp = ''
            if isinstance(probe, list) and len(probe) > 0:
                probes = 0
                for item in probe:
                    probes += 1
                    tmp += 'PROBE' + str(probes) + '_DIFF\t'

            with open(args.outfile + '.tsv', 'w') as out:
                out.write('GENOME_ID\tPASS/FAIL\tFORWARD_DIFF\tREVERSE_DIFF\t' + tmp + ('AMPLICON_DIFF\t' if args.amplicon is not None else '') + 'NOTES\n')
                out.write('\n'.join(log) + '\n')

            if args.json:
                with open(args.outfile + '.json', 'w') as outfile:
                    json.dump(genomes, outfile)

            if os.path.exists(forward):
                os.remove(forward)
            if os.path.exists(reverse):
                os.remove(reverse)
            if os.path.exists(file):
                os.remove(file)
