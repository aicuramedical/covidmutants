#!/usr/bin/env python3
# coding=utf-8
# title           :typetempsort.py
# description     :Monthly parsing according to variants (very slow)
# date            :20211214
# version         :1.1.0
# copyright       :Michael Bekaert
# ==============================================================================

import argparse
import codecs
import datetime
import os.path
import tarfile

from Bio import SeqIO

release = '1.1.0'

date_time_obj = '%Y-%m-%d'
variants = []

parser = argparse.ArgumentParser(description='Parse GISAID bulk files into monthly parcels.')
parser.add_argument('-m', '--meta', dest='metapath', type=str, required=True,
                    help='File path to the GISAID "metadata_tsv_20xx_xx_xx.tar.xz" file')
parser.add_argument('-f', '--fasta', dest='sequencepath', type=str, required=True,
                    help='File path to the GISAID "sequences_fasta_20xx_xx_xx.tar.xz" file')
parser.add_argument('-s', '--start', dest='start_date', type=str, required=True,
                    help='Start date (inclusive): eg, 2020-05 (for May 2020)')
parser.add_argument('-e', '--end', dest='end_date', type=str, required=True,
                    help='End date (inclusive): eg, 2021-05 (for May 2021)')
parser.add_argument('-v', '--virus', dest='virus',
                    help='Virus subtype to be extracted: eg, Mu, Alpha, Beta, Gamma, Delta, GH/490R, Lambda, Omicron or ALL')
args = parser.parse_args()

if not os.path.exists(args.metapath) or not os.path.exists(args.sequencepath):
    parser.print_help()
    exit(1)

try:
    starts = datetime.datetime.strptime(args.start_date, '%Y-%m')
    ends = datetime.datetime.strptime(args.end_date, '%Y-%m')
except:
    parser.print_help()
    exit(2)

if args.virus == 'ALL':
    args.virus = None
    variants = ['Mu', 'Alpha', 'Beta', 'Gamma', 'Delta', 'GH/490R', 'Lambda', 'Omicron']  # True as 2021-12-09

print('Pre-select sequences...', end='\r')
keep = {}
data_file = set()
nb_sequences = 0
with tarfile.open(args.metapath, "r:xz") as tar:
    for member in tar.getmembers():
        if member.name == 'metadata.tsv':
            for sequence in tar.extractfile(member).readlines():
                dat = sequence.decode().split('\t')
                if len(dat) > 14 and dat[3].startswith('20'):
                    nb_sequences += 1
                    try:
                        timestamp = datetime.datetime.strptime(dat[3], date_time_obj)
                        if (args.virus is None or args.virus in dat[13]) and timestamp <= ends and timestamp >= starts:
                            if len(variants) > 0:
                                for ext in variants:
                                    if ext in dat[13]:
                                        keep[dat[0]] = [ext.replace('/', '_') + '_' + timestamp.strftime('%Y-%m'),
                                                        '[' + str(dat[10]) + '] [' + str(dat[11]) + '] [' + str(
                                                            dat[13]) + ']']
                                        data_file.add(ext.replace('/', '_') + '_' + timestamp.strftime('%Y-%m'))
                            else:
                                keep[dat[0]] = [timestamp.strftime('%Y-%m'),
                                                '[' + str(dat[10]) + '] [' + str(dat[11]) + '] [' + str(dat[13]) + ']']
                                data_file.add(timestamp.strftime('%Y-%m'))
                    except:
                        pass

print('Pre-select sequences... found ' + str(len(keep.keys())) + ' out of ' + str(nb_sequences) + ' sequences')

if len(keep.keys()) > 0:
    print('Accessing sequences...', end='\r')
    nb_fastas = 0
    loop = round(nb_sequences / 100)

    with tarfile.open(args.sequencepath, "r:xz") as tar:
        for member in tar.getmembers():
            if member.name == 'sequences.fasta':

                # open destination files
                file_handles = {}
                for file in data_file:
                    file_handles[file] = open(file + '.fasta', 'w')

                print('Accessing sequences (that may take a while)... 0%', end='\r')

                with tar.extractfile(member) as f:
                    ft = codecs.getreader("utf-8")(f)

                    for record in SeqIO.parse(ft, "fasta"):
                        nb_fastas += 1
                        if nb_fastas % loop == 0:
                            print('Accessing sequences (that may take a while)... ' + str(
                                round(nb_fastas / nb_sequences * 100)) + '%', end='\r')
                        name = record.id.split('|')
                        if name[0] in keep.keys():
                            file_handles[keep[name[0]][0]].write(
                                '>' + str(record.id) + ' ' + str(keep[name[0]][1]) + '\n' + str(record.seq) + '\n')

                for file in file_handles:
                    file_handles[file].close()

print('Accessing sequences (that may take a while)... done')

# Get metadata_tsv_20xx_xx_xx.tar.xz
# Get sequence_fasta_20xx_xx_xx.tar.xz
# Run the script
#    typetempsort.py --meta metadata_tsv_20xx_xx_xx.tar.xz --fasta sequences_fasta_20xx_xx_xx.tar.xz --start 2021-10 --end 2021-10 --virus Delta
