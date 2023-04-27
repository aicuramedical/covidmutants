#!/usr/bin/env python3
# coding=utf-8
# title           :typetempsort.py
# description     :Monthly parsing according to variants (very slow)
# date            :20230427
# version         :1.3.2
# copyright       :Micha\"el Bekaert
# ==============================================================================

import argparse
import codecs
import datetime
import gzip
import os.path
import tarfile

from Bio import SeqIO

release = '1.3.2'

date_time_obj = '%Y-%m-%d'
variants = []

parser = argparse.ArgumentParser(description='Parse GISAID bulk files into monthly parcels.')
parser.add_argument('-i', '--input', dest='inputpath', type=str, required=True,
                    help='File path to the GISAID "Downloads > Genomic epidemilogy - nextregions" export file "hcov_xxxxxx_xxxx-xx-xx_xx-xx.tar.gz" file')
parser.add_argument('-s', '--start', dest='start_date', type=str, required=True,
                    help='Start date (inclusive): eg, 2020-05 (for May 2020)')
parser.add_argument('-e', '--end', dest='end_date', type=str, required=True,
                    help='End date (inclusive): eg, 2023-04 (for April 2023)')
parser.add_argument('-v', '--virus', dest='virus',
                    help='Covid19 virus subtype to be extracted: eg, Alpha, Beta, Gamma, Delta, GH/490R, Epsilon, Eta, Iota,Kappa, Lambda, Mu, Theta, Zeta, Omicron or ALL')
args = parser.parse_args()

if not os.path.exists(args.inputpath):
    parser.print_help()
    exit(1)

try:
    starts = datetime.datetime.strptime(args.start_date, '%Y-%m')
    ends = datetime.datetime.strptime(args.end_date, '%Y-%m')
except:
    parser.print_help()
    exit(2)

if args.virus is None or args.virus.upper() == 'ALL':
    args.virus = None
    # variants = ['ALPHA', 'BETA', 'DELTA', 'GAMMA', 'OMICRON', 'EPSILON', 'ETA', 'IOTA', 'KAPPA', 'LAMBDA', 'MU', 'THETA', ' ZETA', 'GH/490R']  # True as of 2022-10-20
    variants = ['ALPHA', 'BETA', 'DELTA', 'GAMMA', 'OMICRON', 'LAMBDA', 'MU', 'GH/490R']  # True as of 2023-04-27
else:
    variants.append(args.virus.upper())

keep = {}
data_file = set()
nb_sequences = 0
with tarfile.open(args.inputpath, "r:gz") as tar:
    print('Pre-select sequences...', end='\r')
    for member in tar.getmembers():
        if member.name.endswith('.tsv'):
            meta_header = {'strain': -1, 'gisaid_epi_isl': -1, 'date': -1, 'variant': -1}
            for sequence in tar.extractfile(member).readlines():
                dat = sequence.decode().strip().split('\t')
                if len(dat) >= 10 and meta_header['gisaid_epi_isl'] == -1 and 'gisaid_epi_isl' in dat:
                    for item_num, item in enumerate(dat):
                        for key in meta_header.keys():
                            if item.lower() == key:
                                meta_header[key] = item_num
                elif len(dat) >= 10 and meta_header['gisaid_epi_isl'] != -1:
                    nb_sequences += 1
                    try:
                        timestamp = datetime.datetime.strptime(dat[meta_header['date']], date_time_obj)
                        if timestamp <= ends and timestamp >= starts:
                            for ext in variants:
                                if ' ' + ext in dat[meta_header['variant']].upper():
                                    keep[dat[meta_header['strain']]] = [ext.replace('/', '_') + '_' + timestamp.strftime('%Y-%m'),
                                                                        '[' + str(dat[meta_header['gisaid_epi_isl']]) + '] [' + str(dat[meta_header['date']]) + ']']
                                    data_file.add(ext.replace('/', '_') + '_' + timestamp.strftime('%Y-%m'))
                    except:
                        pass
            print('Pre-select sequences... found ' + str(len(keep.keys())) + ' out of ' + str(nb_sequences) + ' sequences')

    if len(keep.keys()) > 0:
        print('Accessing sequences...', end='\r')
        nb_fastas = 0
        loop = round(nb_sequences / 100)

        for member in tar.getmembers():
            if member.name.endswith('.fasta'):

                # open destination files
                file_handles = {}
                for file in data_file:
                    file_handles[file] = gzip.open(file + '.fasta.gz', 'wt')

                print('Accessing sequences (that may take a while)... 0%', end='\r')

                with tar.extractfile(member) as f:
                    ft = codecs.getreader('utf-8')(f)

                    for record in SeqIO.parse(ft, 'fasta'):
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

# Get hcov_xxxxxx_xxxx-xx-xx_xx-xx.tar.gz
# Run the script
#    typetempsort.py --input hcov_xxxxxx_xxxx-xx-xx_xx-xx.tar.gz --start 2020-05 --end 2023-04 --virus Delta
