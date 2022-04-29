#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import csv

class FastaFilterError(Exception): pass

#-f /data/mus_musculus/GCF_000001635.26_GRCm38.p6_genomic.fna -i /tools/HybridScaffold/scripts/data/mm_fasta_filter_chrs.txt -o /data/mus_musculus/GCF_000001635.26_GRCm38.p6_genomic.filtered.fna


def run(fasta, include_chrs, output):
    chrs_to_keep = parse_include_chrs_file(include_chrs)
    in_chr_to_keep = False
    with open(output, 'w') as outfile:
        with open(fasta, 'r') as infile:
            for line in infile:
                line = line.strip()
                is_defline = check_if_defline(line)
                if is_defline:
                    in_chr_to_keep = check_if_defline_indicates_chr_to_keep(line, chrs_to_keep)
                if in_chr_to_keep:
                    write_line_to_output(line, is_defline, chrs_to_keep, outfile)


def parse_include_chrs_file(include_chrs):
    chrs_to_keep = {}
    with open(include_chrs, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        next(csvreader)
        for row in csvreader:
            chrs_to_keep[row[0]] = row[1]
    return chrs_to_keep


def check_if_defline(line):
    return line.startswith('>')


def check_if_defline_indicates_chr_to_keep(line, chrs_to_keep):
    return line.strip('>') in chrs_to_keep


def write_line_to_output(line, is_defline, chrs_to_keep, outfile):
    if is_defline:
        str_to_write = '>' + chrs_to_keep[line.strip('>')]
    else:
        str_to_write = line
    outfile.write(str_to_write + '\n')


def main():
    # command line arguments
    parser = argparse.ArgumentParser(
        description='Filter FASTA file to only the chromosomes/contigs listed in the include_chrs file, and rename chromosomes if needed.',
        epilog='filter_fasta version 1.0?1 (c)2020 Bionano Genomics Inc. all rights reserved. Written by Aliz R Rao')
    parser.add_argument('--fasta', '-f', required=True,
                        help='''Input FASTA file to filter.''')
    parser.add_argument('--include_chrs', '-i', required=True,
                        help='''Tab-delimited text file containing columns CompntName and RenameTo. The column CompntName should contain the names of the contigs in the FASTA file that you wish to keep; these can be copied easily from the key file generated during in-silico digestion of the original (unfiltered) FASTA. Contigs missing from this column will be filtered out. Contigs will be renamed as listed in the RenameTo column.''')
    parser.add_argument('--output', '-o', required=True,
                        help='''Path and filename for filtered FASTA output file.''')

    args = parser.parse_args()

    run(fasta=args.fasta, include_chrs=args.include_chrs, output=args.output)

    return 0

if __name__ == "__main__": sys.exit(main())
