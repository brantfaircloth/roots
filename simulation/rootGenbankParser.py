#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Brant Faircloth on 2009-10-28.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import pdb
import re
from Bio import SeqIO


def main():
    matK_records = []
    for record in SeqIO.parse(open('Kressetal_all_root_barcodes.fasta', 'rU'), 'fasta'):
        regex = re.compile('rbcL')
        search = re.search(regex, record.description)
        if search:
            matK_records.append(record)
    out = open('/Users/bcf/Documents/UCLA/Lab/Organisms/Roots/Kressetal_rbcL_records.fa', 'w')
    SeqIO.write(matK_records, out, 'fasta')
    out.close()

if __name__ == '__main__':
    main()

