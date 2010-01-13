#!/usr/bin/env python
# encoding: utf-8
"""
seq_matcher.py

Created by Brant Faircloth on 2009-12-12.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import pdb
import os
import re
import glob
import numpy
import optparse
import cStringIO
import subprocess
import multiprocessing
from Queue import Empty
from Bio import SeqIO
from Bio import pairwise2
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIStandalone
# make sure we use sqlite supporting foreign keys
from pysqlite2 import dbapi2 as sqlite3


# python ../matcher.py -i ./ -d psbA-trnH.sqlite -e /usr/local/ncbi/blast/bin/blastn -b ../blast/voucher_trnh_rbcl -l 100 -n 6


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"
    p = optparse.OptionParser(usage)
    p.add_option('--input', '-i', dest = 'input', action='store', \
    type='string', default = None, \
    help='Path to the sequences taken from each core', \
    metavar='FILE')
    p.add_option('--database', '-d', dest = 'database', action='store', \
    type='string', default = None, \
    help='Path to the database holding simulated sequence info', \
    metavar='FILE')
    p.add_option('--executable', '-e', dest = 'executable', action='store', \
    type='string', default = None, \
    help='Path to the blast executable', \
    metavar='FILE')
    p.add_option('--blast', '-b', dest = 'blast', action='store', \
    type='string', default = None, \
    help='Path to blast database', \
    metavar='FILE')
    p.add_option('--length', '-l', dest = 'length', action='store', \
    type='int', default = 100, help='The minimum length of reads to process')
    p.add_option('--nproc', '-n', dest = 'nproc', action='store', \
    type='int', default = 1, help='The number of processing cores to run')
    p.add_option('--single', '-s', dest = 'single', action='store', \
    type='string', default = None, help='If running a single locus, the locus to use')
    (options,arg) = p.parse_args()
    if not options.input:
        p.print_help()
        sys.exit(2)
    return options, arg


def get_known_sp(c, id):
    known = set()
    c.execute('SELECT distinct(spp) FROM reads where id = ?', (id,))
    for r in c.fetchall():
        known.add(r[0])
    return known

def record_printer(record, best_match):
    print 'query name=\t\t', record.query
    print 'query length=\t\t', record.query_length
    print 'subj name=\t\t', best_match
    print 'sbjct length=\t\t', record.alignments[0].length
    print 'align length=\t\t', record.alignments[0].hsps[0].align_length
    print record.alignments[0].hsps[0].query
    print record.alignments[0].hsps[0].match
    print record.alignments[0].hsps[0].sbjct
    print '\n'


def worker(input, output, options, db, exe):
    for f in iter(input.get, 'STOP'):
        # make per-process connection to the dbase
        con = sqlite3.connect(options.database)
        c = con.cursor()
        # get the core name
        core_name = f.split('-')[1]
        infile = os.path.join(os.getcwd(), f)
        out, err = subprocess.Popen("%s -query %s -db %s -evalue 1e-50 -dust 'yes' -max_target_seqs 5 -outfmt 5" % (exe, infile, db), shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate()
        blast = NCBIXML.parse(cStringIO.StringIO(out))
        loci = {'rbcl':set(), 'psba':set()}
        problems = {'rbcl':{}, 'psba':{}}
        group = set()
        scores = [] # this records the number of actual reads matching
        for record in blast:
            # TODO:  screen and remove short matches
            #pdb.set_trace()
            if record.alignments:
                #pdb.set_trace()
                query_name = ' '.join(record.query.split('_')[2:4])
                # capture locus name from fasta header
                #pdb.set_trace()
                locus = record.alignments[0].title.split('|')[-2]
                # skip short matches
                if record.alignments[0].hsps[0].align_length >= options.length:
                    best_match = ' '.join(record.alignments[0].title.split('|')[-1].strip(' ').split(' ')[0:2])
                    try:
                        if best_match not in loci[locus]:
                            loci[locus].add(best_match)
                    except:
                        pdb.set_trace()
                    if  query_name != best_match:
                        k = '%s:%s' % (query_name, best_match)
                        #record_printer(record, best_match)
                        if k in problems[locus].keys():
                            problems[locus][k] += 1
                        else:
                            problems[locus].update({k:1})
        # track problems on a per-core basis
        #metagroup[core_name] = problems
        # get the set of known species
        known = get_known_sp(c, core_name)
        # determine the intersection of the two loci - this should be sp.
        if options.single:
            inferred = loci[options.single]
        else:
            # present in both
            inferred = loci['psba'].intersection(loci['rbcl'])
        # show symmetric difference btw. known and inferred
        prnted = False
        results = [[],[]]
        ####
        # Species Level
        ####
        for i in known.intersection(inferred):
            results[0].append([core_name, i, ''])
        for i in known.difference(inferred):
            results[0].append([core_name, i, '-'])
        for i in inferred.difference(known):
            results[0].append([core_name, i, '+'])
        #####
        # Genus Level  
        ####
        known_genus = set([sp.split(' ')[0] for sp in known])
        inferred_genus = set([sp.split(' ')[0] for sp in inferred])
        for i in known_genus.intersection(inferred_genus):
            results[1].append([core_name, i, ''])
        for i in known_genus.difference(inferred_genus):
            results[1].append([core_name, i, '-'])
        for i in inferred_genus.difference(known_genus):
            results[1].append([core_name, i, '+'])
        output.put(results)

def main():
    path_filter = re.compile(".fsa$", re.IGNORECASE)
    # get and parse our command-line options
    options, args = interface()
    db = options.blast
    exe = options.executable
    file_count = 0
    metagroup = {}
    #outinfo = open('summary.out.txt','w')
    files = [f for f in os.listdir(options.input) if path_filter.search(f)]
    my = []
    n_procs = options.nproc
    task_queue = multiprocessing.Queue()
    done_queue = multiprocessing.Queue()
    # submit tasks
    for f in files:
        task_queue.put(f)
    # start some processes
    for i in range(n_procs):
        multiprocessing.Process(target=worker, args = (task_queue, done_queue, options, db, exe)).start()
    # get and print results
    #print 'Unordered resutls:'
    for i in range(len(files)):
        my.append(done_queue.get())
    #tell child processes to stop
    for i in range(n_procs):
        task_queue.put('STOP')
    outfile = open('species.out.txt','w')
    #pdb.set_trace()
    for core in my:
        for r in core[0]:
            outfile.write('%s\t%s\t%s\n' % (r[0], r[1], r[2]))
    outfile.close()
    outfile = open('genus.out.txt','w')
    for core in my:
        for r in core[1]:
            outfile.write('%s\t%s\t%s\n' % (r[0], r[1], r[2]))
    outfile.close()
    

if __name__ == '__main__':
    main()

