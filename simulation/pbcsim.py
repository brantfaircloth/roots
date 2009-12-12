#!/usr/bin/env python
# encoding: utf-8
"""
pbcsim.py

Created by Brant Faircloth on 2009-11-30.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import pdb
import re
import sys
import numpy 
import random
import optparse
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
# make sure we use sqlite supporting foreign keys
from pysqlite2 import dbapi2 as sqlite3


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"
    p = optparse.OptionParser(usage)
    p.add_option('--input', '-i', dest = 'input', action='store', \
    type='string', default = None, \
    help='Path to the sequences from which we will randomly sample', \
    metavar='FILE')
    p.add_option('--output', '-o', dest = 'output', action='store', \
    type='string', default = None, \
    help='Path to the file to hold our results', \
    metavar='FILE')
    p.add_option('--database', '-d', dest = 'database', action='store', \
    type='string', default = 'roots.sqlite', \
    help='Path to the database to hold our results', \
    metavar='FILE')
    p.add_option('--sample', '-s', dest = 'sample', action='store', \
    type='float', default = 15, help='The likely number of species per sample')
    p.add_option('--variance', '-v', dest = 'sample_sd', action='store', \
    type='float', default = 5.0, \
    help='The variance (standard deviation) in number of species/sample')
    p.add_option('--sampling_freq', '-f', dest = 'sampling_freq', action='store', \
    type='float', default = 0.95, help='The frequency with which we are \
    likely to "draw" DNA from the root pool')
    p.add_option('--cores', '-c', dest = 'cores', action='store', \
    type='int', default = 96,
    help='The (total) number of cores samples')
    p.add_option('--reads', '-r', dest = 'reads', action='store', \
    type='int', default = 1000,
    help='The (total) number of reads per core')
    (options,arg) = p.parse_args()
    if not options.input:
        p.print_help()
        sys.exit(2)
    return options, arg

def length(handle):
    '''determine length of a fasta file file'''
    count = 0
    for line in open(handle,'r'):
        if line.startswith('>'):
            count += 1
    return count

def species_per_core(mean, sd, count):
    '''compute the number of species per virtual soil core'''
    #pdb.set_trace()
    return abs(numpy.random.normal(mean, sd, count).round())

def species_sampler(handle, sample):
    '''randomly select some sequences from a fasta file'''
    #pdb.set_trace()
    lines = length(handle)
    rsample = random.sample(xrange(lines), int(sample))
    count = 0
    random_sequences = ()
    for record in SeqIO.parse(open(handle, 'r'), 'fasta'):
        if count in rsample:
            random_sequences += (record,)
        count += 1
    return random_sequences

def root_freq(sample, roots):
    '''compute the true frequency of sequences in solution'''
    return (numpy.random.mtrand.dirichlet([1] * sample) * roots).round()

class Error():
    '''A generic class for generating sequence errors in PCR products
    and sequencing reads'''
    def __init__(self, **kwargs):
        self.kwargs = kwargs
        if 'h_rate' in self.kwargs.keys():
            #self.approx_homo_error_rate = sum(numpy.random.poisson(self.kwargs['h_rate'], 10000) > 0)/10000.
            self._homogen()
        #if rate in self.kwargs.keys():
        self.homo_error_relative = numpy.array([])
        self.homo_error_overall = numpy.array([])
        self.other_error_overall = numpy.array([])
    
    def _homogen(self):
        homo = ''
        for b in ['Aa','Cc', 'Gg', 'Tt']:
            homo += (('[%s]{%s,}|') % (b, self.kwargs['h_length'] + 1))
        homo = homo[:-1]
        #pdb.set_trace()
        self.re_homo = re.compile(homo)
    
    def homopolymer(self, seq):
        matches = tuple(re.finditer(self.re_homo, str(seq)))
        seq_list = list(str(seq)) 
        if matches:
            #for m in matches:
            #    print m.group()
            # is it insertion or deletion (multinomial) with p = 0.5
            insertions = numpy.random.binomial(1,0.5, len(matches))
            # let's treat the homopolymer error creation as a poisson process
            # so that we can decide both whether there is error at a particular
            # homopolymer run and **how** much error there is (i,e. an indel
            # of 0,1,2,3...,n bases).  We'll use the binomial above to determine
            # if we have an insertion or deletion
            errors = numpy.random.poisson(self.kwargs['h_rate'], len(matches))
            # go in reverse
            for pos, m in enumerate(matches[::-1]):
                if insertions[pos] and errors[pos]:
                    #pdb.set_trace()
                    seq_list = seq_list[:m.start()] + list(m.group() + m.group()[0] * errors[pos]) + seq_list[m.end():]
                elif not insertions[pos] and errors[pos]:
                    if m.end() - m.start() > errors[pos]:
                        seq_list = seq_list[:m.start()] + list(m.group()[:len(m.group()) - errors[pos]]) + seq_list[m.end():]
                    else:
                        pass
                else:
                    pass
            #pdb.set_trace()
            self.homo_error_overall = numpy.append(self.homo_error_overall, float(sum(errors))/len(seq_list))
            self.homo_error_relative = numpy.append(self.homo_error_relative, float(sum(errors))/len(matches)) 
        else:
            self.homo_error_overall = numpy.append(self.homo_error_overall, 0.)
            self.homo_error_relative = numpy.append(self.homo_error_relative, 0.)
        return Seq(''.join(seq_list), SingleLetterAlphabet())
        #print '\n'
        #pdb.set_trace()

    def other(self, seq):
        seq_list = list(str(seq))
        # lets treat "regular errors" as a binomial process
        errors = numpy.random.binomial(1, self.kwargs['rate'], len(seq_list))
        if 1 not in errors:
            self.other_error_overall = numpy.append(self.other_error_overall, 0.)
        else:
            bases = {0:'A',1:'C',2:'G',3:'T'}
            for pos, e in enumerate(seq_list[::-1]):
                if errors[pos] == 1:
                    #pdb.set_trace()
                    # see if it's an insertion or deletion
                    if numpy.random.binomial(1,0.5) and errors[pos]:
                        # just pick a random base and insert it
                        new_base = bases[numpy.random.random_integers(0,3)]
                        seq_list = seq_list[:pos] + [new_base] + seq_list[pos:]
                    else:
                        seq_list = seq_list[:pos] + seq_list[pos + 1:]
            #pdb.set_trace()
            self.other_error_overall = numpy.append(self.other_error_overall, float(sum(errors))/len(seq_list))
        return Seq(''.join(seq_list), SingleLetterAlphabet())


def read_lengths(seq, mean, sd):
    '''sample from the entire length of the read, in both
    directions, returning the forward (potentially partial) read and the
    reverse (potentially partial) read (the revcomp)'''
    l,r = numpy.random.normal(mean, sd, 2).round()
    return (seq[:int(l)], seq[len(seq)-int(r):].reverse_complement())


def dna_sample(root_freq, sampling_freq):
    '''since adding DNA to PCR reactions is basically a sampling process,
    recreate that process by sampling the available pool of species - i,e.
    we are probably going to drop some low-count species here, a process at
    which we're interesting in looking'''
    root_pop = numpy.array([], dtype = 'int64')
    # using our relative freqs from the dirichlet, expand those counts to give
    # the root population we have in our virtual sample
    for pos,count in enumerate(root_freq):
        root_pop = numpy.append(root_pop, numpy.array([pos] * count))
    # NOTE:  we could probably just numpy.random.shuffle() and then slice
    # we're assuming we draw sampling_freq of the population when we pipet
    # so, we'll randomly sample, without replacement, the index values of the
    # elements created in the root_pop array that is above
    root_pop_sample = numpy.array(random.sample(xrange(len(root_pop)), int(sampling_freq*len(root_pop))))
    #pdb.set_trace()
    root_sample = root_pop[root_pop_sample]
    # we'll reform the root_freq post sampling
    post_sample_counts = numpy.zeros(len(root_freq))
    for r in root_sample:
        post_sample_counts[r] += 1
    #pdb.set_trace()
    # now, we'll reindex the root_pop array and return our sampled
    # root population
    return root_sample, post_sample_counts

def tables(c):
    '''create sqlite tables for root data'''
    c.execute('''PRAGMA foreign_keys=ON;''')
    try:
        c.execute('''DROP TABLE reads''')
        c.execute('''DROP TABLE cores''') 
    except:
        pass
    c.execute('''create table cores (
    id int PRIMARY KEY,
    mid_tag text,
    true_count int,
    sample_count int,
    homo_error real,
    other_error real,
    all_error real
    )''')
    c.execute('''create table reads (
    id int,
    indiv_id int,
    side text,
    spp text,
    read_length real,
    homo_error real,
    other_error real,
    all_error real,
    FOREIGN KEY(id) REFERENCES cores(id)
    )''')

def insert_read_row(c, core_index, individual_index, spp, reads, h_error, s_error):
    '''docstring for insert_read_row'''
    for i, side in enumerate(['l','r']):
        c.execute('''INSERT INTO reads (id, indiv_id, side, spp, read_length, 
        homo_error, other_error, all_error) VALUES (?,?,?,?,?,?,?,?)''', 
        (core_index, individual_index, side, spp, len(reads[i]), 
        h_error[i].round(3), s_error[i].round(3), 
        (h_error[i]+s_error[i]).round(3)))

def write_reads(i, side, fsa, core_index, individual_index, spp, reads, h_error, s_error):
    '''generator to hold the sequence reads for efficient writing'''
    header = '%s_%s_%s' % (core_index, spp.replace(' ','_'), side)
    yield SeqRecord(reads[i], id = header, name = header, description = header)

def main():
    # get and parse our command-line options
    options, args = interface()
    # create a dbase
    con = sqlite3.connect(options.database)
    c = con.cursor()
    # create some tables
    tables(c)
    # commit the additions
    con.commit()
    # create a file for the generated sequence reads
    fsa = open(options.output, 'w')
    # generate some counts of species in each virtual root core
    cores = species_per_core(options.sample, options.sample_sd, options.cores)
    for core_index, core in enumerate(cores):
        #pdb.set_trace()
        # randomly select some sequences from the group that will be in our virtual
        # soil core
        core_species = species_sampler(options.input, core)
        # use a dirichlet to generate random relative frequencies for the roots
        # in the virtual soil core
        #pdb.set_trace()
        core_true_freq = root_freq(core, options.reads)
        # since adding DNA to PCR reactions is basically a sampling process,
        # recreate that process by sampling the available pool of species - i,e.
        # we are probably going to drop some low-count species here, at which we're
        # interesting in looking
        #
        # Similarly, our work in the lab, separating roots from soil will also
        # mimic a sampling process (perhaps not as random)
        core_sample, core_sample_freq = dna_sample(core_true_freq, options.sampling_freq)
        # we know that the PCR and sequencing processes entail incorporation of
        # some error to each read.  Roche 454 would tell us that it's 1%
        # cumulative from their E. coli work (per Roche rep.) so here, we're going
        # to
        #
        # create PCR and Sequencing error instances
        pcr_error = Error(rate = 2.6e-5)
        sequencing_error = Error(h_length = 3, h_rate = 0.15, rate = 0.01)
        # error rates for each error type are held in:
        # self.homo_error_relative (this is the rate at each homo run)
        # self.homo_error_overall (this is the homo error rate, over all bases)
        # self.other_error_overall (this is the non-homo error rate, over all bp)
        #pdb.set_trace()
        c.execute('''INSERT INTO cores (id, true_count, sample_count) VALUES 
        (?,?,?)''', (core_index, int(sum(core_true_freq > 0)), 
        int(sum(core_sample_freq > 0))))
        iterator = itertools.chain()
        for individual_index, individual in enumerate(core_sample):
            # get the the voucher read for the species
            record  = core_species[int(individual)]
            # add some error to the PCR sequences.  This is likely to be a small
            # rate and this is a vastly oversimplified approximation of a "real"
            # error generation process (which would be exponenential)
            # TODO: make process more real
            pcr_seq = pcr_error.other(record.seq)
            #pdb.set_trace()
            # get read lengths for a particular fragment from both ends
            reads = read_lengths(pcr_seq, 407.5, 126.1)
            # add some error to those reads
            reads = sequencing_error.homopolymer(reads[0]), sequencing_error.homopolymer(reads[1])
            h_error = sequencing_error.homo_error_overall[-2:]
            reads = sequencing_error.other(reads[0]), sequencing_error.other(reads[1])
            s_error = sequencing_error.other_error_overall[-2:]
            # get the spp name
            spp = ' '.join(record.description.split('|')[4].split(' ')[1:3])
            insert_read_row(c, core_index, individual_index, spp, reads, h_error, s_error)
            for i, side in enumerate(['l','r']):
                iterator = itertools.chain(iterator, write_reads(i, side, fsa, core_index, individual_index, spp, reads, h_error, s_error))
        h_err_over = numpy.mean(sequencing_error.homo_error_overall).round(3)
        s_err_over = numpy.mean(sequencing_error.other_error_overall).round(3)
        a_err_over = numpy.mean(sequencing_error.homo_error_overall + sequencing_error.other_error_overall).round(3)
        c.execute('''UPDATE cores set homo_error = ?, other_error = ?, 
        all_error = ? WHERE id = ?''', (h_err_over, s_err_over, a_err_over, core_index))
        con.commit()
        SeqIO.write(iterator, fsa, "fasta")
        pdb.set_trace()
    cur.close()
    con.close()
        
    pdb.set_trace()

if __name__ == '__main__':
    main()
