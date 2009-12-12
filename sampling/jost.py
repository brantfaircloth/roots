#!/usr/bin/env python
# encoding: utf-8
"""
jostD.py

Created by Brant Faircloth on 2009-12-07.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""
   
import pdb
import sqlite3
import numpy

def groupList(cur):
    cur.execute('SELECT distinct(Quadrat) from plot where Quadrat IS NOT NULL')
    return cur.fetchall()

def speciesByQuadrat(cur, q):
    cur.execute('''SELECT Latin, count(Latin) from plot where Quadrat = "%s" and Status="alive" and Stem = "main" group by Latin''' % q)
    return cur.fetchall()
    
def main():
    out = open('QuadratJostDData.txt', 'w')
    con = sqlite3.connect('/Users/bcf/Documents/UCLA/Lab/Organisms/Roots/GIS/BCI2005.sqlite')
    cur = con.cursor()     
    # get list of quadrats
    qdrats = groupList(cur)
    out.write("Quadrat\tShannonsH\tJostsD\n")
    for q in qdrats:
        # select count of all records for the quadrat
        #q_count = allQuadratRecords(cur, q)
        # select the records for the quadrat grouped by species
        sp_count = speciesByQuadrat(cur, q)
        # turn sp_count into an array
        sp_array = []
        for s in sp_count:
            sp_array.append(float(s[1]))
        sp_array = numpy.array(sp_array)
        #pdb.set_trace()
        shannons_h = - sum((sp_array/sum(sp_array) * numpy.log(sp_array/sum(sp_array))))
        jost_d = numpy.exp(shannons_h)
        out.write(('\'%s\'\t%s\t%s\n') % (q[0], shannons_h, jost_d))
    out.close()


if __name__ == '__main__':
    main()

