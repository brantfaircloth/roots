#!/usr/bin/env python
# encoding: utf-8
"""
sp_count.py

Created by Brant Faircloth on 2009-12-20.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import pdb
import numpy
# make sure we use sqlite supporting foreign keys
from pysqlite2 import dbapi2 as sqlite3


def rect_dist(input, dist=2.5):
    return numpy.vstack((input-dist, input+dist))

def main():
    con = sqlite3.connect('/Users/bcf/Documents/UCLA/Lab/Organisms/roots/GIS/BCI2005.sqlite')
    cur = con.cursor()
    # choose a random (x) point in space from [0,999)
    x = numpy.random.uniform(low=0, high=1000, size=500)
    # choose a random (y) point in space from [0,999)
    y = numpy.random.uniform(low=0, high=500, size=500)
    # get x points ± 5m from x
    xs = rect_dist(x)
    # get y points ± 5m from y
    ys = rect_dist(y)
    results = []
    for i in xrange(500):
        x_pos = xs[:, i]
        y_pos = ys[:, i]
        cur.execute('''SELECT count(distinct(Latin)) FROM plot where gx between ? and ? and gy between ? and ? and Status = "alive"''', (x_pos[0], x_pos[1], y_pos[0], y_pos[1]))
        results.append(cur.fetchall()[0][0])
        #print x_pos, y_pos, results[i]
    pdb.set_trace()
    results = numpy.array(results)
    print 'mean= ', results.mean()
    print 'std= ', results.std()


if __name__ == '__main__':
    main()

