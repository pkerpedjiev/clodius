#!/usr/bin/python

from __future__ import print_function

import clodius.fpark as cf
import h5py
import math
import negspy.coordinates as nc
import numpy as np
import os
import os.path as op
import pybedtools as pbt
import random
import sys
import argparse

# all entries are broked up into ((tile_pos), [entry]) tuples
# we just need to reduce the tiles so that no tile contains more than
# max_entries_per_tile entries
# (notice that [entry] is an array), this format will be important when
# reducing to the most important values
def reduce_values_by_importance(entry1, entry2, max_entries_per_tile=100, reverse_importance=False):
    if reverse_importance:
        combined_entries = sorted(entry1 + entry2,
                key=lambda x: -float(x[-1]))
    else:
        combined_entries = sorted(entry1 + entry2,
                key=lambda x: float(x[-1]))
    return combined_entries[:max_entries_per_tile]

def main():
    parser = argparse.ArgumentParser(description="""
    
    python tileBedFileByImportance.py file.bed
""")

    parser.add_argument('bedfile')
    parser.add_argument('--importance-column', type=str,
            help='The column containing information about how important'
            "that row is. If it's absent, then use the length of the region."
            "If the value is equal to `random`, then a random value will be"
            "used for the importance (effectively leading to random sampling)")
    parser.add_argument('--assembly', type=str, default='hg19',
            help='The genome assembly to use')
    parser.add_argument('--max-per-tile', type=int, default=100)
    parser.add_argument('--tile-size', default=1024)
    parser.add_argument('-o', '--output-file', default='/tmp/tmp.hdf5')

    #parser.add_argument('argument', nargs=1)
    #parser.add_argument('-o', '--options', default='yo',
    #					 help="Some option", type='str')
    #parser.add_argument('-u', '--useless', action='store_true', 
    #					 help='Another useless option')

    args = parser.parse_args()

    if op.exists(args.output_file):
        os.remove(args.output_file)

    # The tileset will be stored as an hdf5 file
    print("Writing to: {}".format(args.output_file), file=sys.stderr)
    f = h5py.File(args.output_file, 'w')

    # We neeed chromosome information as well as the assembly size to properly
    # tile this data
    tile_size = args.tile_size
    chrom_info = nc.get_chrominfo(args.assembly)
    assembly_size = chrom_info.total_length
    max_zoom = int(math.ceil(math.log(assembly_size / tile_size) / math.log(2)))


    # store some meta data
    d = f.create_dataset('meta', (1,), dtype='f')

    d.attrs['zoom-step'] = 1            # we'll store data for every zoom level
    d.attrs['max-length'] = assembly_size
    d.attrs['assembly'] = args.assembly
    d.attrs['chrom-names'] = nc.get_chromorder(args.assembly)
    d.attrs['chrom-sizes'] = nc.get_chromsizes(args.assembly)
    d.attrs['chrom-order'] = nc.get_chromorder(args.assembly)
    d.attrs['tile-size'] = tile_size
    d.attrs['max-zoom'] = max_zoom =  math.ceil(math.log(d.attrs['max-length'] / tile_size) / math.log(2))
    d.attrs['max-width'] = tile_size * 2 ** max_zoom

    bed_file = pbt.BedTool(args.bedfile)


    print("max_zoom:", max_zoom)

    all_parts = []

    def line_to_np_array(line):
        '''
        Convert a bed file line to a numpy array which can later
        be used as an entry in an h5py file.
        '''

        if args.importance_column is None:
            importance = line.stop - line.start
        elif args.importance_columns == 'random':
            imporance = random.random()
        else:
            importance = line.fields[importance]

        genome_start = nc.chr_pos_to_genome_pos(str(line.chrom), line.start, args.assembly)
        genome_end = nc.chr_pos_to_genome_pos(line.chrom, line.stop, args.assembly)
        parts = [genome_start, genome_end] + map(str,line.fields[3:]) + [importance]

        return parts


    dset = sorted([line_to_np_array(line) for line in bed_file], key=lambda x: x[0])
    #dset = [line_to_np_array(line) for line in bed_file]

    '''
    for i in range(len(dset[0])):
        col = [d[i] for d in dset]
        f.create_dataset('{}_{}'.format(max_zoom, i), data=col, compression='gzip')
    '''
    f.create_dataset('{}'.format(max_zoom), data=dset, compression='gzip')

    curr_zoom = max_zoom - 1

    tile_width = tile_size

    # each entry in pdset will correspond to the values visible for that tile
    pdset = cf.ParallelData(dset).map(lambda x: (int(x[0] / tile_width), [x]))

    print('pd:', pdset.take(10))


    while curr_zoom >= 0:
        print("x1:", pdset.take(2))
        pdset = pdset.reduceByKey(lambda e1,e2: reduce_values_by_importance(e1, e2, 
            max_entries_per_tile = args.max_per_tile))
        #tile_nums_values = [(int(d[0] / tile_size), d) for d in dset]
        pdset = pdset.map(lambda x: (x[0] / 2, x[1]))

        new_dset = [item for sublist in [d[1] for d in pdset.collect()] for item in sublist]
        print("len(new_dset)", len(new_dset))
        f.create_dataset('{}'.format(curr_zoom), data=new_dset, compression='gzip')

        curr_zoom -= 1
        


    #print("dset:", dset)
    

if __name__ == '__main__':
    main()

