#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
#import struct
import pickle
from datetime import datetime
from graph_tool.all import *
import glob
os.environ["NUMBA_DISABLE_INTEL_SVML"] = "1"
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import time
from numba import jit
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                               AutoMinorLocator)

@jit(nopython=True)
def make_edges_out(P):
    # input: weighted adjacency matrix
    # output: list of edges associated with a "distance" defined based on the flux exchanged between nodes
    edges = []
    for row in range(P.shape[0]):
        for col in range(P.shape[1]):
            if P[row, col]!=0 :
                d = 1/P[row,col]
                edges.append((row, col, d ))
    return edges



def main() :
    openmp_set_num_threads(4)

    # we get all available dates by reading the time stamps of all adjacency matrices)
    all_dates = [file[-19:-5] for file in glob.glob('./data/networks/P_*')]
    # convert those to datetime objects
    all_dates = [datetime.strptime(d, '%Y%m%d%H%M%S') for d in all_dates]
    # sort the dates
    all_dates.sort()

    ndates = len(all_dates)

    print("Number of time steps: ", ndates)
    

    # before creating the plots, check the existence of the directory to store files
    if not(os.path.exists('/data/closeness-fields')):
        try:
            os.makedirs('./data/closeness-fields')
        except OSError:
            print ("Creation of the directory failed.")
        else:
            print ("Creation of the directory successfull.")

    # loop over time steps
    idx = -1
    for d in all_dates:
        idx+=1

        # load weighted adjacency matrix ("out-weights")
        file ='./data/networks/P_'+d.strftime('%Y%m%d%H%M%S')+'.pckl'
        f = open(file, 'rb')
        P = pickle.load(f)
        f.close()

        P = P.toarray()
        Adj = np.where(P>0, 1, 0)

        if idx == 0 :
            print("Network size: ", P.shape[0])
        
        ################################################################
        # test that the number of particles released is sufficiently high
        Kout = np.sum(Adj, axis = 1)
        ALPHA = 0.1
        nb = 0.0
        num_part_released = 900
        for i in range(P.shape[0]):
            if Kout[i] > ALPHA*num_part_released:
                nb += 1.
        if nb/P.shape[0]*100>1.0:
            print("the number of particles released in this step (",d,") is not sufficient!")

        ################################################################
        # construct directed network using graph_tool
        net_A_out = Graph(directed=True)
        vlist_out = net_A_out.add_vertex(P.shape[0])
        # set links according to the adjacency matrix, the "distance" between nodes is defined in the fucntion make_edges_out and depends on the flux between the nodes (!not related to geodesic distance!)
        distances_out = net_A_out.new_edge_property("double")
        eprops = [distances_out]
        edges = make_edges_out(P)
        net_A_out.add_edge_list(edges, eprops=eprops)


        ##############################################################
        # closeness
        print(d, ' -- computing closeness ', end='')
        t0 = time.time()

        harmonic_closeness_out = closeness(net_A_out, weight = distances_out, norm = True, harmonic = True)
        harmonic_closeness_out = harmonic_closeness_out.get_array()
        t1 = time.time()

        print(' -- total time : ', t1 -t0)
        ##############################################################
        # store closeness values
        file_location = "./data/closeness-fields/hcl_"+("%03d" % idx)+"_out.pckl"
        f = open(file_location, 'wb')
        pickle.dump(harmonic_closeness_out, f)
        f.close()

        


if __name__ == "__main__":
    main()
