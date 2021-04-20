#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
#import struct
import pickle
from datetime import datetime
import glob


def main() :

    # we get all available dates by reading the time stamps of all adjacency matrices)
    all_dates = [file[-19:-5] for file in glob.glob('./data/networks/P_*')]
    # convert those to datetime objects
    all_dates = [datetime.strptime(d, '%Y%m%d%H%M%S') for d in all_dates]
    # sort the dates
    all_dates.sort()

    ndates = len(all_dates)

    print("Number of time steps: ", ndates)
    
    # before creating the plots, check the existence of the directory to store files
    if not(os.path.exists('/data/degree-fields')):
        try:
            os.makedirs('./data/degree-fields')
        except OSError:
            print ("Creation of the directory failed.")
        else:
            print ("Creation of the directory successfull.")
 


    # loop over time steps
    idx = -1
    for d in all_dates:
        print("Computing degree - ", d)
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
        num_part_released = 1500
        for i in range(P.shape[0]):
            if Kout[i] > ALPHA*num_part_released:
                nb += 1.
        if nb/P.shape[0]*100>1.0:
            print("the number of particles released in this step (",d,") is not sufficient!")


        ##############################################################
        # store degree values
        file_location = "./data/degree-fields/Kout"+("%03d" % idx)+".pckl"
        f = open(file_location, 'wb')
        pickle.dump(Kout, f)
        f.close()

        


if __name__ == "__main__":
    main()
