#!/usr/bin/env python
# coding: utf-8

#import os
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

    # loop over time steps
    idx = -1
    for d in all_dates:
        print('computing entropy - ', d)
        idx+=1

        # load weighted adjacency matrix ("out-weights")
        file ='./data/networks/P_'+d.strftime('%Y%m%d%H%M%S')+'.pckl'
        f = open(file, 'rb')
        P = pickle.load(f)
        f.close()

        P = P.toarray()
        if idx == 0 :
            print("Network size: ", P.shape[0])
        
        ################################################################
        # compute entropy
        H1out = [-1.*np.sum(P[i,P[i,:]!=0]*np.log(P[i,P[i,:]!=0])) for i in range(P.shape[0])]
        H1out = np.array(H1out)

        ##############################################################
        # store entropy values
        file_location = "./data/entropy-fields/H1out"+("%03d" % idx)+".pckl"
        f = open(file_location, 'wb')
        pickle.dump(H1out, f)
        f.close()

        


if __name__ == "__main__":
    main()
