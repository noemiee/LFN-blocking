#!/usr/bin/env python
# coding: utf-8


import os
import glob
import numpy as np
import struct
import pickle
import time
from scipy.sparse import coo_matrix
os.environ["NUMBA_DISABLE_INTEL_SVML"] = "1"
from numba import jit, set_num_threads, prange



def readpart10(directory,i,numpart,nspec=1):
# read one flexpart output file
# input : directory (file location), i is a string or an index corresponding to a date in the 'dates' file,
#         numpart total number of particles to read in that file, nspec : number of species


    # determine path for the file to read
    if isinstance(i,(int, float)):
        dates_file = open(directory + "/dates", "r")
        dates = dates_file.read()
        dates_file.close()
        dates = dates.splitlines()
        file = "partposit_" + str(dates[i])
    elif isinstance(i,(str)):
        file = "partposit_" + i
    file = directory + "/" + file

    # if numpart not given, calculate numpart (total) and return it
    if isinstance(numpart, (str)) :
        # the following assumes that the file has the exact format written in partouput.f90 of Flexppart v10.4
        # then the size of the file is exactly 4 * [3+(numpart+1)*(14+nspec)] bytes
        nbytes = os.path.getsize(file)
        numpart = round((nbytes/4+3)/(14+nspec)-1);
        output = numpart
        print(int(numpart))
        return

    # initialize arrays to hold return values
    outheader = 0
    npoint = np.zeros((numpart,1), dtype = np.int32)
    xyz = np.zeros((numpart,3), dtype = np.float32)
    #itramem = np.zeros((numpart,1), dtype = np.int32)
    #vars = np.zeros((numpart,7), dtype = np.float32)
    #xmass = np.zeros((numpart,nspec), dtype = np.float32)

    #  read file header (time stamp)
    try :
        particle_file = open(file, "rb")
    except OSError:
        print("Could not open/read file:", file)
        sys.exit()

    dummy = particle_file.read(4) # read first dummy value
    packed_outheader = particle_file.read(4)
    outheader = struct.unpack('@i',packed_outheader)
    dummy = particle_file.read(4) # read next dummy value

    particle_file.close()

    # check that the number of particles is at least 1
    if numpart<1 :
        print('ERROR : number of particles must be at least 1. \n')
        return

    # note : FLEXPART writes the output as
    #         write(unitpartout) npoint(i),xlon,ylat,ztra1(i),
    #    +    itramem(i),topo,pvi,qvi,rhoi,hmixi,tri,tti,
    #    +    (xmass1(i,j),j=1,nspec)
    # written as 4-byte numbers with empty 4 byte value before and afterwards
    # so the number of 4-byte values to read for each particle/ trajectory is 14+nspec
    nvals = 14 + nspec;

    # read all data as 4-byte int values
    particle_file = open(file, "rb")
    dummy = particle_file.read(3*4) # skip header
    data = particle_file.read(nvals*numpart*4)
    dataAsInt = struct.unpack('@'+nvals*numpart*'i', data)
    dataAsInt = np.reshape(dataAsInt, (numpart,nvals))
    particle_file.close()
    
    
    # select values for output
    npoint = dataAsInt[0:numpart,1]
    xyz = dataAsFloat[0:numpart,2:5]
    #itramem = dataAsInt[0:numpart,5]
    #vars = dataAsFloat[0:numpart,6:13]
    #xmass = dataAsFloat[0:numpart,13:(13+nspec)]

    return [outheader, npoint, xyz]#, itramem, vars, xmass]



@jit("i8[:,:](i4, f8[:], f8[:], i4, f8[:,:])", nopython=True, parallel=True)
def evaluate_cells(numpart, x, y, ppc, grid):
    # given the final position of particles (ordered according to their release positions) returns an array of (initial, final) cell indices.
    indices = np.empty((numpart, 2), dtype = np.int64)
    for i in prange(numpart):
        index = np.where((x[i]>=grid[:,0]) & (x[i]<=grid[:,6]) & (y[i]<=grid[:,1]) & (y[i]>= grid[:,7]))[0]
        if len(index) > 0: # if the particle ends up in the domain
            indices[i,0] = int(i/ppc)
            indices[i,1] = index[0]
        else :
            # particle out of gridded domain
            #print('the final grid does not cover the whole domain', x,y)
            indices[i,0] = int(i/ppc)
            indices[i,1] = -1
    return indices


@jit("void(f8[:,:], i4, i8[:,:])", nopython=True)
def makeP(P, numpart, indices):
    for n in range(numpart):
        if indices[n,1] == -1: # ignore particles that went out of the gridded domain
            pass
        else:
            P[indices[n,0], indices[n,1]] += 1.0



@jit("void(f8[:,:], i4)", nopython=True, parallel=True)
def normalize_out(A, Ncells) :
    out_strenghts = np.sum(A, axis = 1)
    for i in prange(Ncells) :
        for j in prange(Ncells) :
                if (A[i,j] != 0) :
                    A[i,j] /= out_strenghts[i]

@jit("void(f8[:,:], i4)", nopython=True, parallel=True)
def normalize_in(A, Ncells) :
    in_strenghts = np.sum(A, axis = 0)
    for i in prange(Ncells) :
        for j in prange(Ncells) :
                if (A[i,j] != 0) :
                    A[i,j] /= in_strenghts[j]


def main():
    
    ##################################################################################################################
    Fp_output_folder = './output/'

    #################################################################################################################
    # load grid
    grid_file_location = '../grid/grid.pckl'
    f = open(grid_file_location, 'rb')
    grid = pickle.load(f)
    f.close()
    ##################################################################################################################

    Ncells = len(grid)
    ppc = 1500 # num of particles released per cell (particles per cell)

    # compute number of particles per zone
    numpart = Ncells*ppc
    print('number of simulated particles : ', numpart)

    ##################################################################################################################
    ###### read x, y positions from flexpart output files ############################################################

    x = np.zeros(numpart)
    y = np.zeros(numpart)

    # file to read
    file =  'end'
    print('reading file : ', file)
    # extract x and y positions
    out = readpart10(Fp_output_folder , file, numpart)
    cells = out[1][:] # corresponds to release number
    x_tmp = out[2][:,0]
    y_tmp = out[2][:,1]
    # make sure that the order is cell by cell (i.e. the ppc first particles <-> first cell, the ppc next particles <-> second cell, etc.)
    I = np.argsort(cells)
    x = x_tmp[I]
    y = y_tmp[I]

    ##################################################################################################################
    ###### compute flow network ######################################################################################

    # evaluate initial and final cell for each trajectory
    start_t = time.time()
    indices = evaluate_cells(numpart, x, y, ppc, grid)  # list of (initial, final) positions
    end_t = time.time()
    print('Time to evaluate indices: ', end_t-start_t)


    # construct matrix of transition probabilities P and adjacency matrix Adj
    start_t = time.time()
    P = np.zeros((Ncells, Ncells), dtype = np.float64) # transport matrix
    makeP(P, numpart, indices)
    F = np.copy(P)
    
    # normalize P
    normalize_out(P, Ncells)
    # normalize F
    normalize_in(F, Ncells)
    # convert matrices to sparse format
    P = coo_matrix(P)
    F = coo_matrix(F)

    end_t = time.time()
    print('Time to evaluate transition and adjacency matrices: ', end_t-start_t)



    ##################################################################################################################
    # strore matrix

    P_file_location = "tmp/P.pckl"
    F_file_location = "tmp/F.pckl"

    f = open(P_file_location, 'wb')
    pickle.dump(P, f)
    f.close()
    f = open(F_file_location, 'wb')
    pickle.dump(F, f)
    f.close()

    ##################################################################################################################
    # delete the content of output folder to save space
    #if os.path.isdir('./output/') and os.listdir('./output/'):
    #    files = glob.glob('./output/*')
    #    for f in files : 
    #        os.remove(f)

if __name__ == "__main__":
    main()
