#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
#import struct
import pickle
from datetime import datetime
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import time
from netCDF4 import Dataset, num2date
import matplotlib.tri as tri


def viz(lon0, lon1, lat0, lat1, Xi, Yi, e_interp, minimum, maximum, lon_r, lat_r, z_data, title, file_name) :

    # create figure 
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (16,8))
    plt.title(title, fontsize=16)

    # create basemap on which the field will be plotted 
    m = Basemap(projection='merc',llcrnrlat=lat0,urcrnrlat=lat1,llcrnrlon=lon0,urcrnrlon=lon1,lat_ts=20,resolution='l', ax = ax)
    m.drawcoastlines(zorder=2, color="black",linewidth=2.2)
    m.drawcountries(zorder=2,color="black", linewidth=2)
    m.drawparallels(np.arange(-90.,90.,10.),labels=[1,0,0,0], linewidth = 0, fontsize=14)
    m.drawmeridians(np.arange(-180.,180.,20.), labels=[0,0,0,1], linewidth = 0, fontsize=14, rotation=45)

    # filled contour plot of the interpolated entropy
    xxi, yyi = m(Xi, Yi)
    levels = np.linspace(minimum, maximum, 80)
    cnt = m.contourf(xxi, yyi, e_interp, levels = levels, cmap="viridis", extend = 'both')#, vmin = minimum, vmax = maximum)
    # set edges of contours the same color as faces
    for c in cnt.collections:
        c.set_edgecolor("face")
    # add colobar
    cbar = m.colorbar(cnt, ax= ax, extend = 'both')

    # plot geopotential height contours
    lons, lats = np.meshgrid(lon_r,lat_r)
    x,y = m(lons,lats)
    cs = m.contour(x, y, z_data/10., linewidths=1.5, colors = 'w') #:1333
    plt.clabel(cs, inline=1, fontsize=8, fmt='%.0f')

    plt.savefig(file_name, format='pdf', dpi=300,  bbox_inches="tight")
    plt.close()






def main() :

    # we get all available dates by reading the time stamps of all adjacency matrices)
    all_dates = [file[-19:-5] for file in glob.glob('./data/networks/P_*')]
    # convert to datetime objects
    all_dates = [datetime.strptime(d, '%Y%m%d%H%M%S') for d in all_dates]
    # sort the dates
    all_dates.sort()

    ndates = len(all_dates)
    
    # load the network grid in order to map each network node to a physical location
    file ='./data/grid/grid.pckl'
    f = open(file, 'rb')
    Grid = pickle.load(f)
    f.close()
    # compute the coordinates of the center of each network cell
    cell_centers_x = []
    cell_centers_y = []
    for cell in range(len(Grid)):
        xmid = (Grid[cell][0]+Grid[cell][2])/2.
        ymid = (Grid[cell][1]+Grid[cell][-1])/2.
        cell_centers_x.append(xmid)
        cell_centers_y.append(ymid)

    cell_centers_x = np.array(cell_centers_x)
    cell_centers_y = np.array(cell_centers_y)

    ###### load geopotential data at 500hPa for plotting

    nc = Dataset('./data/gph/gph2010-500.nc')

    lat = nc.variables['latitude'][:]
    lon = nc.variables['longitude'][:]
    lat_r = lat[(lat>=10) & (lat<=80)]
    lon_r = lon[:]-180.
    time = nc.variables['time']

    z = nc.variables['z'][:, (lat>=10) & (lat<=80), : ]
    z = np.concatenate((z[:,:,int(len(lon_r)/2.):],z[:,:,:int(len(lon_r)/2.)]), axis = 2)

    time = num2date(time[:],time.units)
    
    # before creating the plots, check the existence of the figure directory
    if not(os.path.exists('/figures/entropy')):
        try:
            os.makedirs('./figures/entropy')
        except OSError:
            print ("Creation of the figures/entropy directory failed.")
        else:
            print ("Creation of the figures/entropy directory successfull.")


    # loop over time steps
    idx = -1
    for d in all_dates[:-16]:
        idx+=1
        print('plotting entropy on ', d)

        # compute average entropy field over 4 days (closeness snapshots are available at 00, 06, 12, 18, ie. 4 times per day, we thus need to take the average over 16 snapshots)
        image_number = "%03d" % idx
        file = './data/entropy-fields/H1out'+image_number+'.pckl'
        f = open(file, 'rb')
        e_data= pickle.load(f)
        f.close()
        for j in range(1,16): # add the next 15 time steps ()
            im_num = "%03d" % (idx+j)
            file = './data/entropy-fields/H1out'+im_num+'.pckl'
            f = open(file, 'rb')
            data_tmp= pickle.load(f)
            f.close()
            e_data+=data_tmp
        # average
        e_data = e_data.astype(np.float)/16.
        
       # interpolate 4-days-averaged closeness field for plotting (note that we are interested in latitudes 10-80ºN)
        xi = np.linspace(-180, 180, 100)
        yi = np.linspace(10, 81, 50)

        # Perform linear interpolation of the data on a grid defined by (xi,yi)
        triang = tri.Triangulation(cell_centers_x, cell_centers_y)
        interpolator = tri.LinearTriInterpolator(triang, e_data)
        Xi, Yi = np.meshgrid(xi, yi)
        e_interp = interpolator(Xi, Yi)

        # compute average geopotential height field over 4 days
        z_data = np.zeros((z.shape[1], z.shape[2]))
        dates_within_interval = all_dates[idx:idx+16]
        for dwi in dates_within_interval :
            ind = np.where(time == dwi)[0]
            z_data += z[ind,:,:][0]
        z_data/=16.
        
        # the minimum and maximum values of the plot color scale are fixed so that they are the same in all plots
        maximum = 2.3
        minimum = 1.3


        #############################################################
        # plot entropy
        viz(-175., 175., 15, 75, Xi, Yi, e_interp, minimum, maximum, lon_r, lat_r, z_data, str(all_dates[idx]), './figures/entropy/H1out'+image_number+'.pdf')
        


if __name__ == "__main__":
    main()
