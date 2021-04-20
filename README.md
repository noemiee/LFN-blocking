# LFN-blocking

Code and data for the results presented in:
Characteristic signatures of Northern Hemisphere blocking situations in Lagrangian network representations of atmospheric flow, N. Ehstand et al.

The workflow used to produce the figures is illustrated below. The ERA-Interim [1] data is fed into the model Flexpart [2] which outputs the particles' trajectories required to compute the Lagrangian Flow networks via the file "network.py" (see the methodology section of our paper for complementary information). Once the adjacency matrices of each of the networks are obtained, the measures (degree, entropy and harmonic closeness centrality) are computed and plotted using the "compute-\*.py" and "plot-\*.py" functions. 

![Figure1](https://github.com/noemiee/LFN-blocking/blob/main/Code-organization.png)

The Flexpart code, tutorials as well as information on the ERA-Interim data retrieval and pre-processing is available [here](https://www.flexpart.eu/). 

In this folder we provide the *network.py* file used to read and compute the Lagrangian Flow Networks. Then, the folder *2010* contains all the code and data necessary to reproduce the plots presented in the paper. The sub-folder structure is illustrated below (for the closeness computations). The color of the data folder indicates that it stores the input of the file of the same color. The 'dashed' contour for a folder indicated that the folder does not exist yet but will be created when running the code. 

In order to reproduce the results, the user must have installed python3 and simply follow the folowing steps:
1. uncompress the data present in the *gph* folder
2. run *compute-closeness.py*
3. run *plot-closeness.py*

![Figure2](https://github.com/noemiee/LFN-blocking/blob/main/closeness.png)



In addition, we also provide the networks adjacency matrices for June-July-August 2003 and 2018 (in the respective folders).

[1] Dee et al., “The ERA-interim reanalysis: configuration and performance of the data assimilation system” QuarterlyJournal of the Royal Meteorological Society 137, 553–597 (2011).

[2] Pisso et al., “The  lagrangian particle dispersion model FLEXPART version 10.4” Geoscientific Model Development 12, 4955–4997 (2019).
