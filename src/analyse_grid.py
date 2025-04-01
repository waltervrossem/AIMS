#!/usr/bin/env python

"""
A simple utility for analysing binary AIMS grids
"""

import sys
import dill
import numpy as np
import matplotlib.pyplot as plt

iage = 0  

def print_number_of_models_in_tracks(grid):
    for i in range(len(grid.tracks)):
        print("Track %d contains %d models"%(i,len(grid.tracks[i].names)))

def print_number_of_modes_in_models(grid):
    for i in range(len(grid.tracks)):
        for j in range(len(grid.tracks[i].names)):
            print("Track %d %s %d"%(i,grid.tracks[i].names[j], \
                grid.tracks[i].mode_indices[j+1]-grid.tracks[i].mode_indices[j]))

def plot_number_of_modes(grid):
    nmodes = []
    for track in grid.tracks:
        for i in range(1,len(track.mode_indices)):
            nmodes.append(track.mode_indices[i]-track.mode_indices[i-1])
    nmax = max(nmodes)
    x = np.array(range(nmax+1))
    y = np.zeros((nmax+1,),dtype=np.int)
    for i in x: y[i] = nmodes.count(i)
    plt.plot(x,y,"b-")
    plt.plot(x,y,"bo")
    plt.xlabel(r"Number of modes in a model")
    plt.ylabel(r"Number of models with N modes")
    plt.savefig("nmodes.png")
    plt.cla()
    plt.clf()
    plt.close()

def find_range(grid,i):
    param_values = np.sort(np.unique(np.array([track.params[i] for track in grid.tracks])))
    print("Range on parameter %s: %8.5e to %8.5e, %d values"%(grid.grid_params[i], \
        np.min(param_values),np.max(param_values),len(param_values)))

def find_param_values(grid,i):
    param_values = np.sort(np.unique(np.array([track.params[i] for track in grid.tracks])))
    print("Values for parameter %s: %s"%(grid.grid_params[i],str(param_values)))

def duplicate_ages(track):
    """
    Check to see if you track contains models with duplicate ages.

    :return: ``True`` if there are duplicate age(s)
    :rtype: boolean

    .. warning::
        This method should only be applied after the track has been
        sorted.
    """

    return any(track.glb[i,iage] == track.glb[i+1,iage] for i in range(len(track.names)-1))

def save_track(track,ndx):
    filename = "track%d.txt" % ndx
    with open(filename, "w") as f:
        f.write("%s\n"%(str(track.params)))
        for i in range(len(track.names)):
            f.write("*"*60+"\n")
            f.write("%s\n"%(track.names[i]))
            f.write("Global parameters\n")
            f.write("   %s\n"%(str(track.glb[i])))
            f.write("Modes:\n")
            for j in range(track.mode_indices[i],track.mode_indices[i+1]):
                f.write("   %s\n"%(str(track.modes[j])))

def check_tracks_sorted(grid):
    for i in range(len(grid.tracks)):
        glb = grid.tracks[i].glb
        if (all(glb[i,0] <= glb[i+1,0] for i in range(glb.shape[0]-1))):
            print("Track %d is sorted"%(i))
        else:
            print("Track %d is not sorted"%(i))

if __name__ == "__main__":

    # check number of arguments
    assert (len(sys.argv) > 1), "Usage: test_grid.py binary_grid_file"

    grid = dill.load(open(sys.argv[1],"rb"))
    ntracks = len(grid.tracks)

    print("Model postfix:    "+str(grid.postfix))
    print("Grid parameters:  "+str(grid.grid_params))
    print("User parameters:  "+str(grid.user_params))
    print("Number of dims.:  "+str(grid.ndim+1))
    print("Number of tracks: "+str(ntracks))
    print("Number of models: "+str(sum([len(track.names) for track in grid.tracks])))
    print("Number of modes:  "+str(sum([track.nmodes for track in grid.tracks])))

    nmax = max([len(track.names) for track in grid.tracks])
    print("Number of models in smallest track: "+str(min([len(track.names) for track in grid.tracks])))
    print("Number of models in largest track:  "+str(nmax))

    for track in grid.tracks:
        if (duplicate_ages(track)):
            print("ERROR: the track %s = %s"%(track.grid_params,track.params))
            print("       has models with the same age.  Please remove")
            print("       duplicate models.")
            #for i in range(len(track.models)):
            #    print("Model[%d]: %e"%(i,track.models[i].glb[iage]))

    for i in range(len(grid.grid_params)):
        find_range(grid,i)
        #find_param_values(grid,i)

    hist = np.zeros((nmax+1,))
    for track in grid.tracks: hist[len(track.names)] += 1
    plt.plot(list(range(nmax+1)),hist,"b-")
    plt.xlabel(r"Number of models in a track, N")
    plt.ylabel(r"Number of tracks with N models")
    #plt.show()
    plt.savefig("histogram.png")
    plt.cla()
    plt.clf()
    plt.close()

    plot_number_of_modes(grid)
    #check_tracks_sorted(grid)
    #for i in range(len(grid.tracks)): save_track(grid.tracks[i],i)
    #print_number_of_models_in_tracks(grid)
    #print_number_of_modes_in_models(grid)
