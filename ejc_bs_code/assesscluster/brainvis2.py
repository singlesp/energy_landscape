print('start')
import os
import numpy as np
import scipy.io as sio
import nibabel as nib
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import sys

import mayavi as my
import surfer

##############################
### Get brainmapping files ###
##############################

name_root = sys.argv[1]
numClusters = int(sys.argv[2])
#lausanneScaleBOLD = int(sys.argv[3])
nparc = int(sys.argv[3])
scan = sys.argv[4]
scan_ind = int(sys.argv[5])
homedir = "/Users/sps253"

# index in humanRegionNames.mat corresponding to parcellation scale in lausanneScaleBOLD

#if lausanneScaleBOLD == 60:
#    parc_ind = 1
#if lausanneScaleBOLD == 125:
#    parc_ind = 2
#if lausanneScaleBOLD == 250:
#    parc_ind = 3

# which scan to plot centroids from
if scan == 'R':
    scanlab = 'Rest'
if scan == 'N':
    scanlab = 'nBack'
if scan == 'C':
    scanlab = ['RestComb','nBackComb','Overall']
    scanlab = scanlab[scan_ind-1]

os.environ['SUBJECTS_DIR']='/Applications/freesurfer/7.1.1/subjects'

sys.path.append(homedir + "/Documents/brain_states-master")
from brain_states-master import *

namefile = homedir + "/Dropbox/Cornblath_Bassett_Projects/code/control_fc/restnbackpipeline/analysiscode/human_regionNames.mat"
roinames = sio.loadmat(namefile)['roinames'][0][parc_ind]
nparc = len(roinames)
for i in np.arange(0,nparc):
    roinames[i] = str(roinames[i][0][0])

##################
### Set values ###
##################

views=['lat','med']
cbar=False
clim = 0.5
thrsh = -100
subject_id = 'fsaverage'
surf = 'pial'
hemis = ['lh','rh']
adv_view = {'lh': {'azimuth': 45, 'elevation': 79},'rh': {'azimuth': 135, 'elevation': 79}}

savedir = homedir + "/Dropbox/Cornblath_Bassett_Projects/BrainStateTransitions/brain_states/results/" + \
    name_root + '/centroids/k' + str(numClusters) +'/'
if not os.path.exists(savedir):
    os.makedirs(savedir)

#################
### Load data ###
#################


clusterfile = homedir + "/Dropbox/Cornblath_Bassett_Projects/BrainStateTransitions/brain_states/results/" + \
 name_root + '/transitionprobabilities/' + scanlab + "ClusterCentroids_k" + str(numClusters) + name_root + ".mat"
data = sio.loadmat(clusterfile)
centroids = data['kClusterCentroids']
clusterNames = [x[0][0] for x in data['clusterNames']]
clusterNamesUp = [x[0][0] for x in data['clusterNamesUp']]
clusterNamesDown = [x[0][0] for x in data['clusterNamesDown']]

# yeo 
yeo = sio.loadmat(homedir + "/Dropbox/Cornblath_Bassett_Projects/BrainStateTransitions/brain_states/data/yeo7netlabelsLaus250.mat")
np.array([1,1,7,7,7])
centroids = [(yeo['network7labels'] == x).astype(float) for x in np.array([1,1,7,7,7])]
centroids = np.squeeze(np.asarray(centroids)).T[:nparc,:]
centroids[centroids == 0] = -1
centroids[:,numClusters - 1] = np.ones(nparc)
###################
### Plot brains ###
###################

os.chdir(savedir)
for hemi in hemis:
    aparc_file = homedir + "/Dropbox/Cornblath_Bassett_Projects/code/brainmapping2/el_fsavg/" + hemi + ".myaparc_" + str(lausanneScaleBOLD) + ".annot"          
    for K in np.arange(numClusters):
        #centroid = stats.zscore(centroids[:,K])
        centroid = centroids[:,K]
        centroid[np.logical_and(centroid < 0.5,centroid > -0.5)] = -101
        for view in views:  
            vtx_state = getvertdata_lausanne(vtx_data = np.expand_dims(centroid,axis = 1), roinames = roinames, aparc_file = aparc_file, hemi = hemi)
            fig = my.mlab.figure(size=(340,340))
            fig = my.mlab.gcf()
            brain = surfer.Brain(subject_id, hemi, surf,figure=fig,views=view,background='white', alpha = 1)
            brain.add_data(vtx_state, min = -clim, max = clim, thresh = thrsh, colormap="plasma", alpha=.8,colorbar=cbar)       
            #if view == 'lat':
            #    brain.show_view(adv_view[hemi])
            fname = clusterNames[K] + view + hemi + str(K) + '.png'
            my.mlab.savefig(figure=fig,filename=fname,size=(1,1))
            my.mlab.close()

#Arrange brains into grid

plt.figure(figsize = [numClusters,2])
arial = {'fontname':'Arial'}

for K in np.arange(numClusters):
    plt.figure(figsize = [1.7,1.7])
    ttl = clusterNames[K] + '\n ' + clusterNamesUp[K] + ' / ' + clusterNamesDown[K]
    plt.suptitle(ttl,fontsize=8,fontweight='bold',**arial)        
    for H,hemi in enumerate(hemis):
        for V, view in enumerate(views):    #Arrange brains into grid
            fname = clusterNames[K] + view + hemi + str(K) + '.png'
            img = mpimg.imread(fname)
            plt.subplot(2,2,(H + 2*V+1))
            imgplot = plt.imshow(img,aspect='auto')
            plt.axis('off')

    fname = clusterNames[K] + str(lausanneScaleBOLD) + "ThreshYeo_k" + str(numClusters) + ".png"
    plt.subplots_adjust(hspace=0, wspace=0.05)
    #plt.tight_layout()
    plt.savefig(fname,dpi=500,bbox_inches='tight',pad_inches=0.1)
    print(['saved',fname])
"""
for V in np.arange(1,len(views)+1):
    for K in np.arange(1,numClusters+1):
        fname = clusterNames[K-1] + views[V-1] + hemi + str(K-1) + '.png'
        img = mpimg.imread(fname)
        plt.subplot(2,numClusters,K + numClusters*(V-1))
        imgplot = plt.imshow(img,aspect='auto')
        if V == 1:
            ttl = clusterNames[K-1] + '\n ' + clusterNamesUp[K-1] + ' / ' + clusterNamesDown[K-1]
            plt.title(ttl,fontsize=8,fontweight='bold',**arial)
        plt.axis('off')

fname = scanlab + hemi + "ClusterCentroids_k" + str(numClusters) + name_root + ".png"
plt.subplots_adjust(hspace=0, wspace=0.05)
#plt.tight_layout()
plt.savefig(fname,dpi=500,bbox_inches='tight',pad_inches=0)
print(['saved',fname])
"""
#delete intermediate files
for V in np.arange(1,len(views)+1):
    for K in np.arange(1,numClusters+1):
        fname = clusterNames[K-1] + views[V-1] + hemi + str(K-1) +  '.png'
        os.remove(fname)
print('Deleted intermediate files')
