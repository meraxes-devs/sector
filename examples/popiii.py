import sector as mc
import sys, os
import h5py # read Meraxes' output directly
import numpy as np

# cosmology to be consistent with the model (i.e. SAM and N-body)
# though they can also be found in Meraxes' output
cosmo = {'omega_M_0' : 0.3121, 
'omega_lambda_0' : 0.6879, 
'omega_b_0' : 0.049043394, 
'omega_n_0' : 0.0,
'N_nu' : 0,
'h' : 0.6751,
'n' : 0.968,
'sigma_8' : 0.815
}

target_snap = int(sys.argv[1]) # The snapshot of your interested target
target_index = int(sys.argv[2]) # The index of your interested target
sedPath = sys.argv[3] # The sed library you would like to use
# How much time in units of Myr to observe gals after the beginning 
# of the snapshot. If negative, put it at the end of the snapshot
DeltaT = float(sys.argv[4]) 

fmeraxes = '/DIRECTORY/meraxes.hdf5'
save_dir = '/DIRECTORY/'
population = 3 # specify that you are assuming the entire SF history is popIII (can be inconsistent with Meraxes!)

# save SFH
mc.save_star_formation_history(fmeraxes, target_snap, [target_index,],
                cosmo['h'], prefix='sfh_index%d_snap'%target_index,outPath=save_dir, population=population)

# get the spectrum (dust-free because they are metal free)
mc.composite_spectra(fmeraxes, snapList=[target_snap,] , gals=[save_dir+'/sfh_index%d_snap_%03d.bin'%(target_index, target_snap),],
                h=cosmo['h'], Om0=cosmo['omega_M_0'],DeltaT=DeltaT, population=population,
                sedPath=sedPath,
                outType="sp", obsFrame=False, IGM=False,
                betaBands = [], restBands = [[1600, 100],], obsBands = [],
                prefix='/spectra_DeltaT%.3fMyr_index%d_snap'%(DeltaT, target_index),
                outPath=save_dir)

# get the photometry
mc.composite_spectra(fmeraxes, snapList=[target_snap,] , gals=[save_dir+'/sfh_index%d_snap_%03d.bin'%(target_index, target_snap),],
                h=cosmo['h'], Om0=cosmo['omega_M_0'],DeltaT=DeltaT, population=population,
                sedPath='/home/yqin/bitbucket/meraxes_manu2/input/photometric_tables/%s.hdf5'%sed,
                outType="ph", obsFrame=False, IGM=False,
                betaBands = [], restBands = [[1600, 100],], obsBands = [],
                prefix=sed+'/photometry_DeltaT%.3fMyr_index%d_snap'%(DeltaT, target_index),
                outPath=save_dir)

