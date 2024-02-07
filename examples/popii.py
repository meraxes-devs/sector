import sector as mc
import sys, os
from dragons import meraxes # Use functions from the dragons package
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

meraxes.io.set_little_h(cosmo['h'])
target_snap = int(sys.argv[1]) # The snapshot of your interested target
target_index = int(sys.argv[2]) # The index of your interested target

fmeraxes = '/DIRECTORY/meraxes.hdf5'
save_dir = './'
redshift = meraxes.io.read_snaplist(fmeraxes,h=cosmo['h'])[1][target_snap] # get the redshift

# save SFH
mc.save_star_formation_history(fmeraxes, [target_snap,], [target_index,],
                cosmo['h'], prefix='sfh_snap', outPath=save_dir)

# get the spectrum (dust-free)
mc.composite_spectra(fmeraxes, snapList=[target_snap,], gals=[save_dir+'/sfh_snap_%03d.bin'%(target_snap),],
                h=cosmo['h'], Om0=cosmo['omega_M_0'],
                sedPath=mc.STARBURST99_Kroupa, # we also provide STARBURST99_Salpeter
                outType="sp", obsFrame=True,
                betaBands = [], restBands = [[1600, 100],], obsBands = [],
                prefix='spectra_snap', nThread=1,
                outPath=save_dir)


# calculate dust factors 
gals =meraxes.io.read_gals(fmeraxes,sim_props=False,indices=target_index,
                                    props=("ColdGas", "MetalsColdGas","Spin","Rvir"), 
                                    snapshot=target_snap,quiet=True)

Z = gals['MetalsColdGas'] /  gals['ColdGas']
Z[np.isnan(Z)] = 0
Z[Z<0] = 0
Z[Z>1] = 1

factor = (Z / 0.02)**0.65 * gals['ColdGas'] * cosmo['h'] * (gals['Spin'] * gals['Rvir'] * cosmo['h']/ np.sqrt(2) * 1e3) ** (-2.0) * np.exp(-0.35 * redshift)

# get the spectrum (dust-attenuated)
mc.composite_spectra(fmeraxes, snapList=[target_snap,], gals=[save_dir+'/sfh_snap_%03d.bin'%(target_snap),],
                        h=cosmo['h'], Om0=cosmo['omega_M_0'],
                        sedPath=mc.STARBURST99_Kroupa,
                        dust=np.array([13.5 * factor,
                            np.zeros_like(factor)-1.6, 
                            381.3 * factor, 
                            np.zeros_like(factor)-1.6, 
                            np.zeros_like(factor)+10000000]).T.reshape([1, len(target_index), 5]),
                        outType="sp", obsFrame=True,
                        betaBands = [], restBands = [[1600, 100],], obsBands = [],
                        prefix='spectra_dust_snap', nThread=1, 
                        outPath=save_dir)
