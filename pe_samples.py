import h5py
import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sn

from scipy.interpolate import interp1d
from astropy import cosmology, units as u


'''
This scripts have some helper functions that can be used to read posterior samples files 
from PE runs. The expected file-type is HDF5. 
'''

def getMasses(q, mc):
    """
    Given chirp mass and mass-ratio, this function finds the omponent masses
    """
    m1 = mc * (1 + q)**(1/5) * (q)**(-3/5)
    m2 = mc * (1 + q)**(1/5) * (q)**(2/5)
    return (m1, m2)
   
   
def read(posterior_samples_file):
    """
    This function accepts posterior samples file in HDF5 format and returns a pandas
    DataFrame for the samples. This helps in readability of the samples when printed
    in the command line, Jupyter Notebook.
    """
    data = h5py.File(posterior_samples_file)
    engine = list(data['lalinference'].keys())
    samples = data['lalinference'][engine[0]]['posterior_samples'][:]
    samples_dataframe = pd.DataFrame(samples)
    return samples_dataframe
    
    
def get_redshifts(posterior_samples_file=None, distances=None, N=10000):
    """
    This function computes the redshifts using the Planck15 cosmology.
    It accepts as input the posterior samples from a PE file and then
    computes redshift by interpolating the distance-redshift function.

    Parameters
    ----------
    distances: float or numpy.ndarray
              distance(s) in Mpc

    N: int
      Number of steps for the computation of the interpolation function

    Example
    -------
    >>> pe_samples.get_redshifts(posterior_samples_file)
    array([0.00225566, 0.00450357, 0.00674384, 0.00897655,
           0.01120181, 0.0134197 , 0.01563032, 0.01783375
           0.02003009, 0.02221941])
    """
    if posterior_samples_file:
        df = read(posterior_samples_file)
        distances = df['dist'].values
    function = cosmology.Planck15.luminosity_distance
    min_dist = np.min(distances)
    max_dist = np.max(distances)
    z_min = cosmology.z_at_value(func=function, fval=min_dist*u.Mpc)
    z_max = cosmology.z_at_value(func=function, fval=max_dist*u.Mpc)
    z_steps = np.linspace(z_min - (0.1*z_min), z_max + (0.1*z_min), N)
    lum_dists = cosmology.Planck15.luminosity_distance(z_steps)
    s = interp1d(lum_dists, z_steps)
    redshifts = s(distances)
    return redshifts
    

def massdist(posterior_samples_file):
    """
    This function plots the mass distributions given the posterior 
    samples file. It creates the plots for the masses in detector
    frame and the source frame.
    """
    df = read(posterior_samples_file)
    distances = df['dist'].values
    z = get_redshifts(distances=distances)
    mc = df['mc'].values
    mc_source = mc/(1 + z)
    q = df['q'].values
    m1, m2 = getMasses(q, mc)
    m1_source, m2_source = getMasses(q, mc_source)
    
    pl.figure(1)
    pl.hist(mc, 100, density=True, alpha=0.3)
    sn.kdeplot(mc, shade=True, label='Detector-Frame')
    pl.hist(mc_source, 100, density=True, alpha=0.3)
    sn.kdeplot(mc_source, shade=True, label='Source-Frame')
    pl.legend()
    pl.xlabel('$\\mathcal{M}~(M_{\\odot})$')

    pl.figure(2)
    pl.subplot(1,2,1)
    sn.kdeplot(m1, m2, shade=True, label='Detector-Frame')
    pl.scatter(m1, m2, marker='.', s=1, alpha=0.4)
#     pl.legend()
    pl.title('Detector-Frame')
    pl.xlabel('$m_1~(M_{\\odot})$')
    pl.xlabel('$m_2~(M_{\\odot})$')

    pl.subplot(1,2,2)
    sn.kdeplot(m1_source, m2_source, shade=True, label='Source-Frame')
    pl.scatter(m1_source, m2_source, marker='.', s=1, alpha=0.4)
#     pl.legend()
    pl.title('Source-Frame')
    pl.xlabel('$m_1~(M_{\\odot})$')
    
    
    








