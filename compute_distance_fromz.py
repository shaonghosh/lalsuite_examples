import numpy as np

c = 299792458.0
def get_lum_dist(z, N=1e6, H0=67.9, OmegaM=0.3065, OmegaL=0.6935):
    """
    Given redshift value this function computes the luminosity
    distance by simple application of integral by Riemann sum
    """
    z_grid = np.linspace(0, z, int(N))
    dz = np.mean(np.diff(z_grid))
    I = 1/np.sqrt(OmegaM*(1 + z_grid)**3 + OmegaL)
    S = np.sum(I*dz)
    return (c/1000)*(1+z)*S/H0

def get_comoving_dist(z, N=1e6, H0=67.9, OmegaM=0.3065, OmegaL=0.6935):
    """
    Given redshift value this function computes the comoving 
    distance by simple application of integral by Riemann sum
    Ref: https://ned.ipac.caltech.edu/level5/Hogg/Hogg4.html
    """
    d_H = c/H0/1000
    z_grid = np.linspace(0, z, int(N))
    dz = np.mean(np.diff(z_grid))
    E = np.sqrt(OmegaM*(1 + z_grid)**3 + OmegaL)
    Dc = d_H*np.sum(dz/E)
    return Dc

def get_redshift(D_L, tol=0.001, H0=67.9, OmegaM=0.3065, OmegaL=0.6935):
    """
    Given a luminosity distance this function computes the redshift
    using method of bisection.
    """
    z1 = 1
    z2 = 10
    D1 = get_lum_dist(z1, N=1e6, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
    D2 = get_lum_dist(z2, N=1e6, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
    while np.sign(D_L - D1) == np.sign(D_L - D2):
        z1 /= 2
        z2 *=2
        D1 = get_lum_dist(z1, N=1e6, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
        D2 = get_lum_dist(z2, N=1e6, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)

    z = (z1 + z2)/2
    D = get_lum_dist(z, N=1e6, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
    while np.abs(D_L - D) > tol:
        if D < D_L:
            z1 = z
        elif D > D_L:
            z2 = z
        z = (z1 + z2)/2
        D = get_lum_dist(z, N=1e6, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)

    return z

def get_redshifted_mass(m, D_L, H0=67.9, OmegaM=0.3065, OmegaL=0.6935):
    """
    For a source at a given luminosity distance this function computes
    the redshifted mass value
    """
    z = get_redshift(D_L, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
    return (1 + z)*m




