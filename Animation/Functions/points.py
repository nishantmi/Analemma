import numpy as np
from analemma_gen import sph_cart

# Returns multiple numbers, precision too low of data
# Depricated function
#def peri_apop(data):
#    # Assuming data to be in shape (365,3)
#    r = np.sqrt(data[:, 0]**2 + data[:, 1]**2 + data[:, 2]**2) # Distance from sun at all points
#
#    apopapsis = np.where(r==r.min())
#    apopapsis = np.squeeze(apopapsis)
#
#    periapsis = np.where(r==r.max())
#    periapsis = np.squeeze(periapsis)
#    return periapsis,apopapsis
##    return data[periapsis],data[apopapsis]

def equinox(radec):
    # assuming dec as vectors
    dec = radec[:,1]
    ra = radec[:,0]
    positive = dec[np.where(dec>0)]
    negative = dec[np.where(dec<0)]

    eq1 = np.where(positive == positive.min()) # The two data points closest to zero
    eq2 = np.where(negative == negative.max())
    ra1 = ra[np.where(dec==positive[eq1])]
    ra1 = np.squeeze(ra1)
    ra2 = ra[np.where(dec==negative[eq2])]
    ra2 = np.squeeze(ra2)

    c1 = np.array(sph_cart(20,np.pi/2,ra1)).T
    c2 = np.array(sph_cart(20,np.pi/2,ra2)).T
    return c1,c2 # return type cartesian

def solstice(radec):
    #  Dec as vector
    dec = radec[:,1]
    ra = radec[:,0]
    sol1 = ra[np.where(dec==dec.max())]
    sol2 = ra[np.where(dec==dec.min())]
    cs1 = np.array(sph_cart(20,np.pi/2,sol1)).T
    cs2 = np.array(sph_cart(20,np.pi/2,sol2)).T
    return cs1,cs2 # return type cartesian

if __name__ == '__main__':
#    Minor Testing script
    a = np.array([3,4,5,6,7,8,-3,-4,-2])
    p,a = equinox(a)
