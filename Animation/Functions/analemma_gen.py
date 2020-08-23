import numpy as np
#import anomaly # This can be added if you want to access functions to convert between anomalies
from anomaly import true_anomaly

def time_array(n,period,lon,time_obs=1200,lon0=0,mean_solar_day=24):
    #n=no. of points, period = period of revolution, lon = longitude, time_obs = 24 hr format
    #lon0 = longitude facing earth at periapsis
    t0 = (int(time_obs/100) + (time_obs%100)/60)/mean_solar_day
    delta_t = (lon0-lon)/(2*np.pi) + (t0-0.5)
    t = np.linspace(0,period,n) + delta_t
    return t

def get_dec(f,i,omega=np.radians(85.901)):
    #f = true anomaly, i = inclination, omega = argument of periapsis)
    dec = np.arcsin(np.sin(i)*np.cos(f+np.pi/2-omega))
    return dec

def get_ra(f,i,omega=np.radians(85.901)):
    #f = true anomaly, i = inclination, omega = argument of periapsis)
    ra = np.arccos(np.sin(f+np.pi/2-omega)/np.cos(np.arcsin(np.sin(i)*np.cos(f+np.pi/2-omega))))
    #Following make ra in the proper range of [0,2pi]
    mask = np.ones(np.shape(ra), dtype = bool)
    case = mask * ((f+np.pi/2-omega<=np.pi/2)|(f+np.pi/2-omega>3*np.pi/2))
    ra[case] = (2*np.pi) - ra[case]
    return ra

def hour_angle(period,ra,t):
    #Period of revolution, RA and time array
    ra_obs = (np.modf(t*(1+1/period))[0])*2*np.pi - np.pi/2
    ha = (ra_obs - ra)%(2*np.pi)
    return ha

def analemma_hour_angle(period,ra,t,omega=np.radians(85.901)):
    #omega = argument of periapsis
    mean_anomaly = (2*np.pi*t/period) % (2*np.pi)
    ra_hyp_obs = (mean_anomaly - omega) % (2*np.pi)
    ana_ha = ra_hyp_obs - ra
    return ana_ha

def altitutde(ha,dec,lat):
    #Hour angle, declination, Latitude
    alt = np.arcsin(np.cos(lat)*np.cos(dec)*np.cos(ha)+np.sin(lat)*np.sin(dec))
    return alt

def azimuth(dec,alt,ha,lat):
    #Declination, altitude, hour angle, latitude
    az = (np.arctan2(np.cos(dec)*np.sin(ha)/np.cos(alt),((np.sin(lat)*np.cos(dec)*np.cos(ha))-(np.cos(lat)*np.sin(dec)))/np.cos(alt)))%(2*np.pi)
    return az

def sph_cart(r,theta,phi):
    #r=radius,phi=angle with X-axis, theta=Angle with Z-axis
    X=r*np.sin(theta)*np.cos(phi)
    Y=r*np.sin(theta)*np.sin(phi)
    Z=r*np.cos(theta)
    return [X,Y,Z]

def analemma_gen(e,i,period,n,lon,lat,omega = 85.901,time_obs=1200,mean_solar_day = 24,lon0 = 0,mode='altaz',r_cart=1):
    #e = eccentricity, i = inclination, period = period of revolution, n = number of data points in time array,
    #lon = longitude, lat = latitude, omega = argument of periapsis, time_obs = time of observance in hhmm format,
    #mean_solar_day in hours, lon0 = longitude facing sun at periapsis, mode = unit in which you want the data points
    if (type(e)==float)|(type(e)==int):
        e = np.array([e])
    if (type(i)==float)|(type(i)==int):
        i = np.array([i])
    
    e = e.reshape((len(e),1))
    i = i.reshape((len(i),1,1))
    omega = np.radians(omega)
    period = period #days
    n = n #no. of values
    lon0 = np.radians(lon0)
    lon = np.radians(lon)
    lat = np.radians(lat)
    incl = np.radians(i)
    time_of_obs = time_obs
    
    t = time_array(n,period,lon,time_of_obs,lon0,mean_solar_day)
    ecc = np.ravel(e)
    f = np.array([])
    for item in ecc:
        f = np.append(f,true_anomaly(item,period,t))
    f = f.reshape((len(e),n))
    r = get_ra(f,incl,omega)
    d = get_dec(f,incl,omega)
    ha = analemma_hour_angle(period,r,t,omega)
    alt = altitutde(ha,d,lat)
    az = azimuth(d,alt,ha,lat)
    cart = sph_cart(r_cart,np.pi/2 - alt,az)
    if mode=='radec':
        return [r,d]
    elif mode=='altaz':
        return [alt,az]
    elif mode=='cartesian':
        return cart
    else:
        raise Exception('Please enter correct mode')
    
if __name__ == '__main__':    
	#In cartesian it will give coordinates with radius 1 unit.
	print(analemma_gen(e = np.linspace(0,1,10,endpoint=False),
	                i = np.linspace(0,90,10,endpoint=False),
	                period = 365.2422,                          #no. of mean solar days
	                n = 365,                                    #no. of points
	                lon = 0,                                    #deg
	                lat = 0,                                    #deg
	                omega = 85.901,                             #argument of periapsis
	                time_obs = 1200,                            #24 hr format
	                mean_solar_day = 24,                        #hours
	                lon0 = 0,                                   #deg
	                mode = 'cartesian',							#'radec' / 'altaz' / 'cartesian'
	                r_cart = 10))                         			# radius of sphere for cartesian only
