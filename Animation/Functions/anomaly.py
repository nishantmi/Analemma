import numpy as np
from scipy.optimize import fsolve

def mean2eccentric(mean,eccentricity):      #mean is the array of mean anomalies to be converted to eccentric anomaly
    def mean_eccentric_relation(eccentric):
        # The following realtion between is satisfied when dummy = 0
        dummy = eccentric - eccentricity*np.sin(eccentric)- mean
        return dummy
    eccentric_guess = np.zeros((1,len(mean)))
    return fsolve(mean_eccentric_relation,eccentric_guess)

def mean2true(mean,eccentricity):           #mean is the array of mean anomalies to be converted to true anomaly
    e = eccentricity
    E = mean2eccentric(mean,e)
    return 2*np.arctan((((1+e)/(1-e))**0.5)*np.tan(E/2))

def eccentric2mean(eccentric,eccentricity): #eccentric is the array of eccentric anomalies to be converted to mean anomaly
    return eccentric - eccentricity*np.sin(eccentric)

def eccentric2true(eccentric,eccentricity): #eccentric is the array of eccentric anomalies to be converted to true anomaly
    E = eccentric
    e = eccentricity
    f = 2*np.arctan((((1+e)/(1-e))**0.5)*np.tan(E/2))
    f = (f + 2*np.pi)%(2*np.pi)
    return f

def true2eccentric(true, eccentricity):     #true is the array of true anomalies to be converted to eccentric anomaly
    e = eccentricity
    E = 2*np.arctan((((1-e)/(1+e))**0.5)*np.tan(true/2))
    E = (E+2*np.pi)%(2*np.pi)
    return E

def true2mean(true, eccentricity):          #true is the array of true anomalies to be converted to mean anomaly
    e = eccentricity
    E = true2eccentric(true,e)
    M = eccentric2mean(E,e)
    return M

def true_anomaly(eccentricity,period,t):    #period = period of revolution, t = time array
    e = eccentricity
    M = 2*np.pi*t/period
    E = mean2eccentric(M,e)
    return eccentric2true(E,e)
