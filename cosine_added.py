"""
Created on Fri Dec 15 14:13:12 2017

Copyright (c) 2017, Sam Lambrick 
All rights reserverd.
Subject to the MIT licence.

This script makes use of the elliptical_solid_angle module to calculate solid 
angles and uses the functions from that module to add a cosine distribution to 
model the proportion of an emitter of gas that passes through an elliptical 
appeture.

Again two functions are included, one that performs a double numerical integral
and the other that has had the theta integral performed analytically. The 
double integral is included for completness, and as a demonstartion of how the 
inclusion of distributions can be achieved.

The script calculates the solid angle and 'intensity' caught by the apeture for
a fixed perpendicular distance h for a range of parrallel distances p then plots
and saves the results using matplotlib.pyplot. The particular values or p,h,a,b
chosen are related to research work that I do in helium atom microscopy.
"""

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import elliptical_solid_angle as esa

# Use the ggplot style because I like it
plt.style.use("ggplot")

def intensity(p, h, a, b):
    """Calculates the intensity of a cosine distribution that enters an elliptical
    appeture.
    """
    
    if p < a:
        th_2 = lambda y : esa.theta_2(p, h, y, a, b)
        integrand = lambda y : np.cos(2*th_2(y))
        I = np.pi/2 - 0.25*integrate.quad(integrand, 0, 2*np.pi)[0]
    elif p == a:
        th_3 = lambda y : esa.theta_3(p, h, y, a, b)
        integrand = lambda y : np.cos(2*th_3(y))
        I = np.pi/4 - 0.5*integrate.quad(integrand, 0, np.pi/2)[0]
    else:
        phi_max = esa.phi_max_calc(p, a, b)
        th_4 = lambda y : esa.theta_4(p, h, y, a, b)
        th_5 = lambda y : esa.theta_5(p, h, y, a, b)
        integrand = lambda y : np.cos(2*th_4(y)) - np.cos(2*th_5(y))
        I = 0.5*integrate.quad(integrand, 0, phi_max)[0]
    
    return(I)

def intensity2(p, h, a, b):
    """Calculates the intensity of a cosine distribution that enters an elliptical
    appeture. Uses a double integral.
    """
    
    if p == 0:
        th_1 = lambda y : esa.theta_1(h, y, a, b)
        integrand = lambda y,z : np.sin(y)*np.cos(y)
        I = 2*integrate.dblquad(integrand, 0, np.pi, lambda x : 0, th_1)[0]
    elif p < a:
        th_2 = lambda y : esa.theta_2(p, h, y, a, b)
        integrand = lambda y,z : np.sin(y)*np.cos(y)
        I = 2*integrate.dblquad(integrand, 0, np.pi, lambda x : 0, th_2)[0]
    elif p == a:
        th_3 = lambda y: esa.theta_3(p, h, y, a, b)
        integrand = lambda y,z : 2*np.sin(y)*np.cos(y)
        I = integrate.dblquad(integrand, 0, np.pi/2, lambda x : 0, th_3)[0]
    else:
        phi_max = esa.phi_max_calc(p, a, b)
        th_4 = lambda y : esa.theta_4(p, h, y, a, b)
        th_5 = lambda y : esa.theta_5(p, h, y, a, b)
        integrand = lambda y,z : 2*np.sin(y)*np.cos(y)
        I = integrate.dblquad(integrand, 0, phi_max, th_4, th_5)[0]
    
    return(I)

if __name__ == "__main__": 
    # Constants
    aa = np.sqrt(2)/2
    bb = 0.5
    hh = 2.121
    
    # Number of points to plot
    n = 51
    
    # Variables
    ps = np.linspace(0, 5, num=n)
    omegas = np.zeros(n)
    omegas2 = np.zeros(n)
    intensities = np.zeros(n)
    intensities2 = np.zeros(n)
    
    # Calculations
    for i in range(len(ps)):
        omegas[i] = esa.solid_angle_calc(ps[i], hh, aa, bb)
        intensities[i] = intensity(ps[i], hh, aa, bb)
        omegas2[i] = esa.solid_angle_calc2(ps[i], hh, aa, bb)
        intensities2[i] = intensity2(ps[i], hh, aa, bb)
    
    # Normailse the two sets of data so that the values are the proportion of
    # the distribution (uniform or cosine) that from the hemisphere enter the 
    # ellipse.
    omegas = omegas/(2*np.pi)
    omegas2 = omegas2/(2*np.pi)
    intensities = intensities/np.pi
    intensities2 = intensities2/np.pi
    
    # Plotting
    plt.plot(ps, omegas2, '.', label="Double solid angle")
    plt.plot(ps, intensities2, '.', label="Double cosine distrbution")
    plt.plot(ps, omegas, label="Solid angle")
    plt.plot(ps, intensities, label="Cosine distrbution")
    plt.xlabel("Position form the ellipse centre")
    plt.ylabel("Proportion into ellipse")
    plt.legend(loc="upper right")
    
    plt.savefig("simple_ellipse_plot.eps")


