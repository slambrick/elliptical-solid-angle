"""
Created on Fri Dec  8 16:58:00 2017

Copyright (c) 2017, Sam Lambrick 
All rights reserverd.
Subject to the MIT licence.

Contains the nessacery functions to perform on-axis analytical solid angle
integrals for ellipses. Also contains a function for integrating and calculating
the solid angle.

The derivation of the functions was performed by Abbas et al. 2015, for an 
understanding of the equations used it is recommened to look there.
doi:10.1016/j.nima.2014.10.061

The theta integral from Abbas et al. can be performed analytically. This has 
been used here along with the double integral form.

Variables used here, using the same naming convention as used by Abbas:
    p - the distance between the centre of the ellipse and the point of interest
        in the plane of the ellipse
    a - the semi-axis of the ellipse along which the the point of interest lies
    b - the semi-axis of the ellipse perpendicula
    h - the vertical perpendicular distance between the point of interest and 
        the plane of the ellpise

Solid angle functions:
    solid_angle_calc  - performs a single integration to calculate the solid 
                        angle
    solid_angle_calc2 - performs a double inetegration, using the form directly
                        from Abbas et al., to calculate the solid angles

Neither of these functions return the error estimate on the integrals, however 
it is easy to modify them to do so.

This module makes use of numeric python (numpy) and scientific python (scipy). 
"""

import numpy as np
import scipy.integrate as integrate

def sec(theta):
    """Calculates the sectant, 1/cos(theta)."""
    return(1/np.cos(theta))

def r_1(phi, a, b):
    """Calculates the distance from the point of interest to the edge of the
    ellipse, in the plane of the ellipse. For the case p = 0."""
    
    denominator = (a*np.sin(phi))**2 + (b*np.cos(phi))**2
    result = a*b/np.sqrt(denominator)
    return(result)

def r_2(p, phi, a, b):
    """Calculated the distance from the point of interest to the edge of the 
    ellipse, in the plane of the ellipse. For the case p < a."""
    
    numerator1 = -p*b*b*np.cos(phi)
    numerator2 = a*b*np.sqrt((b*np.cos(phi))**2 + (a*a - p*p)*(np.sin(phi))**2)
    denominator = (a*np.sin(phi))**2 + (b*np.cos(phi))**2
    result = (numerator1 + numerator2)/denominator
    return(result)

def r_3(p, phi, a, b):
    """Calculates the distance from the point of interest to the edge of the 
    ellipse, in the plane of the ellipse. For the case of p = a."""
    
    numerator = 2*a*b*b*np.cos(phi)
    denominator = (a*np.sin(phi))**2 + (b*np.cos(phi))**2
    result = numerator/denominator
    return(result)

## This is the negative of what it is in the paper
def r_4(p, phi, a, b):
    """Calculates the distance from the point of interest to the edge of the
    ellipse projected in the plane of the ellipse. For the case p > a."""
    
    numerator1 = a*a*p* (np.sin(phi))**2 * sec(phi)
    numerator2 = a*b*np.sqrt((b*np.cos(phi))**2 + (a*a-p*p)*(np.sin(phi))**2)
    demoninator = (a*np.sin(phi))**2 + (b*np.cos(phi))**2
    result = ((numerator1 + numerator2)/demoninator) - p*sec(phi)
    return(-result)

def r_5(p, phi, a, b):
    """Calculates the distance between the two interceps of the ellipse on the 
    line of integration."""
    
    numerator = 2*a*b*np.sqrt((b*np.cos(phi))**2 + (a*a-p*p)*(np.sin(phi))**2)
    denominator = a**2 * (np.sin(phi))**2 + b**2 * (np.cos(phi))**2
    result = numerator/denominator
    return(result)

def theta_1(h, phi, a, b):
    """Upper limit of the theta integral for the case p = 0."""
    result = np.arctan(r_1(phi, a, b)/h)
    return(result)

def theta_2(p, h, phi, a, b):
    """Upper limit of the theta integral for the case p < a."""
    result = np.arctan(r_2(p, phi, a, b)/h)
    return(result)

def theta_3(p, h, phi, a, b):
    """Upper limit of the theta integral for the case p = a."""
    result = np.arctan(r_3(p, phi, a, b)/h)
    return(result)

def theta_4(p, h, phi, a, b):
    """Lower limit of the theta integral for the case p > a."""
    result = np.arctan(r_4(p, phi, a, b)/h)
    return(result)

def theta_5(p, h, phi, a, b):
    """Upper limit of the theta integral for the case p > a."""
    result = np.arctan((r_4(p, phi, a, b) + r_5(p, phi, a, b))/h)
    return(result)

def integrand_1(phi, h, a, b):
    """Solid angle integrand for the case p = 0."""
    result = np.cos(theta_1(h, phi, a, b))
    return(result)

def integrand_2(phi, p, h, a, b):
    """Solid angle integrand for the case p < a."""
    result = np.cos(theta_2(p, h, phi, a, b))
    return(result)

def integrand_3(phi, p, h, a, b):
    """Solid angle integrand for the case p = a."""
    result = 2*np.cos(theta_3(p, h, phi, a, b))
    return(result)

def integrand_4(phi, p, h, a, b):
    """Solid angle integrand for the case p > a."""
    result = 2*(np.cos(theta_4(p,h,phi,a,b)) - np.cos(theta_5(p,h,phi,a,b)))
    return(result)

def phi_max_calc(p,a,b):
    """Upper limit of the phi integral for the case p > a."""
    y = np.sqrt(b**2/(p**2 - a**2))
    result = np.arctan(y)
    return(result)

def solid_angle_calc(p, h, a, b):
    """Calculates the solid angle subtended by an ellipse at a point on axis of
    the semi-axis of the ellipse. This function is not vectorized as it uses 
    scipy.integrate.quad.
    
    Inputs:
     p - The distance between the point of interest and the centre of the 
         ellipse, projected into the plane of the ellipse. Should not be 
         negative.
     h - The perpendicular distance between the point of interest and the plane
         of the ellipse. Should not be negative.
     a - The semi-axis of the ellipse along which the point of interest lies 
         when it is projected into the plane of the ellipse.
     b - The other semi-axis of the ellipse.
    
    Output:
     omega - The solid angle subtended by the ellipse from the point of 
             interest
    """
    
    if p == 0:
        omega = np.pi - integrate.quad(integrand_1, 0, np.pi, args=(h,a,b,))[0]
        omega = 2*omega
    elif p < a:
        omega = np.pi - integrate.quad(integrand_2, 0, np.pi, args=(p,h,a,b,))[0]
        omega = 2*omega
    elif p == a:
        omega = np.pi - integrate.quad(integrand_3, 0, np.pi/2, args=(p,h,a,b,))[0]
    else:
        phi_max = phi_max_calc(p, a, b)
        omega = integrate.quad(integrand_4, 0, phi_max, args=(p,h,a,b,))[0]
        
    return(omega)

def solid_angle_calc2(p, h, a, b):
    """Calculates the solid angle subtended by an ellipse at a point on axis of
    the semi-axis of the ellipse. Uses a double integration, recommended to use
    solid_angle_calc instead. This function is not vectorized as it uses 
    scipy.integrate.quad.
    
    Inputs:
     p - The distance between the point of interest and the centre of the 
         ellipse, projected into the plane of the ellipse. Should not be 
         negative.
     h - The perpendicular distance between the point of interest and the plane
         of the ellipse. Should not be negative.
     a - The semi-axis of the ellipse along which the point of interest lies 
         when it is projected into the plane of the ellipse.
     b - The other semi-axis of the ellipse.
    
    Output:
     omega - The solid angle subtended by the ellipse from the point of 
             interest
    """
    
    if p == 0:
        th_1 = lambda y : theta_1(h, y, a, b)
        integrand = lambda y,z : np.sin(y)
        omega = 2*integrate.dblquad(integrand, 0, np.pi, lambda x : 0, th_1)[0]
    elif p < a:
        th_2 = lambda y : theta_2(p, h, y, a, b)
        integrand = lambda y,z : np.sin(y)
        omega = 2*integrate.dblquad(integrand, 0, np.pi, lambda x : 0, th_2)[0]
    elif p == a:
        th_3 = lambda y: theta_3(p, h, y, a, b)
        integrand = lambda y,z : 2*np.sin(y)
        omega = integrate.dblquad(integrand, 0, np.pi/2, lambda x : 0, th_3)[0]
    else:
        phi_max = phi_max_calc(p, a, b)
        th_4 = lambda y : theta_4(p, h, y, a, b)
        th_5 = lambda y : theta_5(p, h, y, a, b)
        integrand = lambda y,z : 2*np.sin(y)
        omega = integrate.dblquad(integrand, 0, phi_max, th_4, th_5)[0]
    
    return(omega)
