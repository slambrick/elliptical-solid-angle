"""
Created on Fri Dec 15 14:13:12 2017

Copyright (c) 2017-22, Sam Lambrick
All rights reserverd.
Subject to the MIT licence.

This script makes use of the elliptical_solid_angle module to calculate solid
angles and uses the functions from that module to add a cosine distribution to
model the proportion of an emitter of gas that passes through an elliptical
appeture.

The script calculates the solid angle and 'intensity' caught by the apeture for
a fixed perpendicular distance h for a range of parrallel distances p then
plots and saves the results using matplotlib.pyplot. The particular values or
p,h,a,b chosen are related to research work that I do in helium atom
microscopy.
"""

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import elliptical_solid_angle as esa

# Use the ggplot style because I like it
plt.style.use("ggplot")



def displace_ellipse():
    """Consideres the chaning cosine intensity and solid angle as the point of
    interest is dispaced from the ellipse."""

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
        omegas[i] = esa.solid_angle_calc(ps[i], 0, hh, aa, bb)
        intensities[i] = esa.cosine_intensity(ps[i], hh, aa, bb)
        omegas2[i] = esa.solid_angle_calc2(ps[i], 0, hh, aa, bb)
        intensities2[i] = esa.cosine_intensity2(ps[i], hh, aa, bb)

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


if __name__ == '__main__':
    displace_ellipse()
