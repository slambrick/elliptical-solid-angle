#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 27 15:51:10 2018

Copyright (c) 2019-22, Sam Lambrick
All rights reserverd.
Subject to the MIT licence.

Using the solid angle module to map the solid angles and compare the three
different appraoches.

@author: sam
"""

import numpy as np
import matplotlib.pyplot as plt
import elliptical_solid_angle as esa
from matplotlib.patches import Ellipse
import time


def solid_map_plot(omegas, prange, qrange, title_str, output_file):
    # Create a map of the solid angles as a function of p and q
    fig1, ax1 = plt.subplots(1, 1)
    plt.imshow(omegas, cmap='jet', interpolation='nearest',
               extent=(pmin, pmax, qmin, qmax))
    plt.xlabel('x/mm')
    plt.ylabel('y/mm')
    plt.title(title_str)
    plt.colorbar()

    # Add an ellipse showing to the plot
    ellipse = Ellipse(xy=(0, 0), width=2*a, height=2*b, edgecolor='black',
                      fc='None', lw=2)
    ax1.add_patch(ellipse)

    # Save the plot to an eps file
    plt.savefig(output_file)


if __name__ == '__main__':
    # Semi-axes of the detector aperture
    a = 1
    b = np.sqrt(2)/2

    # Sample pinhole plate distance
    h = 2.121

    # The range of p and q
    pmin = 0.
    pmax = 2
    qmin = 0
    qmax = 1.5

    # The number of values to calculate in each direction
    n_p = 41
    n_q = 31

    # Create numpy arrays of the ps and qs
    ps = np.linspace(pmin, pmax, n_p)
    qs = np.linspace(qmin, qmax, n_q)

    # Create a variable to put the solid angle results in
    omegas1 = np.zeros([n_q, n_p])
    omegas2 = np.zeros([n_q, n_p])
    omegas3 = np.zeros([n_q, n_p])

    # Loop through all the ps and qs
    start = time.time()
    for i in range(n_p):
        for j in range(n_q):
            omegas1[j, i] = esa.solid_angle_calc_general(ps[i], qs[j], h, a, b)
    end = time.time()
    print('Vector potential approach (s):')
    print(end - start)

    start = time.time()
    for i in range(n_p):
        for j in range(n_q):
            omegas2[j, i] = esa.solid_angle_calc2(ps[i], qs[j], h, a, b)
    end = time.time()
    print('Double integral geometric approach (s):')
    print(end - start)

    start = time.time()
    for i in range(n_p):
        for j in range(n_q):
            omegas3[j, i] = esa.solid_angle_calc(ps[i], qs[j], h, a, b)
    end = time.time()
    print('Single integral geometric approach (s):')
    print(end - start)

    # Flip the data for plotting
    omegas1 = np.flipud(omegas1)
    omegas2 = np.flipud(omegas2)
    omegas3 = np.flipud(omegas3)

    solid_map_plot(omegas1, [pmin, pmax], [qmin, qmax], 'Vector potential',
                   'map_solidAngle1.eps')
    solid_map_plot(omegas2, [pmin, pmax], [qmin, qmax],
                   'Double geometric integral', 'map_solidAngle2.eps')
    solid_map_plot(omegas3, [pmin, pmax], [qmin, qmax],
                   'Single geometric integral', 'map_solidAngle3.eps')
