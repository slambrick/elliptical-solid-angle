# elliptical-solid-angle
A python module to calculate the solid angle subtended by an ellipse from a point.

## Synopsis

Analytic solid angle integrals for an ellipse exsist and have been derived by Abbas et al. 2015, doi:10.1016/j.nima.2014.10.061, and by J.T. Conway 2010, doi:10.1016/j.nima.2009.11.075.

This module includes functions to perform those integrals. Three functions are provided, one which directly uses the forms of the integrals from Abbas et al. another that has had the theta integral performed analytically first, and a third using the form derived by J.Y. Conway. It is recommended that the third of these functions is used as it only requires a single integation and is quicker to evaluate than the integral by Abbas et al., however all three are provided for completness.

## Requirments

The module use `numpy` and `scipy.integrate`. The example scripts `cosine_added.py` and `solid_map.py` also makes use of `matplotlib.pyplot` and `time`. The code was written for python3.

## Code example

The three functions are `solid_angle_calc`, `solid_angle_calc2`, `solid_angle_calc_general`, all of which take four arguments.

`solid_angle_calc(p, q, h, a, b)`

returns the solid angle subtended by an ellipse with semi-axes of `a` and `b` (the point of interest lying along the axis `a` when projected into the plane of the ellipse). `p` is the distance between the point and the centre of the ellipse along the axis of `a`; `q` is the distane between the point and the centre of the ellipse along the axis of `b`. `h` is the perpendicular distance between the poin and the plane of the ellipse.

The script `cosine_added.py` uses the two solid angle functions to calculate the solid angle for a fixed `h` for a number of `p` values. It also modulates the results by a cosine distribution (centred normally to the plane of the ellipse), [link](https://en.wikipedia.org/wiki/Lambert%27s_cosine_law). For both the raw solid anlge and the modulated cosine functions the script makes use of pyplot to plot the results and save them to a `.eps` file.

The script `solid_map.py` maps out the solid angle subtended as a function of distance from the ellipse for all three cases and times them.

## Licence

This code is provided with the MIT licence.
