# elliptical-solid-angle
A python module to calculate the solid angle subtended by an ellipse from a point.

## Synopsis

Analytic solid angle integrals for an ellipse exsist and have been derived by Abbas et al. 2015, doi:10.1016/j.nima.2014.10.061.

This module includes functions to perform those integrals. Two functions are provided, one which directly uses the forms of the integrals from Abbas et al. and another that has had the theta integral performed analytically first. It is recommended that the second of these functions is used as it only requires a single integation and not a double, however both are provided for completness.

## Code example

The two functions are solid_angle_calc and solid_angle_calc2, both of which take four arguments.
    solid_angle_calc(p, h, a, b)
returns the solid angle subtended by an ellipse with semi-axes of a and b (the point of interest lying along the axis a when projected into the plane of the ellipse). p is the distance between the point and the centre of the ellipse when the point is projected into the plane of the ellipse. h is the perpendicular distance between the poin and the plane of the ellipse.

## Licence

This code is provided with the MIT licence.