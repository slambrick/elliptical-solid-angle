# Copyright (c) 2022, Sam Lambrick
# All rights reserverd.
# Subject to the MIT licence.
#
# Contains the nessacery functions to perform analytical solid angle
# integrals for ellipses. Also contains a function for integrating and
# calculating the solid angle.
#
# The derivation of the functions was performed by Abbas et al. 2015, for an
# understanding of the equations used it is recommened to look there. The
# variables follow the convention used there.
# doi:10.1016/j.nima.2014.10.061
#
# The derivation of the simpler general function (using vector potentials) was
# done by John T. Conway 2010. His formual should be used in general, however the
# geometric derivation of Abbas et al. can be combined more readily with other
# functions.
# doi:10.1016/j.nima.2009.11.075
#
# The theta integral from Abbas et al. can be performed analytically. This has
# been used here along with the double integral form. The arguments p and q are
# distances and should be positive.
#
# Variables used here, using the same naming convention as used by Abbas:
#     p - the distance between the centre of the ellipse and the point of
#         interest in the plane of the ellipse along the major axis of the
#         ellipse
#     q - the distance between the centre of the ellipse and the point of
#         interest in the plane of the ellipse along the minor axis of the
#         ellipse
#     a - the semi-major axis of the ellipse
#     b - the semi-minor axis of the ellipse
#     h - the vertical perpendicular distance between the point of interest and
#         the plane of the ellpise
#     rho - sqrt(p^2 + q^2), the distance from the centre of the ellipse to the
#           point of interest projected into the plane of the ellipse
#
# Solid angle functions:
#     calc_omega        - performs a single integration to calculate the
#                         solid angle using integral forms from John T. Conway
#     solid_angle_calc  - performs a single integration to calculate the solid
#                         angle using the integral forms from Abbas et al. with
#                         the theta integral performed
#     solid_angle_calc2 - performs a double inetegration, using the form directly
#                         from Abbas et al., to calculate the solid angles
#
# Calling syntax:
#     solid_angle_calc...(p, q, a, b, h)
#
# None of these functions return the error estimate on the integrals, however
# it is easy to modify them to do so.

library(pracma)


# Calculates the sectant, 1/cos(theta).
sec <- function(theta) 1/cos(theta)


## Geometric functions for Abbas ---------------------------------------------#

# Calculates the distance from the point of interest to the edge of the
# ellipse, in the plane of the ellipse. For the case p = 0."""
r_1<- function(phi, a, b) {
    denominator = (a*sin(phi))^2 + (b*cos(phi))^2
    result = a*b/sqrt(denominator)
    return(result)
}

# Calculated the distance from the point of interest to the edge of the
# ellipse, in the plane of the ellipse. For the case p < a."""
r_2 <- function (p, phi, a, b) {
    numerator1 = -p*b*b*cos(phi)
    numerator2 = a*b*sqrt((b*cos(phi))^2 + (a*a - p*p)*(sin(phi))^2)
    denominator = (a*sin(phi))^2 + (b*cos(phi))^2
    result = (numerator1 + numerator2)/denominator
    return(result)
}

# Calculates the distance from the point of interest to the edge of the
# ellipse, in the plane of the ellipse. For the case of p = a."""
r_3 <- function(p, phi, a, b) {
    numerator = 2*a*b*b*cos(phi)
    denominator = (a*sin(phi))^2 + (b*cos(phi))^2
    result = numerator/denominator
    return(result)
}

# This is the negative of what it is in the paper
# Calculates the distance from the point of interest to the edge of the
# ellipse projected in the plane of the ellipse. For the case p > a."""
r_4 <- function(p, phi, a, b) {
    numerator1 = a*a*p*(sin(phi))^2 * sec(phi)
    numerator2 = a*b*sqrt((b*cos(phi))^2 + (a*a-p*p)*(sin(phi))^2)
    demoninator = (a*sin(phi))^2 + (b*cos(phi))^2
    result = ((numerator1 + numerator2)/demoninator) - p*sec(phi)
    return(-result)
}

#Calculates the distance between the two interceps of the ellipse on the
# line of integration."""
r_5 <- function(p, phi, a, b) {
    numerator = 2*a*b*sqrt((b*cos(phi))^2 + (a*a-p*p)*(sin(phi))^2)
    denominator = a^2 * (sin(phi))^2 + b^2 * (cos(phi))^2
    result = numerator/denominator
    return(result)
}

# Calculates the distance between the point of interest to the edge of the
# ellipse for the case rho < a, rho < b."""
r_6 <- function(p, q, phi, a, b) {
    tmp = (a*sin(phi))^2 + (b*cos(phi))^2
    numerator1 = p*b*b*cos(phi) + q*a*a*sin(phi)
    numerator2 = a*b*sqrt(tmp - (p*tan(phi) - q)^2 * (cos(phi))^2)
    denominator = tmp
    result = (-numerator1 + numerator2)/denominator
    return(result)
}

# Calculates the distance between the point of interest and the edge of
# the ellipse (projected into the plane of the ellipse) along the line of
# integration."""
r_8 <- function(p, q, phi, a, b) {
    tmp = (a*sin(phi))^2 + (b*cos(phi))^2
    numerator1 = a*a*sin(phi)*(p*tan(phi) - q)
    numerator2 = a*b*sqrt(tmp - (p*sin(phi) - q*cos(phi))^2)
    denominator = tmp
    result = (numerator1 + numerator2)/denominator - p*sec(phi)
    return(-result)
}

# Calculated the distance from the edge of the ellipse to the other edge
# along the line of integration for the case rho > a, rho > b, for a specific
# value of phi."""
r_9 <- function(p, q, phi, a, b) {
    tmp = (a*sin(phi))^2 + (b*cos(phi))^2
    numerator = 2*a*b*sqrt(tmp - (p*sin(phi) - q*cos(phi))^2)
    denominator = tmp
    result = numerator/denominator
    return(result)
}

## Integration limits for Abbas -----------------------------------------------#

# Upper limit of the theta integral for the case p = 0."""
theta_1 <- function(h, phi, a, b) arctan(r_1(phi, a, b)/h)

# Upper limit of the theta integral for the case p < a."""
theta_2 <- function(p, h, phi, a, b) arctan(r_2(p, phi, a, b)/h)

# Upper limit of the theta integral for the case p = a."""
theta_3 <- function(p, h, phi, a, b) arctan(r_3(p, phi, a, b)/h)

# Lower limit of the theta integral for the case p > a."""
theta_4 <- function(p, h, phi, a, b) arctan(r_4(p, phi, a, b)/h)

# Upper limit of the theta integral for the case p > a."""
theta_5 <- function(p, h, phi, a, b) arctan((r_4(p, phi, a, b) + r_5(p, phi, a, b))/h)

# Upper limit of integration for the case rho < a, rho < b."""
theta_6 <- function(p, q, h, phi, a, b) arctan(r_6(p, q, phi, a, b)/h)

# Lower limit of integration for the case rho > a, rho > b."""
theta_8 <- function(p, q, h, phi, a, b) arctan(r_8(p, q, phi, a, b)/h)

# Upper limit of integration for the case rho > a, rho > b."""
theta_9 <- function(p, q, h, phi, a, b) arctan((r_8(p, q, phi, a, b) + r_9(p, q, phi, a, b))/h)


## Integrands for Abbas ------------------------------------------------------#

# Solid angle integrand for the case p = 0."""
integrand_1 <- function(phi, h, a, b) cos(theta_1(h, phi, a, b))

# Solid angle integrand for the case p < a."""
integrand_2 <- function(phi, p, h, a, b)cos(theta_2(p, h, phi, a, b))

# Solid angle integrand for the case p = a."""
integrand_3 <- function(phi, p, h, a, b) 2*cos(theta_3(p, h, phi, a, b))

# Solid angle integrand for the case p > a."""
integrand_4 <- function(phi, p, h, a, b) 2*(cos(theta_4(p, h, phi, a, b)) - cos(theta_5(p, h, phi, a, b)))

# Solid angle integrand for the case rho < a, rho < b."""
integrand_5 <- function(phi, p, q, h, a, b) cos(theta_6(p, q, h, phi, a, b))

# Solid angle integrand for the case rho > a, rho > b."""
integrand_6 <- function(phi, p, q, h, a, b) cos(theta_8(p, q, h, phi, a, b)) - cos(theta_9(p, q, h, phi, a, b))

# Upper limit of the phi integral for the case p > a."""
phi_max_calc <- function(p, a, b) {
    y = sqrt(b^2/(p^2 - a^2))
    arctan(y)
}

# Lower limit of the phi integral for the case rho > a, rho > b."""
phi_prime_max_calc <- function(p, q, a, b) {
    arg_numerator = p*q - sqrt(-a*a*b*b + b*b*p*p + a*a*q*q)
    arctan(arg_numerator/(p*p - a*a))
}

# Upper limit of the phi integral for the case rho > a, rho > b."""
phi_primeprime_max_calc <- function(p, q, a, b) {
    arg_numerator = p*q + sqrt(-a*a*b*b + b*b*p*p + a*a*q*q)
    arctan(arg_numerator/(p*p - a*a))
}

## Abbas solid angle calculations --------------------------------------------#

#Calculates the solid angle subtended by an ellipse at a point on axis of
#    the semi-axis of the ellipse. This function is not vectorized as it uses
#    scipy.integrate.quad.
#
#    Inputs:
#     p - The distance between the point of interest and the centre of the
#         ellipse, projected into the plane of the ellipse along the major axis.
#         Should not be negative.
#     q - The distance between the point of interest and the centre of the
#         ellipse, projected into the plane of the elllipse along the minor
#         axis. Should not be negative
#     h - The perpendicular distance between the point of interest and the plane
#         of the ellipse. Should not be negative.
#     a - The semi-axis of the ellipse along which the point of interest lies
#         when it is projected into the plane of the ellipse.
#     b - The other semi-axis of the ellipse.
#
#    Output:
#     omega - The solid angle subtended by the ellipse from the point of
#             interest
solid_angle_calc <- function(p, q, h, a, b) {
    if (p == 0 && q != 0) {
        # No need to preform the more complicated integral if the point of
        # interest lies along either axis of the ellipse
        p = q
        q = 0
        tmp = a
        a = b
        b = tmp
    }
    if (q == 0) {
        # The cases of the point of interest lying along an axis of the ellipse
        if (p == 0) {
            f = function(th) integrand_1(th, h, a, b)
            omega = pi - integrate.quad(f, 0)
            omega = 2*omega
        } else if (p < a) {
            f <- function(th) integrand_2(th, p, h, a, b)
            omega = pi - integrate.quad(f, 0, pi)
            omega = 2*omega
        } else if (p == a) {
            f <- function(th) integrand_3(th, p, h, a, b)
            omega = pi - integrate.quad(f, 0, pi/2)
        } else {
            phi_max = phi_max_calc(p, a, b)
            f <- function(th) integrand_4(th, p, h, a, b)
            omega = integrate.quad(f, 0, phi_max)
        }
    } else {
        # The case of the point of interest not lying along an axis of the
        # ellispe
        if ((p^2/a^2 + q^2/b^2) < 1) {
            # Point inside the ellipse
            f <- function(th) integrand_5(th, p, q, h, a, b)
            omega = 2*pi - integrate.quad(g, 0, 2*pi)
        } else if(p > a || q > b) {
            # Point outiside of the bounding box
            if (q > b && p <= a) {
                tmp = q
                q = p
                p = tmp
                tmp = a
                a = b
                b = tmp
            }
            phi_min = phi_prime_max_calc(p, q, a, b)
            phi_max = phi_primeprime_max_calc(p, q, a, b)
            f <- function(th) integrand_6(th, p, q, h, a, b)
            omega = integrate.quad(f, phi_min, phi_max)
        } else {
            # Point outside of the ellipse but inside the bounding box
            # TODO: derive expression for this
            omega = calc_omega(p, q, h, a, b)
        }
    }
    return(omega)
}


def solid_angle_calc2(p, q, h, a, b):
    """Calculates the solid angle subtended by an ellipse at a point on axis of
    the semi-axis of the ellipse. Uses a double integration, recommended to use
    solid_angle_calc instead. This function is not vectorized as it uses
    scipy.integrate.quad.

    Inputs:
     p - The distance between the point of interest and the centre of the
         ellipse, projected into the plane of the ellipse along the major axis
         of the ellipse. Should not be negative.
     q - The distance between the point of interest and the centre of the
         ellipse, projected into the plane of the ellipse, along the minor axis
         of the ellipse. Should not be negative
     h - The perpendicular distance between the point of interest and the plane
         of the ellipse. Should not be negative.
     a - The semi-major-axis of the ellipse.
     b - The semi-minor-axis of the ellipse.

    Output:
     omega - The solid angle subtended by the ellipse from the point of
             interest
    """

if p == 0 and q != 0:
    # No need to preform the more complicated integral if the point of
    # interest lies along either axis of the ellipse
    p = q
q = 0
tmp = a
a = b
b = tmp

if q == 0:
    # The cases of the point of interest lying along an axis of the ellipse
    if p == 0:
    th_1 = lambda y: theta_1(h, y, a, b)
integrand = lambda y,z: sin(y)
omega = 2*integrate.dblquad(integrand, 0, pi, lambda x: 0,
                            th_1)[0]
elif p < a:
    th_2 = lambda y: theta_2(p, h, y, a, b)
integrand = lambda y,z: sin(y)
omega = 2*integrate.dblquad(integrand, 0, pi, lambda x: 0,
                            th_2)[0]
elif p == a:
    th_3 = lambda y: theta_3(p, h, y, a, b)
integrand = lambda y,z: 2*sin(y)
omega = integrate.dblquad(integrand, 0, pi/2, lambda x: 0,
                          th_3)[0]
else:
    phi_max = phi_max_calc(p, a, b)
th_4 = lambda y: theta_4(p, h, y, a, b)
th_5 = lambda y: theta_5(p, h, y, a, b)
integrand = lambda y,z: 2*sin(y)
omega = integrate.dblquad(integrand, 0, phi_max, th_4, th_5)[0]
else:
    # The case of the point of interest not lying along an axis of the
    # ellispe
    if (p^2/a^2 + q^2/b^2) < 1:
    # Point inside the ellipse
    th_6 = lambda y: theta_6(p, q, h, y, a, b)
integrand = lambda y,z: sin(y)
omega = integrate.dblquad(integrand, 0, 2*pi, lambda x: 0,
                          th_6)[0]
elif p > a or q > b:
    # Point outiside of the bounding box
    if q > b and p <= a:
    tmp = q
q = p
p = tmp
tmp = a
a = b
b = tmp
th_8 = lambda y: theta_8(p, q, h, y, a, b)
th_9 = lambda y: theta_9(p, q, h, y, a, b)
phi_min = phi_prime_max_calc(p, q, a, b)
phi_max = phi_primeprime_max_calc(p, q, a, b)
integrand = lambda y,z: sin(y)
omega = integrate.dblquad(integrand, phi_min, phi_max, th_8,
                          th_9)[0]
else:
    # Point outside of the ellipse but inside the bounding box
    # TODO: derive and implement this
    omega = solid_angle_calc_general(p, q, h, a, b)

return(omega)


## Cosine modified integrals -------------------------------------------------#


# Two functions are included, one that performs a double numerical integral
# and the other that has had the theta integral performed analytically. The
# double integral is included for completness, and as a demonstartion of how the
# inclusion of distributions can be achieved.

# Calculates the intensity of a cosine distribution that enters an
# elliptical appeture."""
cosine_intensity <- function(p, h, a, b) {
    if p < a:
        th_2 = lambda y: theta_2(p, h, y, a, b)
    integrand = lambda y: cos(2*th_2(y))
    II = pi/2 - 0.25*integrate.quad(integrand, 0, 2*pi)[0]
    elif p == a:
        th_3 = lambda y: theta_3(p, h, y, a, b)
    integrand = lambda y: cos(2*th_3(y))
    II = pi/4 - 0.5*integrate.quad(integrand, 0, pi/2)[0]
    else:
        phi_max = phi_max_calc(p, a, b)
    th_4 = lambda y: theta_4(p, h, y, a, b)
    th_5 = lambda y: theta_5(p, h, y, a, b)
    integrand = lambda y: cos(2*th_4(y)) - cos(2*th_5(y))
    II = 0.5*integrate.quad(integrand, 0, phi_max)[0]
    
    return(II)
}

def cosine_intensity2(p, h, a, b):
    """Calculates the intensity of a cosine distribution that enters an
    elliptical appeture. Uses a double integral."""

if p == 0:
    th_1 = lambda y : theta_1(h, y, a, b)
integrand = lambda y,z : sin(y)*cos(y)
II = 2*integrate.dblquad(integrand, 0, pi, lambda x: 0, th_1)[0]
elif p < a:
    th_2 = lambda y : theta_2(p, h, y, a, b)
integrand = lambda y,z : sin(y)*cos(y)
II = 2*integrate.dblquad(integrand, 0, pi, lambda x: 0, th_2)[0]
elif p == a:
    th_3 = lambda y: theta_3(p, h, y, a, b)
integrand = lambda y,z : 2*sin(y)*cos(y)
II = integrate.dblquad(integrand, 0, pi/2, lambda x: 0, th_3)[0]
else:
    phi_max = phi_max_calc(p, a, b)
th_4 = lambda y : theta_4(p, h, y, a, b)
th_5 = lambda y : theta_5(p, h, y, a, b)
integrand = lambda y,z : 2*sin(y)*cos(y)
II = integrate.dblquad(integrand, 0, phi_max, th_4, th_5)[0]

return(II)


## Conway method functions ---------------------------------------------------#

# Integrand for the general case from 'John T. Conway' 'Analytic solution
# for the solid angle subtended at a point source radiation vector potential'
# which claims to be general, unlike Abbas et.al. who do not consider the
# case inside the bounding box but outisde the ellipse."""
integrand_general <- function(phi, p, q, h, a, b) {
    tmp = p*p + q*q + h*h + 2*a*p*cos(phi) + 2*b*q*sin(phi) + 
        (a*cos(phi))^2 + (b*sin(phi))^2
    term1 = 1 - h/(sqrt(tmp))
    tmp = tmp - h*h
    term2 = (a*b + p*b*cos(phi) + q*a*sin(phi))/(tmp)
    result = term1*term2
    return(result)
}

# Calculates the solid angle using the formula derived by John T. Conway
# 2010. Is faster than the alternative integral in 'solid_angle_calc'.
#
# Inputs:
#   p - The distance between the point of interest and the centre of the
#       ellipse, projected into the plane of the ellipse along the major axis.
#       Should not be negative.
#   q - The distance between the point of interest and the centre of the
#       ellipse, projected into the plane of the elllipse along the minor
#       axis. Should not be negative
#   h - The perpendicular distance between the point of interest and the plane
#       of the ellipse. Should not be negative.
#   a - The semi-axis of the ellipse along which the point of interest lies
#       when it is projected into the plane of the ellipse.
#   b - The other semi-axis of the ellipse.
#
# Output:
#   omega - The solid angle subtended by the ellipse from the point of
#           interest
calc_omega <- function(p, q, h, a, b) {
    f <- function(th) integrand_general(th, p, q, h, a, b)
    omega = integral(f, 0, 2*pi)
    return(omega)
}
