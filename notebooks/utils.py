'''A collection of miscellaneous utility functions
'''

from typing import Union, Tuple

import numpy as np
import matplotlib.pyplot as plt
from constants import *


def P_to_a(period: Union[float, np.ndarray], m1: Union[float, np.ndarray],
           m2: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    '''Binary separation from a known period

    Parameters
    ----------
    period : `float/array`
       Binary period in days

    m1 : `float/array`
       Mass of primary star in Msun

    m2 : `flota/array`
       Mass of secondary star in Msun

    Returns
    -------
    a : `float/array`
       Binary separation in Rsun
    '''

    period = period * 24e0 * 3600e0  # in sec
    m1 = m1 * Msun; m2 = m2 * Msun  # in g

    separation = np.power(standard_cgrav * (m1+m2) * np.square(period/(2*pi)), one_third)

    return separation / Rsun


def a_to_P(separation: Union[float, np.ndarray], m1: Union[float, np.ndarray],
           m2:Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    '''Orbital period from a known separation

    Parameters
    ----------
    a : `float/array`
       Binary separation in Rsun

    m1: `float/array`
       Mass of primary star in Msun

    m2: `float/array`
       Mass of secondary star in Msun

    Returns
    -------
    P : `float/array`
       Binary period in days
    '''

    separation = separation * Rsun  # in cm
    m1 = m1 * Msun; m2 = m2 * Msun   # in g

    period = np.power(separation*separation*separation / (standard_cgrav * (m1+m2)), 0.5e0)
    period = (2*pi) * period

    return period / (24e0 * 3600e0)


def a_to_f(separation: Union[float, np.ndarray], m1: Union[float, np.ndarray],
           m2: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    '''Converts semi-major axis to orbital frequency

    Parameters
    ----------
    separation : `float/array`
       Semi-major axis

    m1 : `float/array`
       Primary mass

    m2 : `float/array`
       Secondary mass

    Returns
    -------
    f_orb : `float/array`
       Orbital frequency
    '''

    separation = separation * Rsun  # in cm
    m1 = m1 * Msun; m2 = m2 * Msun  # in g

    f_orb = np.power(standard_cgrav * (m1 + m2) / separation**3, 0.5) / (2*pi)

    return f_orb


def binary_orbits_after_kick(a: float, m1:float, m2:float, m1_remnant_mass: float,
                   w: Union[float, np.ndarray], theta: Union[float, np.ndarray],
                   phi:Union[float, np.ndarray], ids:Union[float, np.ndarray],
                   verbose: bool=False) -> Tuple[Union[float, np.ndarray],
                   Union[float, np.ndarray], Union[float, np.ndarray], Union[float, np.ndarray],
                   Union[float, np.ndarray], Union[float, np.ndarray]]:
    '''Function to compute binary orbital parameters after an asymmetric core-collapse

    Assuming an initial circular orbit, this function calculates the binary configuration after a
    SN explosion with an asymmetric random component. Based on the work of Kalogera (1996)

    Parameters
    ----------
    a : `float`
       Pre-SN separation in Rsun.

    m1 : `float`
       Mass of collapsing star pre-SN in Msun.

    m2 : `float`
       Mass of companion in Msun.

    m1_remnant_mass : `float`
       Gravitational mass of compact object in Msun.

    w : `float/array`
       Natal kick velocity in km/s.

    theta : `float/array`
       Polar angle of kick.

    phi : `float/array`
       Azimutal angle of kick velocity.

    verbose : `bool`
       Flag to control additional output to user.

    Returns
    -------
    a_post : `float/array`
       Post-SN separation in Rsun.

    P_post : `float/array`
       Post-SN orbital period in days.

    e : `float/array`
       Eccentricity of binary post-SN.

    cos_i : `float/array`
       Cosine of the inclination between pre & post SN orbits.

    v_sys : `float/array`
       Systemic velocity post-SN in km/s
    '''

    if verbose: print('calculating post core-collapse orbits for {} kicks'.format(len(w)))

    # Input values conversion to cgs
    a = a * Rsun
    m1 = m1 * Msun
    m2 = m2 * Msun
    m1_remnant_mass = m1_remnant_mass * Msun
    w = w * 1e5

    # Velocity pre-SN
    v_pre = np.sqrt(standard_cgrav * (m1 + m2) / a)

    # Kick velocity (w) must be projected to (x,y,z)
    wx = w * np.cos(phi) * np.sin(theta)
    wy = w * np.cos(theta)
    wz = w * np.sin(phi) * np.sin(theta)

    # Eqs. (3), (4) & (5) of Kalogera (1996)
    a_post = standard_cgrav * (m1_remnant_mass + m2) / \
             (2 * standard_cgrav * (m1_remnant_mass + m2)/a - w**2 - v_pre**2 - 2 * wy * v_pre)
    e = np.sqrt(1 - (wz**2 + wy**2 + v_pre**2 + 2 * wy * v_pre) * a**2 /
             (standard_cgrav * (m1_remnant_mass + m2) * a_post))


    # only interested in bounded binaries
    bounded_mask = (a_post > 0) & (e < 1)
    a_post = a_post[bounded_mask]
    e = e[bounded_mask]
    wx = wx[bounded_mask]
    wy = wy[bounded_mask]
    wz = wz[bounded_mask]

    ids_post = ids[bounded_mask]

    if verbose:
        print('\t{} binaries remain bounded ({:5.2f} percent)'.format(len(e), len(e)/len(w)*100))
        print('\t{} binaries become unbounded ({:5.2f} percent)'.format(len(w)-len(e), (len(w)-len(e))/len(w)*100))

    # update natal kick distro after verbose due to use of len(w)
    try:
        w = w[bounded_mask]
    except TypeError:
        w = w
    try:
        theta = theta[bounded_mask]
    except:
        theta = theta
    try:
        phi = phi[bounded_mask]
    except:
        phi = phi

    # Inclination between pre & post SN orbits. Eq. (11) in Kalogera, 1996
    cos_i = (wy + v_pre) / np.sqrt(wz**2 + (wy + v_pre)**2)

    # Systemic velocity post-SN = eq.(34) x eq.(34), of Kalogera, 1996
    v_sys_2 = np.power((m1_remnant_mass*wx), 2) + np.power((m1_remnant_mass*wz), 2) + \
             np.power((m1_remnant_mass*wy - (m1-m1_remnant_mass)*m2/(m1+m2)*v_pre), 2)
    v_sys_2 = v_sys_2 / np.power((m1_remnant_mass+m2), 2)
    v_sys = np.sqrt(v_sys_2) / 1.e5

    # get orbital period of bounded binaries
    P_post = a_to_P(a_post/Rsun, m1_remnant_mass/Msun, m2/Msun)

    return a_post/Rsun, P_post, e, cos_i, v_sys, w/1e5, theta, phi, ids_post
