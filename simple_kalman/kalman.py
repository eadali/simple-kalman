#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 08:39:14 2022

@author: eadali
"""
def RK4(fun, t_span, y0, n):
    """Explicit Runge-Kutta method of order 4.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system. The calling signature is fun(t, y).
    t_span : array_like
        Interval of integration (t0, tf).
    y0 : array_like
        Initial state.
    n : int
        Number of integration steps.

    Returns
    -------
    t : float
        Integration end time.
    y : ndarray
        The integrated value at t
    """
    # Integration initial and final time
    t0, tf = t_span
    t, y = t0, asarray(y0)
    # Calculate step-size
    h = (tf - t0) / n
    for i in range(n):
        # Calculate slopes
        k1 = asarray(fun(t,         y))
        k2 = asarray(fun(t+(h/2.0), y + h * (k1/2.0)))
        k3 = asarray(fun(t+(h/2.0), y + h * (k2/2.0)))
        k4 = asarray(fun(t+h,       y + h * k3))
        # Update time and states
        t = t + h
        y = y + (1.0/6.0) * h * (k1 + 2*k2 + 2*k3 + k4)
    return t, y


class Kalman:
    def __init__(self, A, B, C, D, Q, R):
        self.A, self.B, self.C, self.D = A, B, C, D
        self.Q, self.R = Q, R
