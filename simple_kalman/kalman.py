#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 08:39:14 2022

@author: eadali
"""
from warnings import warn
from numpy import isscalar, array, dot, transpose, zeros
from numpy.linalg import inv
from scipy.linalg import solve_continuous_are as care



def asmatrix(x):
    if isscalar(x):
        x = [[x]]
    return array(x)    
    
    
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
        self.set_model_gains(A, B, C, D, Q, R)
        self.set_initial_value(None, None)

    def __call__(self, t, u, y):
        return self.integrate(t, u, y)

    def set_model_gains(self, A, B, C, D, Q, R):
        # t0, x0 = self.get_initial_value()
        # self.set_initial_value(t0, x0)
        self.A, self.B = asmatrix(A), asmatrix(B)
        self.C, self.D = asmatrix(C), asmatrix(D)
        self.Q, self.R = asmatrix(Q), asmatrix(R)

    def get_model_gains(self):
        return self.A, self.B, self.C, self.D, self.Q, self.R

    def set_initial_value(self, t0, x0):
        self.t0, self.x0 = t0, x0

    def get_initial_value(self):
        return self.t0, self.x0

    def __set_none_value(self, t, x):
        """Set None states for first cycle."""
        t0, x0 = self.get_initial_value()
        if t0 is None:
            t0 = t
        if x0 is None:
            x0 = x
        self.set_initial_value(t0, x0)

    def __check_monotonic_timestamp(self, t0, t):
        """Check timestamp is monotonic."""
        if t < t0:
            msg = 'Current timestamp is smaller than initial timestamp.'
            warn(msg, RuntimeWarning)
            return False
        return True

    def integrate(self, t, u, y):
        self.__set_none_value(t, zeros(self.A.shape[0]))
        t0, x0 = self.get_initial_value()
        # Check monotonic timestamp
        if not self.__check_monotonic_timestamp(t0, t):
            t0 = t
        # Calculate time step
        dt = t - t0
        # ---------------------
        P = care(self.A.T, self.C.T, self.Q, self.R)
        K = dot(dot(P, transpose(self.C)), inv(self.R))
        print(x0)
        input(dot(self.A , x0) )
        dxdt = (dot(self.A, x0) + dot(self.B, u)
                + K * (y - dot(self.C, self.x0) - self.D * u))
        x = x0 + dt * dxdt
        self.set_initial_value(t, x)
        return x
        
    
    
#TODO: not use dot


# a = -1
# b = 1
# c = 1
# d = 0
# q = 0.01
# r = 0.01
# from control import lqe
# print(lqe(a,1,c,q,r))
# kal = Kalman(a, b, c, d, q, r)
# kal.integrate(0, 0, 0)




