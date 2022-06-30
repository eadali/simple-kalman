#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 08:39:14 2022

@author: eadali
"""

class Kalman:
    def __init__(self, A, B, C, D, Q, R):
        self.A, self.B, self.C, self.D = A, B, C, D
        self.Q, self.R = Q, R