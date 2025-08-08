#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 12:12:57 2021

@author: ictja
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from decimal import Decimal, ROUND_HALF_UP

def point_in_polygon(point, polygon):
    x, y = point
    n = len(polygon)
    inside = False

    p1x, p1y = polygon[0]
    for i in range(n + 1):
        p2x, p2y = polygon[i % n]
        if y > min(p1y, p2y) and y <= max(p1y, p2y) and x <= max(p1x, p2x):
            if p1y != p2y:
                xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
            if p1x == p2x or x <= xinters:
                inside = not inside
        p1x, p1y = p2x, p2y

    # print(inside)
    return inside



KW93_Dry = {
    'Power_law': 1,
    'd' : 1e-3,  # Grain size, 1 mm
    'b' : 5e-10,  # Burgers vector, m, 0.5 nm
    'mu' : 8e10, # Shear Modulus, 8e10 Pa, 80 GPa
    # Dislocation creep
    'A_dis': 3.5e22,  # preexponential factor, unit, 1/s
    'n_dis' : 3.5,  # stress exponent
    'm_dis' : 0,  # grain size exponent
    'Ea_dis' : 540e3,  # activation energy, unit, J/mol   
    'Va_dis' : 20.0e-6,  # activation volume, unit, m3/mol   ### 15 - 25, 这里取平均值20
    # Difusion creep
    'A_dif' : 8.7e15,
    'n_dif' : 1.00,
    'm_dif' : 2.5,
    'Ea_dif' : 300e3,
    'Va_dif' : 6.0e-6}

KW93_Wet = {
    'Power_law': 1,
    'd' : 1e-3,  # Grain size, 1 mm
    'b' : 5e-10,  # Burgers vector, m, 0.5 nm
    'mu' : 8e10, # Shear Modulus, 8e10 Pa, 80 GPa
    # Dislocation creep / wet
    'A_dis' : 2.0e18,  # s-1
    'n_dis' : 3.0,  # stress exponent
    'm_dis' : 0,  # stress exponent
    'Ea_dis' : 430e3,  # J/mol
    'Va_dis' : 15e-6,  # activation volume, unit, m3/mol    ### 10 - 20, 这里取平均值15
    # Difusion creep / wet
    'A_dif' : 5.3e15,
    'n_dif' : 1.0,
    'm_dif' : 2.5,
    'Ea_dif' : 240e3,
    'Va_dif' : 5.0e-6}


NE20_Wet = {
    'Power_law': 99,
    'd' : 1e-2,  # Grain size, 1 mm
    # Dislocation creep / wet
    'A_dis' : 2.28e-18,  # Pa -3.5 s-1
    'n_dis' : 3.5,  # stress exponent
    'm_dis' : 0,  # stress exponent
    'Ea_dis' : 480e3,  # J/mol
    'Va_dis' : 11e-6,  # activation volume, unit, m3/mol 
    # Difusion creep / wet
    'A_dif' : 4.7e-16,
    'n_dif' : 1,
    'm_dif' : 3,
    'Ea_dif' : 335e3,
    'Va_dif' : 4.0e-6}



def GetETA(T ,P ,EII ,method='KW93_Dry'):
    etamin=1e+18  # Lower limit
    etamax=1e+25   # Upper limit
    
    T = T + 273.15 # unit, ºC to K
    etamin=1e+18   # Lower limit
    etamax=1e+23   # Upper limit
    
    R = 8.314  # Gas constant, J/mol/K    
    
    if method=='KW93_Dry':
        para = KW93_Dry
    if method=='KW93_Wet':
        para = KW93_Wet
    if method=='NE20_Wet':
        para = NE20_Wet
        
    # Power-law: EPSILONii=A*(b/d)**m*(SIGMAii/mu)**n*exp[-(Ea+Va*P)/RT)
    if para['Power_law'] == 1:
        d = para['d']  # Grain size, 1 mm
        b = para['b']  # Burgers vector, m, 0.5 nm
        mu = para['mu']  # Shear Modulus, 8e10 Pa, 80 GPa  
        # Dislocation creep  
        A = para['A_dis']  # material constant, unit, 1/(Pa**n * s * m**m)
        n = para['n_dis']  # stress exponent
        m = para['m_dis']  # grain size exponent
        Ea = para['Ea_dis']  # activation energy, unit, J/mol
        Va = para['Va_dis']  # activation volume, unit, m3/mol    15 - 25
        eterm=(Ea+Va*P)/(R*T)
        if eterm > 100:
            eterm = 100  # eterm should be <= 100
        # eterm[eterm > 100] = 100  # eterm should be <= 100
        eta_dislocation = 0.5 * A**(-1/n) * mu * (d/b) ** (m/n) * EII**((1-n)/n) * np.exp(eterm/n)
        
        
        # Difusion creep / dry
        A = para['A_dif']
        n = para['n_dif']
        m = para['m_dif']
        Ea = para['Ea_dif']
        Va = para['Va_dif']
        if eterm > 100:
            eterm = 100  # eterm should be <= 100
        # eterm[eterm > 100] = 100  # eterm should be <= 100
        eta_difusion = 0.5 * A**(-1/n) * mu * (d/b) ** (m/n) * EII**((1-n)/n) * np.exp(eterm/n)
        
        eta = 1 / (1 / eta_dislocation + 1 / eta_difusion)
    
    
    if para['Power_law'] == 99:
        d = para['d']  # Grain size, 1 mm
        # Dislocation creep  
        A = para['A_dis']  # material constant, unit, 1/(Pa**n * s * m**m)
        n = para['n_dis']  # stress exponent
        m = para['m_dis']  # grain size exponent
        Ea = para['Ea_dis']  # activation energy, unit, J/mol
        Va = para['Va_dis']  # activation volume, unit, m3/mol    15 - 25
        eta_dislocation = 0.5 * A**(-1/n) * d**(m/n) * EII**(1/n-1) * np.exp((Ea+Va*P)/(n*R*T))
        # Difusion creep / dry
        A = para['A_dif']
        n = para['n_dif']
        m = para['m_dif']
        Ea = para['Ea_dif']
        Va = para['Va_dif']
        eta_difusion = 0.5 * A**(-1/n) * d**(m/n) * EII**(1/n-1) * np.exp((Ea+Va*P)/(n*R*T))
        eta = 1 / (1 / eta_dislocation + 1 / eta_difusion)
    
    
    if eta<etamin:
        eta=etamin
    
    if eta>etamax:
        eta=etamax
    
    return eta


# old scrip 老代码
def GetETA_2022(T ,P ,EII ,method='Karato and Wu. 1993 Science'):
    etamin=1e+18  # Lower limit
    etamax=1e+25   # Upper limit
    
    d = 1e-2  # Grain size, m
    R = 8.314  # gas constant
    # # Dislocation creep / dry
    # A_dis = 3.5e22  # material constant, unit, 1/(Pa**n * s * m**m)
    # n_dis = 3.5  # stress exponent
    # m_dis = 0  # grain size exponent
    # Ea_dis = 540000  # activation energy, unit, J/mol
    # Va_dis = 20.0e-6  # activation volume, unit, m3/mol    15 - 25
    # # Difusion creep / dry
    # A_dif = 8.7e15
    # n_dif = 1.00
    # m_dif = -2.5
    # Ea_dif = 300000
    # Va_dif = 6.0

    # Dislocation creep / wet
    A_dis = 2.28e-18  # Pa m3 s-1
    n_dis = 3.5  # stress exponent
    m_dis = 0  # stress exponent
    Ea_dis = 480000  # J/mol
    Va_dis = 1.1e-5  # activation volume, unit, m3/mol    15 - 25
    # Difusion creep / wet
    A_dif = 4.7e-16
    n_dif = 1
    m_dif = 3
    Ea_dif = 335000
    Va_dif = 4.0e-6

    # Power-law: EPSILONii=AD*SIGMAii**n*exp[-(Ea+Va*P)/RT)
	# Iterate for viscosity
	# First viscosity value
   
    eta_dis = 0.5 * (A_dis ** (-1 / n_dis)) * (d ** (m_dis / n_dis)) * (EII ** (1 / n_dis - 1)) * np.exp((Ea_dis + Va_dis * P) / (n_dis * R * (T + 273.15)))
    eta_dif = 0.5 * (A_dif ** (-1 / n_dif)) * (d ** (m_dif / n_dif)) * (EII ** (1 / n_dif - 1)) * np.exp((Ea_dif + Va_dif * P) / (n_dif * R * (T + 273.15)))

    eta = 1 / (1 / eta_dis + 1 / eta_dif)

    if eta<etamin:
        eta=etamin
    
    if eta>etamax:
        eta=etamax
    
    return eta


