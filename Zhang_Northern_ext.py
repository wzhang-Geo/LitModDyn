#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# """
# Created on Thu Jun 17 18:14:31 2021
# 对代码进行了简化 - 2023-08-24

# @author: wzhang
# """

from __future__ import division # Changing the Division Operator
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
from scipy.interpolate import griddata
import main_lit


# This is northern transect from Zhang et al., (2022) JGR

print("    Model initialization...")
print("    This is examples from Zhang et al., (2022)")
Thickness_Air = 20*1000.0
Thickness_bottom = 200*1000.0
Thickness_l = 200*1000.0
print('    The thickness of sticky air:', Thickness_Air/1000, ' km')
print('    The thickness of bottom layer:', Thickness_bottom/1000, ' km')
print('    The thickness of each side layers:', Thickness_l/1000, ' km')

g=9.8              # Acceleration of Gravity, m/s^2

# file_in = '/home/wzhang/ownCloud/PhD_Fig/NorProfile_Fig/2021-12-13_V_F02_5/New09_T2_BestModel/Best_Model'

file_in = '/home/wzhang/ownCloud/PhD_Fig/NorProfile_Fig/2021-12-13_V_F02_5/New09_T2_BestModel/Case2_No_Slabs'
print(file_in)


measured_topo = np.loadtxt(file_in + '/0Topo.dat')
measured_topo[:,0] = measured_topo[:,0] + Thickness_l/1000 # unit, km
measured_topo[:,1] = measured_topo[:,1] - Thickness_Air

if Thickness_l==0:
    polygon_LAB = np.loadtxt(file_in + '/LAB_py.out')
else:
    polygon_LAB = np.loadtxt(file_in + '/LAB_py_ex.out')
polygon_LAB[:,0] = polygon_LAB[:,0] * 1000 + Thickness_l
polygon_LAB[:,1] = polygon_LAB[:,1] * -1000 + Thickness_Air
plt.figure()
plt.plot(polygon_LAB[:,0],-polygon_LAB[:,1])

data = np.loadtxt(file_in + '/post_processing_output_for_stokes.dat')
data[:,0] = data[:,0] + Thickness_l/1000
data[:,1] = data[:,1] - Thickness_Air/1000
xsize_Lit = (data[-1,0]-data[0,0])*1000 + Thickness_l * 2
ysize_Lit = -data[-1,1]*1000 + Thickness_bottom
xsize = xsize_Lit
ysize = ysize_Lit
#%%


# Defining resolution, 2 or 5 km
xres = 2e3
yres = 2e3
xnum = int(xsize/xres)+1
ynum = int(ysize/yres)+1
print(" The resolution is ", xres/1000, yres/1000, " km")

# Defining gridsteps
xstp=xsize/(xnum-1)
ystp=ysize/(ynum-1)
  
# Defining number of markers and steps between them in the horizontal and vertical direction
xmx=4 #number of markers per cell in horizontal direction
ymy=4 #number of markers per cell in vertical direction
mxnum=(xnum-1)*xmx #total number of markers in horizontal direction
mynum=(ynum-1)*ymy #total number of markers in vertical direction
mxstep=xsize/mxnum #step between markers in horizontal direction   
mystep=ysize/mynum #step between markers in vertical direction

# Creating nodes, center & markers arrays
x0 = np.arange(0 * xstp, xsize + xstp, xstp)
y0 = np.arange(0 * ystp, ysize + ystp, ystp)
NX,NY = np.meshgrid(x0, y0)

x1 = np.arange(0 * xstp, xsize + 2*xstp, xstp)
y1 = np.arange(0 * ystp, ysize + 2*ystp, ystp)
VX_X,VX_Y = np.meshgrid(x0, y1 - 0.5 * ystp)
VY_X,VY_Y = np.meshgrid(x1 - 0.5 * xstp, y0)

x2 = np.arange(0.5 * xstp, xsize, xstp)
y2 = np.arange(0.5 * ystp, ysize, ystp)
CX,CY = np.meshgrid(x2, y2)

x3 = np.arange(0.5 * mxstep, xsize, mxstep)
y3 = np.arange(0.5 * mystep, ysize, mystep)
MX,MY = np.meshgrid(x3, y3)
del x0, y0, x1, y1, x2, y2, x3, y3


# Creating markers arrays
MI=np.zeros((mynum,mxnum)) # Type
MRHO=np.zeros((mynum,mxnum)) # Density
META=np.zeros((mynum,mxnum)) # viscosity
MTK=np.zeros((mynum,mxnum)) # Temperture
MPR_LitMod  = np.zeros((mynum,mxnum))

MEII  = np.ones((mynum,mxnum))*1.0e-17

# 对measured_topo进行插值
# 进行插值，线性内插和最近点外插
measured_topo_interp = griddata(measured_topo[:,0]*1000, measured_topo[:,1], MX[0,:], method='linear')
measured_topo_interp_nearest = griddata(measured_topo[:,0]*1000, measured_topo[:,1], MX[0,:], method='nearest')

# 对于线性插值产生的NaN，使用最近点方法的结果填充
mask = np.isnan(measured_topo_interp)
measured_topo_interp[mask] = measured_topo_interp_nearest[mask]

plt.figure()
plt.plot(measured_topo[:,0],measured_topo[:,1],'r.',  label='Elevation')
plt.plot(MX[0,:]/1000,measured_topo_interp, label='Elevation by interp')
plt.legend(ncol=1)

#%% Method II: 基于多边形识别区域
LitModXY = np.column_stack((data[:, 0], -1 * data[:, 1])) * 1000

Method = 'nearest'
MI = griddata(LitModXY, data[:, 7], (MX, MY), method=Method)
MRHO = griddata(LitModXY, data[:, 6], (MX, MY), method=Method)
MTK = griddata(LitModXY, data[:, 2], (MX, MY), method=Method)
MPR_LitMod = griddata(LitModXY, data[:, 3], (MX, MY), method=Method)


# 注意，这里是初始数值，判断类型，尤其是在细小的间隙中，Lit_mant: MI=60
for xm in range(mxnum):
    for ym in range(mynum):
        if MY[ym,xm]<0 - measured_topo_interp[xm] : # * 0:	     # sticky air
            MI[ym,xm] = -10
            MRHO[ym,xm] = 1000 # 1000
            MTK[ym,xm] = 0            
            META[ym,xm] = 1.0e18  
        elif MI[ym,xm]<11:    # sediment & crust
            META[ym,xm] = 1.0e21
        # elif main_lit.point_in_polygon((MX[ym,xm],MY[ym,xm]), polygon_LAB):# No slabs 没有板片
        elif main_lit.point_in_polygon((MX[ym,xm],MY[ym,xm]), polygon_LAB) or MI[ym,xm]<91:
            META[ym,xm] = 1.0e22  
            
            MI[ym,xm] = 60

        else:
            META[ym,xm] = 1.0e19 # mantle
            # META[ym,xm] = main_lit.GetETA(MTK[ym,xm],MPR_LitMod[ym,xm],MEII[ym,xm])


#%% 检查结果 MI
var   = (MI,'Type of rocks','Type of rocks')
plt.figure(figsize=(6, 8))
plt.subplot(211)
plt.pcolormesh(MX / 1000, MY / 1000, MI, vmin=-10,vmax=100, cmap='Spectral_r') 

plt.grid(True, linestyle='dotted', linewidth=1)
plt.gca().invert_yaxis()
Cbar = plt.colorbar()
plt.ylim(60,0)
Cbar.set_label(var[1])  # ,fontweight='bold'
plt.xlabel('Distance ($km$)' )
plt.ylabel('Depth ($km$)' )

plt.plot(MX[0,:]/1000,-measured_topo_interp / 1000, 'r-',label='Elevation by interp')

plt.subplot(212)
plt.pcolormesh(MX / 1000, MY / 1000, MI, vmin=-10,vmax=100, cmap='Spectral_r') 

plt.grid(True, linestyle='dotted', linewidth=1)
plt.gca().invert_yaxis()
Cbar = plt.colorbar()
# plt.ylim(400,0)
Cbar.set_label(var[1])  # ,fontweight='bold'
plt.xlabel('Distance ($km$)' )
plt.ylabel('Depth ($km$)' )
    
# plt.plot(measured_topo[:,0],measured_topo[:,1],'r.',  label='Elevation')
plt.plot(MX[0,:]/1000,-measured_topo_interp / 1000, 'r-',label='Elevation by interp')

# subfig_pcolor(ax, var[1], var[2])
# ax.set_title(var[1])
# ax.set_ylim(60,0)
# ax.plot(MX[:, iM_Px1], MY[:, iM_Px1], line_style_Px1)
# ax.plot(MX[:, iM_Px2], MY[:, iM_Px2], line_style_Px2)
# ax.plot(MX[:, iM_Px3], MY[:, iM_Px3], line_style_Px3)
# subfig_quiver(MX, MY, MVX, -MVY, 20, 50)


#########################################
# META=np.ones((mynum,mxnum))*1.0e21 # All viscosity

# for xm in range(mxnum):
#     for ym in range(mynum):
#         if MY[ym,xm]<0 - measured_topo_interp[xm]:	     # stick air
#             MI[ym,xm] = -10
#             MRHO[ym,xm] = 1000
#             MTK[ym,xm] = 0            
#             META[ym,xm] = 1.0e18  
#         elif MI[ym,xm]<11:
#             META[ym,xm] = 1.0e21    # sediment & crust
#         elif MI[ym,xm]<99:
#             META[ym,xm] = 1.0e22    # lithosphere mantle
#         else:
#             # META[ym,xm] = 1.0e19 # mantle
#             META[ym,xm] = 1.0e22 # mantle
#             # META[ym,xm] = main_lit.GetETA(MTK[ym,xm],MPR_LitMod[ym,xm],MEII[ym,xm])

# =============================================================================
# #%% Method I: scipy.interpolate.griddata 基于地幔组分识别区域
# # nearest 结果可以插值但是板片与岩石圈地幔有间隔，linear插值出现Nan
# 
# LitModXY = np.column_stack((data[:, 0], -1 * data[:, 1])) * 1000
# 
# Method = 'nearest'
# # Method = 'linear'
# MI = griddata(LitModXY, data[:, 7], (MX, MY), method=Method)
# # MI[0:3,:] = griddata(LitModXY, data[:, 7], (MX[0:3,:], MY[0:3,:]), method='nearest')
# 
# MRHO = griddata(LitModXY, data[:, 6], (MX, MY), method=Method)
# # MRHO[0:3,:] = griddata(LitModXY, data[:, 6], (MX[0:3,:], MY[0:3,:]), method='nearest')
# 
# MTK = griddata(LitModXY, data[:, 2], (MX, MY), method=Method)
# # MTK[0:3,:] = griddata(LitModXY, data[:, 2], (MX[0:3,:], MY[0:3,:]), method='nearest')
# 
# MPR_LitMod = griddata(LitModXY, data[:, 3], (MX, MY), method=Method)
# # MPR_LitMod[0:3,:] = griddata(LitModXY, data[:, 3], (MX[0:3,:], MY[0:3,:]), method='nearest')
# 
# # META=np.ones((mynum,mxnum)) * 1.0e+21
# 
# for xm in range(mxnum):
#     for ym in range(mynum):
#         if MI[ym,xm]<0:     # stick air
#             META[ym,xm] = 1.0e18   # stick air
#         elif MI[ym,xm]<11:
#             META[ym,xm] = 1.0e21    # sediment & crust
#         elif MI[ym,xm]<99:
#             META[ym,xm] = 1.0e22    # lithosphere mantle
#             
#         else:
#             META[ym,xm] = 1.0e19 # mantle
#             # META[ym,xm] = main_lit.GetETA(MTK[ym,xm],MPR_LitMod[ym,xm],MEII[ym,xm])
# =============================================================================
        
