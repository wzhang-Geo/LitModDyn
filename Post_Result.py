#!/usr/bin/env python
# The script is for calculating dynamic topography
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
from os import system
from scipy.io import loadmat

plt.style.use('my.mplstyle')


def subfig_pcolor(ax, Title, Label):
#     ax.set_title(Title)
    ax.grid(True, linestyle='dotted', linewidth=1)
    ax.invert_yaxis()
    # font_tick_params(ax)
    Cbar = plt.colorbar(subfig)
#     ax.axis('scaled')
    # ax.set_xlim(0,1000)
    # ax.set_ylim(400,0)
    Cbar.set_label(Label)  # ,fontweight='bold'
    ax.set_ylabel('Depth ($km$)' )
    # ax.set_xlabel('Distance ($km$)'  )
    # subfig_fileBodies(ax)
    # subfig_quiver(MX[0:30,:], MY[0:30,:], MVX[0:30,:], MVY[0:30,:], 5, 10)
    # subfig_quiver(MX, MY, MVX, MVY, 30, 30)
    
def subfig_quiver(X, Y, VX, VY, K1, K2):    
#     K = 10
    xnum, ynum = X.shape
    x1, y1 = np.arange(0,xnum,K1), np.arange(0,ynum,K2)
    xx, yy = np.meshgrid(x1, y1)
    v_x = X[xx, yy]
    vx = VX[xx, yy]
    v_y = Y[xx, yy]
    vy = VY[xx, yy]
    Qkey = float(format(np.max((VX,VY)), '.2g'))
    Q = plt.quiver(v_x, v_y, vx, vy, units='xy', color='black',)
    plt.quiverkey(Q, 0.0, -0.15, Qkey, str(Qkey) + r' $m/s$', labelpos='E', color='black', coordinates='axes')
#     Qkey_cmy = Qkey * 100 * 365.25 * 24 * 3600
#     plt.quiverkey(Q, 0.0, -0.1, Qkey, str(Qkey_cmy) + r' $cm/y$', labelpos='E', color='red', coordinates='axes')

def subfig_vertical_section(ax, Xlabel):
    ax.legend(ncol=1)
#     ax.set_title(Title, fontdict=font())
    ax.set_xlabel(Xlabel)
    
    ax.invert_yaxis()
    ax.grid(True) # ax.grid(True, linestyle='--', linewidth=2)
#     ax.grid(True, linestyle='--', linewidth=2)
    ax.yaxis.set_ticks_position('both')  

#%% 
K = 1000.0
Unit = 'km'

#%% Load results
# filename = 'NorTransect/NorT_c'

# filename = 'NorTransect/NorT_ShorterSlab'

# filename = 'Zhang_Northern_b2'
filename = 'Zhang_etal_2022_NorthPro/N_Ext_9'
print('filename is ', filename)

# filename = 'Zhang_Southern_ext10'


data = loadmat(filename)


CX = data['CX']/K
CY = data['CY']/K
NX = data['NX']/K
NY = data['NY']/K
MX = data['MX']/K
MY = data['MY']/K
MVX = data['MVX']
MVY = data['MVY']


etan1 = data['etan1']

ystp = data['ystp']
g = data['g']

vy1 = data['vy1']

MI = data['MI']
MRHO = data['MRHO']
MTK = data['MTK']
META = data['META']

measured_topo = data['measured_topo']

#%% Figure setting


Px1 = NX.max()*0.3 # 0.51  # vertical_section
iM_Px1 = (abs(MX[0,:]-Px1)).argmin()
line_style_Px1 = 'r-'
# legend_Px1  = 'x=' + str(MX[0, iM_Px1])    + ' ' + length_unit + ' Markers'
legend_Px1  = 'x=' + str(MX[0, iM_Px1])    + ' ' + Unit

Px2 = NX.max()*0.6
iM_Px2 = (abs(MX[0,:]-Px2)).argmin()
line_style_Px2 = 'b-'
legend_Px2  = 'x=' + str(MX[0, iM_Px2])    + ' ' + Unit

Px3 = NX.max()*0.85
iM_Px3 = (abs(MX[0,:]-Px3)).argmin()
line_style_Px3 = 'g-'
legend_Px3  = 'x=' + str(MX[0, iM_Px3])    + ' ' + Unit

grid = plt.GridSpec(nrows=6, ncols=2, wspace=0.2, hspace=0.2)
plt.figure(figsize = (12, 15))  

#%%
var   = (MI,'Type of rocks','Type of rocks')
# 第一个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[0, :])
subfig = ax.pcolormesh(MX, MY, var[0], vmin=-10,vmax=100) 
subfig_pcolor(ax, var[1], var[2])
ax.set_title(var[1])
ax.set_ylim(60,0)
ax.plot(MX[:, iM_Px1], MY[:, iM_Px1], line_style_Px1)
ax.plot(MX[:, iM_Px2], MY[:, iM_Px2], line_style_Px2)
ax.plot(MX[:, iM_Px3], MY[:, iM_Px3], line_style_Px3)
subfig_quiver(MX, MY, MVX, -MVY, 20, 50)

# 第二个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[1:3, :])  
subfig = ax.pcolormesh(MX, MY, var[0], vmin=-10,vmax=100)
subfig_pcolor(ax, var[1], var[2])
ax.plot(MX[:, iM_Px1], MY[:, iM_Px1], line_style_Px1)
ax.plot(MX[:, iM_Px2], MY[:, iM_Px2], line_style_Px2)
ax.plot(MX[:, iM_Px3], MY[:, iM_Px3], line_style_Px3)
ax.axis('scaled')
ax.set_ylim(NY.max(),NY.min())
subfig_quiver(MX, MY, MVX, -MVY, 30, 30)


# 第三个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[3:6, 0])  
ax.plot(var[0][:, iM_Px1], MY[:, iM_Px1], line_style_Px1, label=legend_Px1)
ax.plot(var[0][:, iM_Px2], MY[:, iM_Px2], line_style_Px2, label=legend_Px2)
ax.plot(var[0][:, iM_Px3], MY[:, iM_Px3], line_style_Px3, label=legend_Px3)
subfig_vertical_section(ax, var[1])

ax.set_ylim(60,0)
# 第四个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[3:6, 1]) 
ax.plot(var[0][:, iM_Px1], MY[:, iM_Px1], line_style_Px1, label=legend_Px1)
ax.plot(var[0][:, iM_Px2], MY[:, iM_Px2], line_style_Px2, label=legend_Px2)
ax.plot(var[0][:, iM_Px3], MY[:, iM_Px3], line_style_Px3, label=legend_Px3)
subfig_vertical_section(ax, var[1])
ax.set_xlim(80,100)

plt.savefig(filename + '_MI',bbox_inches='tight', pad_inches=0.1)

#%%
grid = plt.GridSpec(nrows=6, ncols=2, wspace=0.2, hspace=0.2)
plt.figure(figsize = (12, 15))  


var   = (MRHO,'Density ($kg/m^3$)','$kg/m^3$')
# 第一个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[0, :])
subfig = ax.pcolormesh(MX, MY, var[0]) 
subfig_pcolor(ax, var[1], var[2])
ax.set_title(var[1])
ax.set_ylim(60,0)
ax.plot(MX[:, iM_Px1], MY[:, iM_Px1], line_style_Px1)
ax.plot(MX[:, iM_Px2], MY[:, iM_Px2], line_style_Px2)
ax.plot(MX[:, iM_Px3], MY[:, iM_Px3], line_style_Px3)

# 第二个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[1:3, :])  
subfig = ax.pcolormesh(MX, MY, var[0], vmin=3250,vmax=3600)
subfig_pcolor(ax, var[1], var[2])
ax.plot(MX[:, iM_Px1], MY[:, iM_Px1], line_style_Px1)
ax.plot(MX[:, iM_Px2], MY[:, iM_Px2], line_style_Px2)
ax.plot(MX[:, iM_Px3], MY[:, iM_Px3], line_style_Px3)
ax.axis('scaled')
# ax.set_ylim(400,0)

# 第三个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[3:6, 0])  
ax.plot(var[0][:, iM_Px1], MY[:, iM_Px1], line_style_Px1, label=legend_Px1)
ax.plot(var[0][:, iM_Px2], MY[:, iM_Px2], line_style_Px2, label=legend_Px2)
ax.plot(var[0][:, iM_Px3], MY[:, iM_Px3], line_style_Px3, label=legend_Px3)
subfig_vertical_section(ax, var[1])

ax.set_ylim(60,0)
# 第四个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[3:6, 1]) 
ax.plot(var[0][:, iM_Px1], MY[:, iM_Px1], line_style_Px1, label=legend_Px1)
ax.plot(var[0][:, iM_Px2], MY[:, iM_Px2], line_style_Px2, label=legend_Px2)
ax.plot(var[0][:, iM_Px3], MY[:, iM_Px3], line_style_Px3, label=legend_Px3)
subfig_vertical_section(ax, var[1])
ax.set_xlim(3250,3650)
plt.savefig(filename + '_MRHO',bbox_inches='tight', pad_inches=0.1)


#%% 设置变量及单位
var  = (MTK, 'Temperature ($°C$)','$°C$')

# 声明一个GridSpec对象实例，创建的是2行3列的图像布局。
grid = plt.GridSpec(nrows=6, ncols=2, wspace=0.2, hspace=0.2)

# 设置整个图像大小。
plt.figure(figsize = (12, 15))  

# 第一个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[0, :])
subfig = ax.pcolormesh(MX, MY, var[0], vmin=0,vmax=1200) 
subfig_pcolor(ax, var[1], var[2])
ax.set_title(var[1])
ax.set_ylim(60,0)
ax.plot(measured_topo[:,0], -measured_topo[:,1]/1000, 'w--')
ax.plot(MX[:, iM_Px1], MY[:, iM_Px1], line_style_Px1)
ax.plot(MX[:, iM_Px2], MY[:, iM_Px2], line_style_Px2)
ax.plot(MX[:, iM_Px3], MY[:, iM_Px3], line_style_Px3)

# 第二个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[1:3, :])  
subfig = ax.pcolormesh(MX, MY, var[0], vmin=0,vmax=1600)
subfig_pcolor(ax, var[1], var[2])
ax.axis('scaled')
# ax.set_ylim(400,0)
ax.plot(MX[:, iM_Px1], MY[:, iM_Px1], line_style_Px1)
ax.plot(MX[:, iM_Px2], MY[:, iM_Px2], line_style_Px2)

# 第三个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[3:6, 0])  
ax.plot(var[0][:, iM_Px1], MY[:, iM_Px1], line_style_Px1, label=legend_Px1)
ax.plot(var[0][:, iM_Px2], MY[:, iM_Px2], line_style_Px2, label=legend_Px2)
ax.plot(var[0][:, iM_Px3], MY[:, iM_Px3], line_style_Px3, label=legend_Px3)
subfig_vertical_section(ax, var[1])
ax.set_ylim(60,0)
# 第四个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[3:6, 1]) 
ax.plot(var[0][:, iM_Px1], MY[:, iM_Px1], line_style_Px1, label=legend_Px1)
ax.plot(var[0][:, iM_Px2], MY[:, iM_Px2], line_style_Px2, label=legend_Px2)
ax.plot(var[0][:, iM_Px3], MY[:, iM_Px3], line_style_Px3, label=legend_Px3)
subfig_vertical_section(ax, var[1])
ax.set_xlim(0,1600)
plt.savefig(filename + '_MTK',bbox_inches='tight', pad_inches=0.1)


#%% 设置变量及单位
var = (np.log10(META),'$Log_{10}$(Viscosity)','$Pa·s$')

# 声明一个GridSpec对象实例，创建的是2行3列的图像布局。
grid = plt.GridSpec(nrows=6, ncols=2, wspace=0.2, hspace=0.2)

# 设置整个图像大小。
plt.figure(figsize = (12, 15))  

# 第一个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[0, :])
subfig = ax.pcolormesh(MX, MY, var[0]) 
subfig_pcolor(ax, var[1], var[2])
ax.set_title(var[1])
ax.set_ylim(60,0)
# ax.plot(measured_topo[:,0], -measured_topo[:,1]/1000, 'w--')
ax.plot(MX[:, iM_Px1], MY[:, iM_Px1], line_style_Px1)
ax.plot(MX[:, iM_Px2], MY[:, iM_Px2], line_style_Px2)
ax.plot(MX[:, iM_Px3], MY[:, iM_Px3], line_style_Px3)

# 第二个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[1:3, :])  
subfig = ax.pcolormesh(MX, MY, var[0], vmin=19,vmax=20, cmap='viridis') 
subfig_pcolor(ax, var[1], var[2])

mask = MI < 89  
Z_Lit = np.where(mask, var[0], np.nan)
ax.pcolormesh(MX, MY, Z_Lit)


ax.axis('scaled')
# ax.set_ylim(400,0)
ax.plot(MX[:, iM_Px1], MY[:, iM_Px1], line_style_Px1)
ax.plot(MX[:, iM_Px2], MY[:, iM_Px2], line_style_Px2)
ax.plot(MX[:, iM_Px3], MY[:, iM_Px3], line_style_Px3)


# 第三个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[3:6, 0])  
ax.plot(var[0][:, iM_Px1], MY[:, iM_Px1], line_style_Px1, label=legend_Px1)
ax.plot(var[0][:, iM_Px2], MY[:, iM_Px2], line_style_Px2, label=legend_Px2)
ax.plot(var[0][:, iM_Px3], MY[:, iM_Px3], line_style_Px3, label=legend_Px3)
subfig_vertical_section(ax, var[1])
ax.set_ylim(60,0)
# 第四个子图的具体排列位置为(0,0)。
ax = plt.subplot(grid[3:6, 1]) 
ax.plot(var[0][:, iM_Px1], MY[:, iM_Px1], line_style_Px1, label=legend_Px1)
ax.plot(var[0][:, iM_Px2], MY[:, iM_Px2], line_style_Px2, label=legend_Px2)
ax.plot(var[0][:, iM_Px3], MY[:, iM_Px3], line_style_Px3, label=legend_Px3)
subfig_vertical_section(ax, var[1])
# ax.set_xlim(0,1600)
plt.savefig(filename + '_META',bbox_inches='tight', pad_inches=0.1)
