#!/usr/bin/env python
# The script is for calculating dynamic topography
# coding: utf-8
# Post_DyTopo_pro.py 的简单版本，只显示重要信息 - 2023-12-19
# 

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
from matplotlib import colors
from os import system
from scipy.io import loadmat
from scipy.interpolate import griddata

plt.style.use('../my.mplstyle')


def subfig_pcolor(ax, Title, Label):
#     ax.set_title(Title)
    ax.grid(True, linestyle='dotted', linewidth=1)
    ax.invert_yaxis()
    # font_tick_params(ax)
    Cbar = plt.colorbar(subfig,orientation='horizontal')
#     ax.axis('scaled')
    # ax.set_xlim(0,1000)
    # ax.set_ylim(400,0)
    Cbar.set_label(Label)  # ,fontweight='bold'
    ax.set_ylabel('Depth ($km$)' )
    ax.set_xlabel('Distance ($km$)'  )
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

#%% Setting
K = 1000.0
Unit = 'km'

Min_SYY, Max_SYY = -10.0e6, 10.0e6
Min_DyTopo,  Max_DyTopo= -500,500
# norm = colors.Normalize()
norm_META = colors.Normalize(vmin=18, vmax=22)
norm_MI = colors.Normalize(vmin=-10, vmax=100)
norm_MRHO = colors.Normalize(vmin=3250, vmax=3600)
#%% Dynamic topography

# filename = 'Ext_NE20_R2km/N_Ext_R2km_10.mat'
# savename = 'Ext_NE20_R2km'
filename = 'No_slab_Sub_mantle_5e22/N_Ext_0.mat'
savename = 'No_slab_Sub_mantle_5e22'



print('filename is ', filename)
data = loadmat(filename)

CX = data['CX']/K
CY = data['CY']/K
NX = data['NX']/K
NY = data['NY']/K
MX = data['MX']/K
MY = data['MY']/K

etan1 = data['etan1']

ystp = data['ystp']
g = data['g']

vy1 = data['vy1']

MI = data['MI']

MVX = data['MVX']
MVY = data['MVY']
META = data['META']
MRHO = data['MRHO']
pr1 = data['pr1']
# EYY
EYY=(vy1[1:,:]-vy1[:-1,:])/ystp

# Topography
Thickness_Air = 20.0 # 20.0
Thickness_l = 0.0 # 200.0

measured_topo = data['measured_topo']
measured_topo[:,0] = measured_topo[:,0] + Thickness_l
measured_topo[:,1] = measured_topo[:,1] + Thickness_Air*1000.0

# 抛弃两侧的点
SYY = 2*etan1*EYY[:,1:-1] 
# SYY = 2*etan1*EYY[:,1:-1] + data['pr1'] 
SYY = 2*etan1*EYY[:,1:-1] - 20000*1000*g * 0 - 1000*2400*g/2 *0

# Dynamic topography
Rho_Mantle = 3300.0
Rho_Crust = 2700.0 # Crust
# Rho_Crust = 2400.0 # Sediment
print('Rho_Crust is ',Rho_Crust)

dRho_Continent = Rho_Crust - 1.0
dRho_Ocean = Rho_Crust - 1000.0
Topo_C = -SYY/dRho_Continent/g
Topo_O = -SYY/dRho_Ocean/g
Topo_Mix = np.zeros_like(Topo_C)
Rho_Mix = np.zeros_like(Topo_C)

measured_topo_interp=griddata(measured_topo[:,0], measured_topo[:,1], CX[0,:])
# measured_topo_interp=griddata(measured_topo[:,0], measured_topo[:,1], CX[0,:], method='nearest')

for index in range(len(measured_topo_interp)):
    # print(measured_topo_interp[index],Topo_C[10,index])
    if measured_topo_interp[index] > 0:
        Topo_Mix[:,index] = -SYY[:,index]/dRho_Continent/g
        Rho_Mix[:,index] = dRho_Continent
    else:
        Topo_Mix[:,index] = -SYY[:,index]/dRho_Ocean/g
        Rho_Mix[:,index] = dRho_Ocean



# 可视化
grid = plt.GridSpec(nrows=7, ncols=2, wspace=0.2, hspace=0.2)
plt.figure(figsize = (12, 15/4*7))  

# 第一个子图的具体排列位置 - Topography
ax = plt.subplot(grid[0, :])
Label = ' '
ax.plot(measured_topo[:,0], measured_topo[:,1]/K*1000, '-')#,label=Label)
plt.axhline(y=0, color='black', linestyle='-', linewidth=1, label='y=400')


### 
file_in = '../New09_T2_BestModel/Best_Model'
topo_LitMod = np.loadtxt(file_in + '/topo_out.dat')

ax.plot(topo_LitMod[:,0], topo_LitMod[:,1]/K*1000, '-')#,label=Label)

ax.set_ylabel('Topography ($m$)')
# ax.legend(ncol=1)
ax.grid(True)
ax.set_xlim(NX[0,0],NX[0,-1])
# ax.set_ylim(-1.0, 1.0)

iM = 10 #4 # 10
print('Thickness_Air is ', Thickness_Air)
print('CX is ', CY[iM,0:4],)
# 第一个子图的具体排列位置 - Dynamic topograhy
ax = plt.subplot(grid[2, :])

ax.plot(CX[iM,:], Topo_Mix[iM,:], '-',label='DyTopo', color='#1f77b4ff')

# ax.plot(CX[iM,:], Topo_C[iM,:], 'y-',label=Label + ' Continent')
# ax.plot(CX[2,:], Topo_C[2,:], '-',label=Label)
# ax.plot(CX[iM,:], Topo_O[iM,:], 'r-',label=Label + ' Ocean')
# font_tick_params(ax)
# ax.set_xlabel('Distance ($km$)')
ax.set_ylabel('Dynamic topo ($m$)')
ax.legend(ncol=1)
ax.grid(True)
ax.set_xlim(NX[0,0],NX[0,-1])
ax.set_ylim(Min_DyTopo,Max_DyTopo)



# plt.show()
# plt.savefig(filename+'_s',bbox_inches='tight', pad_inches=0.1)
plt.savefig(savename+'_s',bbox_inches='tight', pad_inches=0.1)


ans = CX[iM,:], Topo_Mix[iM,:], Rho_Mix[iM,:]

# 为了保存到文本文件，我们需要将ans转换为一个合适的二维数组
# 我们可以横向堆叠这两行
data_to_save = np.column_stack(ans)

# 保存到文本文件，格式为.txt
np.savetxt(filename+"_DynaTopo.txt", data_to_save, fmt='%.2f')  # fmt指定保存时的数值格式


