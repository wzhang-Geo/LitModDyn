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

plt.style.use('StokesPY/my.mplstyle')


def font_tick_params(ax):
#     ax.tick_params(axis='y', labelsize=14)
    ax.tick_params(direction='in',length=6,labelsize=16)
    LW=2
    ax.spines['bottom'].set_linewidth(LW) ###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(LW) ####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(LW) ###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(LW) ####设置上部坐标轴的粗细
    
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
    
    # 获取colorbar的位置和尺寸
    pos = Cbar.ax.get_position()
    print(pos)
    # 假设dpi为默认值，通常为100，这意味着每英寸100个点，2厘米约为0.787英寸
    move_distance = 0.025  # 将厘米转换为英寸然后转换为图形坐标系的比例
    # 修改colorbar的位置
    Cbar.ax.set_position([pos.x0, pos.y0 + move_distance, pos.width, pos.height])
    
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

def subfig_quiver_Qkey_Constant(X, Y, VX, VY, K1, K2, Qkey_cmy):    
#     K = 10
    xnum, ynum = X.shape
    x1, y1 = np.arange(0,xnum,K1), np.arange(0,ynum,K2)
    xx, yy = np.meshgrid(x1, y1)
    v_x = X[xx, yy]
    vx = VX[xx, yy]
    v_y = Y[xx, yy]
    vy = VY[xx, yy]
    # Qkey = float(format(np.max((VX,VY)), '.2g'))
    Q = plt.quiver(v_x, v_y, vx, vy, units='xy', color='black',)
    # plt.quiverkey(Q, 0.0, -0.15, Qkey, str(Qkey) + r' $m/s$', labelpos='E', color='black', coordinates='axes')
    # Qkey_cmy = Qkey * 100 * 365.25 * 24 * 3600
    Qkey = Qkey_cmy / 100 / 365.25 / 24 / 3600
    plt.quiverkey(Q, 0.0, -0.15, Qkey, str(Qkey_cmy) + r' $cm/y$', labelpos='E', color='black', coordinates='axes')

#%% Setting
K = 1000.0
Unit = 'km'

Min_SYY, Max_SYY = -10.0e6, 10.0e6
Min_DyTopo,  Max_DyTopo= -500,500
# norm = colors.Normalize()
norm_META = colors.Normalize(vmin=18, vmax=22)
norm_MI = colors.Normalize(vmin=-10, vmax=100)
norm_MRHO = colors.Normalize(vmin=3250, vmax=3600)



#%%
# file_in = '/home/ictja/ownCloud/PhD_Fig/NorProfile_Fig/2021-12-13_V_F02_5/New09_T2_BestModel/Best_Model'
file_in = 'D:\wzhang\StokesPY\Zhang_NorthernPro_2024-05-05\Best_Model'
# 可视化
grid = plt.GridSpec(nrows=7*3, ncols=2, wspace=0.2, hspace=0.2)
# plt.figure(figsize = (12, 15/4*7))  
# plt.figure(figsize = (12, 20))  
plt.figure(figsize = (12, 25*1.5))  
################### 第一个子图的具体排列位置 - Topography
ax = plt.subplot(grid[0:2, :])
Label_measured = r"Observed Topography, $h_{obs}$"
label_iso_with = 'Local isostasy with mantle anomalies'
label_iso_without = 'Local isostasy without mantle anomalies'
label_Regional_iso = r"Regional isostasy, $h_{isost}$"


label_dyTopo_with=r"Total, $h_{total}$"
label_dyTopo_without=r"Purely lithospheric model, $h_{lit}$"

label_dyTopo_Delta=r"Solely mantle anomalies, $h_{anom} = h_{total} - h_{lit}$"
label_dyTopo_res = r"Residual topography, $h_{res} = h_{obs} - h_{isost}$"

measured_topo = np.loadtxt(file_in + '/0Topo.dat')
x = measured_topo[:,0]
y = measured_topo[:,1]
yerr = measured_topo[:,2]

# 计算误差上下界
y_upper = y + yerr
y_lower = y - yerr
# 定义颜色
color = 'gray'  # 可以使用任何matplotlib支持的颜色名称或代码

# 使用阴影表示误差，颜色与主线一致
# plt.plot(x, y, color=color)  # 使用统一的颜色
plt.fill_between(x, y_lower, y_upper, color=color, alpha=0.5, label=Label_measured)
# plt.axhline(y=0, color='black', linestyle='-', linewidth=1)


# LitMod topography + Te
topo_LitMod = np.loadtxt(file_in + '/topo_out.dat')
# ax.plot(topo_LitMod[:,0], topo_LitMod[:,1], '-', color='blue', label=label_iso_with)
# ax.plot(topo_LitMod[:,0], topo_LitMod[:,2], '-', color='black', label=label_iso_without)
topo_Te = np.loadtxt(file_in + '/Deflexiotopotao_TeB.dat')
ax.plot(topo_Te[:,0], topo_Te[:,2], 'g-',label=label_Regional_iso)

ax.set_xlim(np.min(measured_topo[:,0]), np.max(measured_topo[:,0]))
ax.set_ylim(-1000.0, 1600.0)
ax.set_ylabel('Elevation ($m$)')
ax.legend(ncol=1, fontsize="16", loc ="lower right", frameon=False)

# ax.grid(True)

ax.text(0.5, 1.00, ' Northern profile ', transform=plt.gca().transAxes, fontsize=24, va='center', ha='center', bbox=dict(facecolor='white', edgecolor='white'))
ax.text(0.01, 0.95, 'a) Topography $vs$ Regional isostasy', transform=plt.gca().transAxes, fontsize=18, va='top')

ax.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=1)

ax.tick_params(
    axis='x',
    bottom=True,      # 显示下方刻度
    top=True,        # 上方刻度（可选）
    labelbottom=False, # 显示下方数字
    labeltop=False    # 不显示上方数字
)


#--------------------------------------------------------------------------------
################### 第二个子图的具体排列位置 - Residual topograhy
# ax = plt.subplot(grid[1, :])


# ax.plot(topo_LitMod[:,0], measured_topo[:,1]-topo_LitMod[:,1], '-', color='blue', label=label_iso_with)
# plt.fill_between(topo_LitMod[:,0], measured_topo[:,1]-topo_LitMod[:,1] - yerr, measured_topo[:,1]-topo_LitMod[:,1] + yerr, color='blue', alpha=0.5, label=label_iso_with)

# ax.plot(topo_LitMod[:,0], measured_topo[:,1]-topo_LitMod[:,2], '-', color='black', label=label_iso_without)
# plt.fill_between(topo_LitMod[:,0], measured_topo[:,1]-topo_LitMod[:,2] - yerr, measured_topo[:,1]-topo_LitMod[:,2] + yerr, color='black', alpha=0.5, label=label_iso_without)

topo_Te_interp_nearest = griddata(topo_Te[:,0], topo_Te[:,2], measured_topo[:,0], method='nearest')
y_upper = measured_topo[:,1]-topo_Te_interp_nearest + yerr
y_lower = measured_topo[:,1]-topo_Te_interp_nearest - yerr
# ax.plot(measured_topo[:,0], measured_topo[:,1]-topo_Te_interp_nearest, 'g--',label=label_Regional_iso)
# plt.fill_between(measured_topo[:,0], y_lower, y_upper, color='green', alpha=0.5, label=label_Regional_iso)



# ax.set_xlim(np.min(measured_topo[:,0]), np.max(measured_topo[:,0]))
# # ax.set_ylim(-250.0, 250.0)
# ax.set_ylim(-1500.0, 2000.0)
# ax.set_ylabel('Residual Elevation ($m$)')
# ax.legend(ncol=1, fontsize="10", loc ="upper right")
# ax.grid(True)
# ax.text(0.01, 0.95, 'b) Residual topography', transform=plt.gca().transAxes, fontsize=18, fontweight='bold', va='top')


#--------------------------------------------------------------------------------
################### 第三个子图的具体排列位置 - Dynamic topograhy
ax = plt.subplot(grid[2:5, :])

Dytopo = np.loadtxt('D:\wzhang\StokesPY\Zhang_NorthernPro_2024-05-05/' + 'Ext_NE20_R2km/N_Ext_R2km_10.mat_DynaTopo.txt') 
Dytopo_1 = np.loadtxt('D:\wzhang\StokesPY\Zhang_NorthernPro_2024-05-05/' + 'Ext_NE20_R2km/N_Ext_R2km_10.mat_DynaTopo.txt') 

ax.plot(Dytopo_1[:,0] - 200, Dytopo_1[:,1], '-', color='blue', linewidth=1.2, label=label_dyTopo_with)


Dytopo_2 = np.loadtxt('D:\wzhang\StokesPY\Zhang_NorthernPro_2024-05-05/' + 'Ext_NE20_R2km_NoSlab/N_Ext_4.mat_DynaTopo.txt')
ax.plot(Dytopo_2[:,0] - 200, Dytopo_2[:,1], '-', color='black', linewidth=1.2, label=label_dyTopo_without)


ax.plot(Dytopo[:,0] - 200, Dytopo_1[:,1] - Dytopo_2[:,1], '-', color='red', linewidth=1.5, label=label_dyTopo_Delta)


# Dytopo = np.loadtxt('Ext_NE20_R2km_Crust=10E23PaS/N_Ext_4.mat_DynaTopo.txt') 
# ax.plot(Dytopo[:,0] - 200, Dytopo[:,1], '-', color='red', label='$η_{crust}=10e23\ Pa\ s$')

# ax.plot(measured_topo[:,0], measured_topo[:,1]-topo_Te_interp_nearest, 'g--',label='Regional isotasy, Te=20 km')

# Dytopo = np.loadtxt('Ext_KW93_R2km/N_Ext_R2km_4.mat_DynaTopo.txt') 
# ax.plot(Dytopo[:,0] - 200, Dytopo[:,1], '-', color='yellow', label='$η_{crust}=10e23\ Pa\ s$')


# font_tick_params(ax)
# ax.set_xlabel('Distance ($km$)')
ax.set_ylabel('Elevation ($m$)')

# ax.grid(True)
ax.set_xlim(np.min(measured_topo[:,0]), np.max(measured_topo[:,0]))
ax.set_ylim(-1000.0, 1000.0)
ax.fill_between(measured_topo[:,0], y_lower, y_upper, color='green', alpha=0.5, label=label_dyTopo_res)

ax.legend(ncol=1, fontsize="14", loc ="lower right", frameon=False)
# ax.legend(ncol=1,
#           fontsize=16,
#           loc='center left',
#           bbox_to_anchor=(1.02, 0.5),
#           frameon=False)

ax.text(0.01, 0.95, 'b) Dynamic $vs$ Residual topography', transform=plt.gca().transAxes, fontsize=18, va='top')

ax.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=1)
ax.tick_params(
    axis='x',
    bottom=True,      # 显示下方刻度
    top=True,        # 上方刻度（可选）
    labelbottom=False, # 显示下方数字
    labeltop=False    # 不显示上方数字
)


plt.savefig('Fig_Northern_iso_dyn_res_topo',bbox_inches='tight', pad_inches=0.1, dpi=600)
plt.savefig('Fig_Northern_iso_dyn_res_topo.pdf',bbox_inches='tight', pad_inches=0.1, dpi=1200)


# import sys
# sys.exit()
#--------------------------------------------------------------------------------

grid = plt.GridSpec(nrows=7*3, ncols=2, wspace=0.2, hspace=0.2)
# plt.figure(figsize = (12, 15/4*7))  
# plt.figure(figsize = (12, 20))  
plt.figure(figsize = (12, 25*1.5)) 

################### 第四个子图的具体排列位置 - Res Vs Dynamic topograhy
label_Regional_iso = r"Regional isostasy, $h_{isost}$"
ax = plt.subplot(grid[0:2, :])
file_in = 'StokesPY/Zhang_SouthernPro_2024-05-05/2023-05-31_15_24_47_Model_F_13_b11'
measured_topo = np.loadtxt(file_in + '/0Topo.dat')
x = measured_topo[:,0]
y = measured_topo[:,1]
yerr = measured_topo[:,2]

# 计算误差上下界
y_upper = y + yerr
y_lower = y - yerr
color = 'gray'  # 可以使用任何matplotlib支持的颜色名称或代码

# 使用阴影表示误差，颜色与主线一致
# plt.plot(x, y, color=color)  # 使用统一的颜色
plt.fill_between(x, y_lower, y_upper, color=color, alpha=0.5, label=Label_measured)


# LitMod topography + Te
topo_LitMod = np.loadtxt(file_in + '/topo_out.dat')
# ax.plot(topo_LitMod[:,0], topo_LitMod[:,1], '-', color='blue', label=label_iso_with)
# ax.plot(topo_LitMod[:,0], topo_LitMod[:,2], '-', color='black', label=label_iso_without)
topo_Te = np.loadtxt(file_in + '/flexure_tao/Deflexiotopotao_Te30km.dat')
ax.plot(topo_Te[:,0], topo_Te[:,2], 'g-',label=label_Regional_iso)

ax.set_xlim(np.min(measured_topo[:,0]), np.max(measured_topo[:,0]))
ax.set_ylim(-4000.0, 2500.0)
ax.set_ylabel('Elevation ($m$)')
ax.legend(ncol=1, fontsize="16", loc ="lower right", frameon=False,  bbox_to_anchor=(0.85, 0.0))
# ax.grid(True)
ax.tick_params(
    axis='x',
    bottom=True,      # 显示下方刻度
    top=True,        # 上方刻度（可选）
    labelbottom=False, # 显示下方数字
    labeltop=False    # 不显示上方数字
)

pos = ax.get_position()

# 缩短宽度（比如变成原来的 80%）
new_width = pos.width * 1255/1070

# 可以同时控制左边距（居中或偏移）
ax.set_position([pos.x0, pos.y0, new_width, pos.height])

ax.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=1)
ax.text(0.5, 1.00, ' Southern profile ', transform=plt.gca().transAxes, fontsize=24, va='center', ha='center', bbox=dict(facecolor='white', edgecolor='white'))
ax.text(0.01, 0.95, 'a) Topography $vs$ Regional isostasy', transform=plt.gca().transAxes, fontsize=18, va='top')

############################################
ax = plt.subplot(grid[2:5, :])

topo_Te_interp_nearest = griddata(topo_Te[:,0], topo_Te[:,2], measured_topo[:,0], method='nearest')
# ax.plot(measured_topo[:,0], measured_topo[:,1]-topo_Te_interp_nearest, 'g--',label=label_Regional_iso)

# Regional isotasy, Te
y_upper = measured_topo[:,1]-topo_Te_interp_nearest + yerr
y_lower = measured_topo[:,1]-topo_Te_interp_nearest - yerr


# Dytopo = np.loadtxt('Ext_NE20_R2km/N_Ext_R2km_10.mat_DynaTopo.txt') 
Dytopo = np.loadtxt('StokesPY/Zhang_SouthernPro_2024-05-05/'+'Ext_NE20_R2km/S_ext_4.mat_DynaTopo.txt') 
ax.plot(Dytopo[:,0] - 200, Dytopo[:,1], '-', color='blue', linewidth=1.2, label=label_dyTopo_with)


Dytopo = np.loadtxt('StokesPY/Zhang_SouthernPro_2024-05-05/'+'Ext_NE20_R2km_NoSlab/S_ext_4.mat_DynaTopo.txt') 
ax.plot(Dytopo[:,0] - 200, Dytopo[:,1], '-', color='black', linewidth=1.2, label=label_dyTopo_without)

#Dytopo = np.loadtxt('Ext_NE20_R2km_AttachSlab/S_ext_4.mat_DynaTopo.txt') 
#ax.plot(Dytopo[:,0] - 200, Dytopo[:,1], '-', color='orange', label=label_dyTopo_BothAttach)


Dytopo_1 = np.loadtxt('StokesPY/Zhang_SouthernPro_2024-05-05/'+'Ext_NE20_R2km/S_ext_4.mat_DynaTopo.txt') 
Dytopo_2 = np.loadtxt('StokesPY/Zhang_SouthernPro_2024-05-05/'+'Ext_NE20_R2km_NoSlab/S_ext_4.mat_DynaTopo.txt') 
ax.plot(Dytopo[:,0] - 200, Dytopo_1[:,1] - Dytopo_2[:,1], '-', color='red', label=label_dyTopo_Delta)


# Dytopo = np.loadtxt('Ext_NE20_R2km_Crust=10E23PaS/N_Ext_4.mat_DynaTopo.txt') 
# ax.plot(Dytopo[:,0] - 200, Dytopo[:,1], '-', color='red', label='$η_{crust}=10e23\ Pa\ s$')

# ax.plot(measured_topo[:,0], measured_topo[:,1]-topo_Te_interp_nearest, 'g--',label='Regional isotasy, Te=20 km')

# Dytopo = np.loadtxt('Ext_KW93_R2km/N_Ext_R2km_4.mat_DynaTopo.txt') 
# ax.plot(Dytopo[:,0] - 200, Dytopo[:,1], '-', color='yellow', label='$η_{crust}=10e23\ Pa\ s$')


# font_tick_params(ax)
# ax.set_xlabel('Distance ($km$)')
# ax.set_ylabel('Dynamic topo ($m$)')
# ax.legend(ncol=1, fontsize="10", loc ="lower right")

#--------------------------------------------------------------------------------
################### 第四个子图的具体排列位置 - Res Vs Dynamic topography


Dytopo_1 = np.loadtxt('StokesPY/Zhang_SouthernPro_2024-05-05/'+'Ext_NE20_R2km/S_ext_4.mat_DynaTopo.txt') 
Dytopo_2 = np.loadtxt('StokesPY/Zhang_SouthernPro_2024-05-05/'+'Ext_NE20_R2km_NoSlab/S_ext_4.mat_DynaTopo.txt') 
# ax.plot(Dytopo[:,0] - 200, (Dytopo_1[:,1] - Dytopo_2[:,1])*1, '-', color='red', label=label_dyTopo_Delta)

# Dytopo = np.loadtxt('Ext_NE20_R2km_Crust=10E23PaS/N_Ext_4.mat_DynaTopo.txt') 
# ax.plot(Dytopo[:,0] - 200, Dytopo[:,1], '-', color='red', label='$η_{crust}=10e23\ Pa\ s$')

# ax.plot(measured_topo[:,0], measured_topo[:,1]-topo_Te_interp_nearest, 'g--',label=label_Regional_iso)
plt.fill_between(measured_topo[:,0], y_lower, y_upper, color='green', alpha=0.5, label=label_dyTopo_res)

ax.set_ylabel('Elevation ($m$)')


ax.set_xlim(np.min(measured_topo[:,0]), np.max(measured_topo[:,0]))
ax.set_ylim(-800.0, 1200.0)

ax.legend(ncol=1, fontsize="14", loc ="upper right", frameon=False,  bbox_to_anchor=(0.90, 1.025))
# ax.legend(ncol=1,
#           fontsize=16,
#           loc='center left',
#           bbox_to_anchor=(1.02, 0.5),
#           frameon=False)

pos = ax.get_position()

# 缩短宽度（比如变成原来的 80%）
new_width = pos.width * 1255/1070

# 可以同时控制左边距（居中或偏移）
ax.set_position([pos.x0, pos.y0, new_width, pos.height])

ax.text(0.01, 0.95, 'b) Dynamic $vs$ Residual topography', transform=plt.gca().transAxes, fontsize=18, va='top')


ax.set_xlabel('Distance ($km$)')

ax.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=1)
ax.tick_params(
    axis='x',
    bottom=True,      # 显示下方刻度
    top=True,        # 上方刻度（可选）
    labelbottom=False, # 显示下方数字
    labeltop=False    # 不显示上方数字
)



#--------------------------------------------------------------------------------
################### End
# plt.show()
# plt.savefig(filename+'_s',bbox_inches='tight', pad_inches=0.1)
plt.savefig('Fig_Southern_iso_dyn_res_topo',bbox_inches='tight', pad_inches=0.1, dpi=600)
plt.savefig('Fig_Southern_iso_dyn_res_topo.pdf',bbox_inches='tight', pad_inches=0.1, dpi=1200)


# ans = CX[iM,:], Topo_Mix[iM,:]

# # 为了保存到文本文件，我们需要将ans转换为一个合适的二维数组
# # 我们可以横向堆叠这两行
# data_to_save = np.column_stack(ans)

# 保存到文本文件，格式为.txt
# np.savetxt(filename+"_DynaTopo.txt", data_to_save, fmt='%.2f')  # fmt指定保存时的数值格式


