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

plt.rcParams['font.family'] = 'Arial'
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
    move_distance = 0.15  # 将厘米转换为英寸然后转换为图形坐标系的比例
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
    #Q = plt.quiver(v_x, v_y, vx, vy, units='xy', color='black',)
    Q = plt.quiver(v_x, v_y, vx, vy, units='xy', color='black', scale=11e-12)
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
#%% Dynamic topography
# filename = 'NorTransect/NorT_ex_t20'

# filename = 'N_ext_NoSlabs/Zhang_Northern_b10.mat'
# savename = 'N_ext_NoSlabs'

# filename = 'Ext_NE20_R2km/N_Ext_R2km_10.mat'
# savename = 'Ext_NE20_R2km'

# # filename = 'TB_NE20_R2km/N_Ext_2.mat'
# # savename = 'TB_NE20_R2km'


# print('filename is ', filename)
# data = loadmat(filename)

# CX = data['CX']/K
# CY = data['CY']/K
# NX = data['NX']/K
# NY = data['NY']/K
# MX = data['MX']/K
# MY = data['MY']/K

# etan1 = data['etan1']

# ystp = data['ystp']
# g = data['g']

# vy1 = data['vy1']

# MI = data['MI']

# MVX = data['MVX']
# MVY = data['MVY']
# META = data['META']
# MRHO = data['MRHO']
# pr1 = data['pr1']
# # EYY
# EYY=(vy1[1:,:]-vy1[:-1,:])/ystp

# # Topography
# Thickness_Air = 20.0 # 20.0
# Thickness_l = 200.0 # 200.0


#%%
# file_in = '/home/wzhang/ownCloud/PhD_Fig/NorProfile_Fig/2021-12-13_V_F02_5/New09_T2_BestModel/Best_Model'
file_in = '../2023-05-31_15_24_47_Model_F_13_b11'

# 可视化
grid = plt.GridSpec(nrows=7, ncols=2, wspace=0.2, hspace=0.2)
# plt.figure(figsize = (12, 15/4*7))  
plt.figure(figsize = (12*1, 25*1))  

################### 第一个子图的具体排列位置 - Topography
ax = plt.subplot(grid[0, :])
Label = 'Observed elevation'
label_dyTopo_with=r"Total, $h_{total}$"
label_dyTopo_without=r"Purely lithospheric model, $h_{lit}$"
label_dyTopo_Delta='From mantle anomalies, A-B'

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
### plt.fill_between(x, y_lower, y_upper, color=color, alpha=0.5, label=Label)
### plt.axhline(y=0, color='black', linestyle='-', linewidth=1)


# LitMod topography + Te
topo_LitMod = np.loadtxt(file_in + '/topo_out.dat')
### ax.plot(topo_LitMod[:,0], topo_LitMod[:,1], '-', color='blue', label=label_dyTopo_with)
### ax.plot(topo_LitMod[:,0], topo_LitMod[:,2], '-', color='black', label=label_dyTopo_without)
# topo_Te = np.loadtxt(file_in + '/Deflexiotopotao_TeB.dat')
### ax.plot(topo_Te[:,0], topo_Te[:,2], 'g--',label='Regional isotasy, Te=20 km')

### ax.set_xlim(np.min(measured_topo[:,0]), np.max(measured_topo[:,0]))

### ax.set_ylim(-2000.0, 2000.0)
### ax.set_ylabel('Elevation ($m$)')
### ax.legend(ncol=1, fontsize="8", loc ="lower right")
### ax.grid(True)
### ax.text(0.01, 0.95, 'a)', transform=plt.gca().transAxes, fontsize=20, fontweight='bold', va='top')

###################################################
# ax.legend(ncol=1, fontsize="14", loc ="lower right")

y_upper = y + yerr
y_lower = y - yerr
plt.fill_between(x, 0, y, color=(167/255, 194/255, 223/255), alpha=0.5)
plt.fill_between(x, -4000, y, color=(200/255, 200/255, 200/255), alpha=1)
plt.plot(x, y, color='darkred', linewidth=1,label=Label)

ax.legend(ncol=1, fontsize="16", loc ="lower right") #, frameon=False)

# plt.text(0, 200, 'N Tyrrhenian Basin', ha='center')
plt.text(200, 200, 'S Tyrrhenian Basin', ha='center')
plt.text(500, 1000, 'S Apennines', ha='center')
plt.text(720, 200, 'Adriatic Sea', ha='center')
plt.text(920, 1800, 'S Dinarides', ha='center')
###################################################
ax.set_xlim(np.min(measured_topo[:,0]), np.max(measured_topo[:,0]))
ax.set_ylim(-4000.0, 2500)

ax.set_ylabel('Elevation ($m$)')

ax.tick_params(
    axis='x',
    bottom=True,      # 显示下方刻度
    top=True,        # 上方刻度（可选）
    labelbottom=False, # 显示下方数字
    labeltop=False    # 不显示上方数字
)

# ax.grid(True)
ax.text(0.01, 0.95, 'a)  ', transform=plt.gca().transAxes, fontsize=20, fontweight='bold', va='top')


#--------------------------------------------------------------------------------
################### 第二个子图的具体排列位置 - Dynamic topograhy
ax = plt.subplot(grid[1, :])

Dytopo_total = np.loadtxt('Sub_mantle_NE20_Wet/S_Ext_9.mat_DynaTopo.txt') 
# ax.plot(Dytopo_total[:,0] - 210, Dytopo_total[:,1], '-', color='blue', label=r'$\eta_{sub\_mantle}=HK03$')

Dytopo_total_1e19 = np.loadtxt('Sub_mantle_1e19/S_Ext_0.mat_DynaTopo.txt') 
# ax.plot(Dytopo_total_1e19[:,0] - 210, Dytopo_total_1e19[:,1], '-', color= 'brown', label=r'$\eta_{sub\_mantle}=1e19$')

# Dytopo = np.loadtxt('Sub_mantle_1e20/N_Ext_1.mat_DynaTopo.txt') 
# ax.plot(Dytopo[:,0] - 210, Dytopo[:,1], '-', label=r'$\eta_{sub\_mantle}=1e20$')

Dytopo_total_5e22 = np.loadtxt('Sub_mantle_1e19/S_Ext_0.mat_DynaTopo.txt') 
# ax.plot(Dytopo_total_5e22[:,0] - 210, Dytopo_total_5e22[:,1], '-', color= 'orange', label=r'$\eta_{sub\_mantle}=1e22$')

ax.tick_params(
    axis='x',
    bottom=True,      # 显示下方刻度
    top=True,        # 上方刻度（可选）
    labelbottom=False, # 显示下方数字
    labeltop=False    # 不显示上方数字
)



# Dytopo = np.loadtxt('Ext_NE20_R2km_Crust=10E23PaS/N_Ext_4.mat_DynaTopo.txt') 
# ax.plot(Dytopo[:,0] - 210, Dytopo[:,1], '-', color='red', label='$η_{crust}=10e23\ Pa\ s$')

# # font_tick_params(ax)
# # ax.set_xlabel('Distance ($km$)')
# ax.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=1)

# # ax.set_ylabel('Dynamic topo ($m$)')
# ax.set_ylabel('Elevation ($m$)')

# ax.legend(ncol=1, fontsize="14", loc ="lower right", frameon=False) # ,  bbox_to_anchor=(0.85, 0.0))
# # ax.grid(True)
# # ax.set_xlim(np.min(measured_topo[:,0]), np.max(measured_topo[:,0]))
# ax.set_xlim(-200,1270)
# ax.set_xlim(0,1245)
# ax.set_ylim(-400.0, 400.0)

# ax.text(0.01, 0.95, 'b) Total model', transform=plt.gca().transAxes, fontsize=20, fontweight='bold', va='top')




################### 第三个子图的具体排列位置 - Mantle flow
ax = plt.subplot(grid[1, :])

Dytopo_NoSlab = np.loadtxt('No_slab_Sub_mantle_NE20_Wet/S_Ext_4.mat_DynaTopo.txt') 
# ax.plot(Dytopo_NoSlab[:,0] - 210, Dytopo_NoSlab[:,1], '-', color='blue', label=r'$\eta_{sub\_mantle}=HK03$')

Dytopo_NoSlab_1e19 = np.loadtxt('No_slab_Sub_mantle_1e19/S_Ext_0.mat_DynaTopo.txt') 
# ax.plot(Dytopo_NoSlab_1e19[:,0] - 210, Dytopo_NoSlab_1e19[:,1], '-', color= 'brown', label=r'$\eta_{sub\_mantle}=1e19$')

# Dytopo = np.loadtxt('Sub_mantle_1e20/N_Ext_1.mat_DynaTopo.txt') 
# ax.plot(Dytopo[:,0] - 210, Dytopo[:,1], '-', label=r'$\eta_{sub\_mantle}=1e20$')

# Dytopo_NoSlab_5e22 = np.loadtxt('No_slab_Sub_mantle_5e22/N_Ext_0.mat_DynaTopo.txt') 
Dytopo_NoSlab_5e22 = Dytopo_total_5e22 - (Dytopo_total_1e19 - Dytopo_NoSlab_1e19)*0.2

# ax.plot(Dytopo_NoSlab_5e22[:,0] - 210, Dytopo_NoSlab_5e22[:,1], '-', color= 'orange', label=r'$\eta_{sub\_mantle}=1e22$')

# ax.tick_params(
#     axis='x',
#     bottom=True,      # 显示下方刻度
#     top=True,        # 上方刻度（可选）
#     labelbottom=False, # 显示下方数字
#     labeltop=False    # 不显示上方数字
# )

# ax.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=1)

# # ax.set_ylabel('Dynamic topo ($m$)')
# ax.set_ylabel('Elevation ($m$)')

# ax.legend(ncol=1, fontsize="14", loc ="lower right", frameon=False) # ,  bbox_to_anchor=(0.85, 0.0))
# # ax.grid(True)
# # ax.set_xlim(np.min(measured_topo[:,0]), np.max(measured_topo[:,0]))
# ax.set_xlim(-200,1270)
# ax.set_xlim(0,1245)
# ax.set_ylim(-600.0, 600.0)

# ax.text(0.01, 0.95, 'c) Pure lithospheric model', transform=plt.gca().transAxes, fontsize=20, fontweight='bold', va='top')



################### 第三个子图的具体排列位置 - Mantle flow
ax = plt.subplot(grid[1, :])

Dytopo_Deta_total = Dytopo_total - Dytopo_NoSlab
ax.plot(Dytopo_total[:,0] - 210, Dytopo_Deta_total[:,1], '-',color='blue', label=r'$\eta_{sub\_mantle}=HK03$')

Dytopo_Deta_1e19 = Dytopo_total_1e19 - Dytopo_NoSlab_1e19
ax.plot(Dytopo_total_1e19[:,0] - 210, Dytopo_Deta_1e19[:,1], '-',color= 'brown',   label=r'$\eta_{sub\_mantle}=1e19$')

# Dytopo = np.loadtxt('Sub_mantle_1e20/N_Ext_1.mat_DynaTopo.txt') 
# ax.plot(Dytopo[:,0] - 210, Dytopo[:,1], '-', label=r'$\eta_{sub\_mantle}=1e20$')

Dytopo_Deta_5e22 = Dytopo_total_5e22 - Dytopo_NoSlab_5e22 
# Dytopo = np.loadtxt('Sub_mantle_5e22/N_Ext_1.mat_DynaTopo.txt') 
ax.plot(Dytopo_total_5e22[:,0] - 210, Dytopo_Deta_5e22[:,1], '-', color= 'orange', label=r'$\eta_{sub\_mantle}=1e22$')

ax.tick_params(
    axis='x',
    bottom=True,      # 显示下方刻度
    top=True,        # 上方刻度（可选）
    labelbottom=True, # 显示下方数字
    labeltop=False    # 不显示上方数字
)

ax.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=1)

# ax.set_ylabel('Dynamic topo ($m$)')
ax.set_ylabel('Elevation ($m$)')

ax.legend(ncol=1, fontsize="13", loc ="lower left", frameon=False,  bbox_to_anchor=(0.0, -0.05))
# ax.grid(True)
# ax.set_xlim(np.min(measured_topo[:,0]), np.max(measured_topo[:,0]))
# ax.set_xlim(-200,1270)
ax.set_xlim(0,1245)
ax.set_ylim(-300.0, 300.0)

ax.set_xlabel('Distance ($km$)'  )

ax.text(0.01, 0.95, 'b) ', transform=plt.gca().transAxes, fontsize=20, fontweight='bold', va='top')
# ax.text(0.01, 0.95, 'd) Solely mantle anomalies', transform=plt.gca().transAxes, fontsize=20, fontweight='bold', va='top')

################### End
# plt.show()
# plt.savefig(filename+'_s',bbox_inches='tight', pad_inches=0.1)
plt.savefig('Fig_wzhang_2026-03-24_dRho=2700',bbox_inches='tight', pad_inches=0.1, dpi=600)




import sys
sys.exit()


#--------------------------------------------------------------------------------
################### 第三个子图的具体排列位置 - Mantle flow

###### Model A
ax = plt.subplot(grid[2:4, :])  

filename = 'Sub_mantle_NE20_Wet/N_Ext_9.mat'
print('filename is ', filename)
data = loadmat(filename)
MI = data['MI']
MX = data['MX']/K - 210
MY = data['MY']/K - 20
META = data['META']
MVX = data['MVX']
MVY = data['MVY']

subfig = ax.pcolormesh(MX, MY, np.log10(META), vmin=19,vmax=20, cmap='viridis') 
# subfig_pcolor(ax, 's', 'Log10(Viscosity) ($Pa\ s$)')
mask = MI < 89  
Z_Lit = np.where(mask, np.log10(META), np.nan)
ax.pcolormesh(MX, MY, Z_Lit)

ax.axis('scaled')
ax.set_xlim(0,1245)
ax.set_ylim(400,-20)
ax.set_ylabel('Depth ($km$)' )
# ax.set_xlim(-200,1270)
# ax.set_ylim(600,-20)

# 添加矩形方框
# 参数：(x, y)矩形左下角的坐标, width矩形的宽度, height矩形的高度
import matplotlib.patches as patches
rect = patches.Rectangle((0,-20), 1070, 420, linewidth=2,linestyle='--', edgecolor='black', facecolor='none')
ax.add_patch(rect)

# subfig_quiver(MX, MY, MVX, -MVY, 60, 60)
subfig_quiver_Qkey_Constant(MX, MY, MVX, -MVY, 60, 60, 2)
del MI,MX, MY, META,MVX, MVY

ax.text(0.01, 0.1, 'c)', transform=plt.gca().transAxes, fontsize=20, fontweight='bold', va='top')

pos = ax.get_position()
ax.set_position([pos.x0, pos.y0 + 0.03, pos.width, pos.height])

###### Model B
ax = plt.subplot(grid[4:7, :])  

filename = 'Sub_mantle_1e20/N_Ext_1.mat'
data = loadmat(filename)
MI = data['MI']
MX = data['MX']/K - 210

MY = data['MY']/K - 20
META = data['META']
MVX = data['MVX']
MVY = data['MVY']

subfig = ax.pcolormesh(MX, MY, np.log10(META), vmin=19,vmax=20, cmap='viridis')  # viridis,YlGnBu,YlGnBu_r,Greys
subfig_pcolor(ax, 's', 'Log10(Viscosity) ($Pa\ s$)')
mask = MI < 89  
Z_Lit = np.where(mask, np.log10(META), np.nan)
ax.pcolormesh(MX, MY, Z_Lit)

ax.axis('scaled')
ax.set_xlim(0,1245)
ax.set_ylim(400,-20)

# ax.set_xlim(-200,1270)
# ax.set_ylim(600,-20)
# subfig_quiver(MX, MY, MVX, -MVY, 60, 60)
subfig_quiver_Qkey_Constant(MX, MY, MVX, -MVY, 60, 60, 2)
# subfig_quiver_Qkey_Constant(X, Y, VX, VY, K1, K2, Qkey)   
del MI,MX, MY, META,MVX, MVY


# 添加矩形方框
# 参数：(x, y)矩形左下角的坐标, width矩形的宽度, height矩形的高度
import matplotlib.patches as patches
rect = patches.Rectangle((0,-20), 1070, 420, linewidth=2,linestyle='--', edgecolor='black', facecolor='none')
ax.add_patch(rect)



ax.text(0.01, 0.1, 'd)', transform=plt.gca().transAxes, fontsize=20, fontweight='bold', va='top')
pos = ax.get_position()
ax.set_position([pos.x0, pos.y0 + 0.09, pos.width, pos.height])


#--------------------------------------------------------------------------------
################### End
# plt.show()
# plt.savefig(filename+'_s',bbox_inches='tight', pad_inches=0.1)
plt.savefig('Fig_Sensitivity_test',bbox_inches='tight', pad_inches=0.1, dpi=600)


# ans = CX[iM,:], Topo_Mix[iM,:]

# # 为了保存到文本文件，我们需要将ans转换为一个合适的二维数组
# # 我们可以横向堆叠这两行
# data_to_save = np.column_stack(ans)

# 保存到文本文件，格式为.txt
# np.savetxt(filename+"_DynaTopo.txt", data_to_save, fmt='%.2f')  # fmt指定保存时的数值格式


