#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# """
# Created on Thu Jun 17 18:14:31 2021
# V4.0 对代码进行了简化 + 修改了MIC 插值函数 - 2023-08-24

# @author: wzhang
# """
#%%
from __future__ import division # Changing the Division Operator
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
from scipy.interpolate import griddata
from datetime import datetime

# import Solver_multigrid as Solver
import Solver as Solver
# from main_lit import BC, model_set, model_extending_set,GetETA,BC_extending #, round_up

starttime = datetime.now()
print(f'Start time:  {starttime}')



prfirst=0.0         # Pressure in the upermost, leftmost (first) cell
bleft, bright, btop, bbottom=1, 1, 1, 1     # [1=free slip -1=no slip ] are implemented from ghost nodes

#%% Model initialization
# from Ex_falling_block import *     # This is examples from Gerya. (2019)
# savename = 'Sinking_Clinder_Free_Slip_'


# from Ex_falling_block_Ana import *   
# savename = 'Sinking_Clinder_Free_Slip_Ana_'

# from Zhang_Northern import * 
# savename = 'Zhang_Northern_a'


# from Zhang_Southern import * 
# savename = 'Zhang_Southern_20240317'


from Zhang_Northern_ext import * 
savename = 'N_Ext_'

# from Zhang_Southern_ext import * 
# savename = 'S_ext_'

print('The name of save file:  ' + savename)

#%% MIC method. Refer from Gerya. (2019)
import Interpolation
typ1 = Interpolation.Griddata((MX, MY), MI, (NX, NY), method='MIC', size=0.5)
rho1 = Interpolation.Griddata((MX, MY), MRHO, (NX, NY), method='MIC')
etas1 = Interpolation.Griddata((MX, MY), META, (NX, NY), method='MIC', size=0.5)
etan1 = Interpolation.Griddata((MX, MY), META, (CX, CY), method='MIC_node_in_center')

# Right parts of equations
RX1=np.zeros((ynum+1,xnum))
RY1=np.zeros((ynum,xnum+1))
RC1=np.zeros((ynum-1,xnum-1))
MVX = np.zeros((mynum,mxnum))
MVY = np.zeros((mynum,mxnum))
MV = np.zeros((mynum,mxnum))

## Compute Stress and strain rate components
EXX=np.zeros((ynum-1,xnum-1)) # Strain rate EPSILONxx, 1/s
SXX=np.zeros((ynum-1,xnum-1)) # deviatoric stress SIGMAxx, Pa
EYY=np.zeros((ynum-1,xnum-1)) # Strain rate EPSILONyy, 1/s
SYY=np.zeros((ynum-1,xnum-1)) # deviatoric stress SIGMAyy, Pa
EXY=np.zeros((ynum,xnum)) # Strain rate EPSILONxy, 1/s
SXY=np.zeros((ynum,xnum)) # deviatoric stress SIGMAxy, Pa
# second invariant of the Strain rate, 1/s, deviatoric stress, Pa
# EII=np.zeros((ynum-1,xnum-1)) # EII=((EXX'2 + EYY'2 + EYX'2 + EXY'2)/2)^0.5
# SII=np.zeros((ynum-1,xnum-1)) # SII=((SXX'2 + SYY'2 + SYX'2 + SXY'2)/2)^0.5


MSXX=np.zeros((mynum,mxnum))  # SIGMAxx - deviatoric normal stress, Pa
MSYY=np.zeros((mynum,mxnum))
MSXY=np.zeros((mynum,mxnum))  # SIGMAyy - shear stress, Pa

MEXX=np.zeros((mynum,mxnum))  # EPSILONxx - normal strain rate, 1/s
MEYY=np.zeros((mynum,mxnum))
MEXY=np.zeros((mynum,mxnum))  # EPSILONyy - shear strain rate, 1/s

MPR=np.zeros((mynum,mxnum))   # Pressure, Pa

# MEII=np.ones((mynum,mxnum))   # EPSILONyy - shear strain rate, 1/s

#%%
stepmax=5
for ntimestep in range(stepmax):
    
    
    ######### <<<<< META  
    MEII0 = MEII.copy()
    Rheology = 'NE20_Wet'
    print('Rheology parameters: ', Rheology)
    if ntimestep > -1:
        import main_lit
        print('ntimestep is ', ntimestep)
        for xm in range(mxnum):
            for ym in range(mynum):
                if MI[ym,xm]<0:	     # sticky air
                    # MI[ym,xm] = -10
                    MRHO[ym,xm] =1000 # 1000
                    MTK[ym,xm] = 0            
                    META[ym,xm] = 1.0e18  
                elif MI[ym,xm]<11:    # sediment & crust
                    META[ym,xm] = 1.0e23
                elif MI[ym,xm]<91:      # Lithopshere mantle
                    META[ym,xm] = 1.0e24 
                else:				# Sublithopshere mantle
                    # META[ym,xm] = 1.0e19 # mantle
                    # META[ym,xm] = main_lit.GetETA(MTK[ym,xm],MPR_LitMod[ym,xm],MEII[ym,xm],method='KW93_Dry')
                    # META[ym,xm] = main_lit.GetETA(MTK[ym,xm],MPR_LitMod[ym,xm],MEII[ym,xm],method='KW93_Wet')
                    META[ym,xm] = main_lit.GetETA(MTK[ym,xm],MPR_LitMod[ym,xm],MEII[ym,xm],Rheology)


        # MRHO[MI<0] = 1000
        # MTK[MI<0] = 0.0
        # META[MI<0] = 1.0e18 
        
        # META[(MI>=0) & (MI<11)] = 1.0e21
        # META[(MI>=11) & (MI<91)] = 1.0e22 
        
        # META[MI>91] = main_lit.GetETA(MTK[MI>91],MPR_LitMod[MI>91],MEII[MI>91])
        
        
    etas1 = Interpolation.Griddata((MX, MY), META, (NX, NY), method='MIC', size=0.5)
    etan1 = Interpolation.Griddata((MX, MY), META, (CX, CY), method='MIC_node_in_center')
    ######### META  >>>>>
    
    # Grid points cycle

    for i in range(ynum):
        for j in range(xnum):
            # Right part of x-Stokes Equation
            if j>0 and i>0 and j<xnum - 1:
                RX1[i,j ]=0

            # Right part of y-Stokes Equation
            if j>0 and i>0 and i<ynum - 1:
                RY1[i,j ]=-g * (rho1[i,j ]+rho1[i,j-1 ])/2
    
    (vx1,resx1,vy1,resy1,pr1,resc1)=Solver.main(prfirst,etas1,etan1,xnum,ynum,xstp,ystp,RX1,RY1,RC1,bleft,bright,btop,bbottom)
    
    # Compute EPSILONxy, SIGMAxy in basic nodes
    for xn in range(xnum):
        for yn in range(ynum):
            # EXY=0.5(dvx/dy+dvy/dx)
            EXY[yn,xn]=0.5*((vx1[yn+1,xn]-vx1[yn,xn])/(VX_Y[yn+1,xn]-VX_Y[yn,xn])
                            +(vy1[yn,xn+1]-vy1[yn,xn])/xstp)
            
    
    # Compute EPSILONxx, SIGMA'xx in pressure nodes
    # deviatoric stress ij = 2*eta*strain rate ij
    for xn in range(xnum-1):
        for yn in range(ynum-1): 
            # EXX=dvx/dx        EYY=dvy/dy
            EXX[yn,xn]=(vx1[yn+1,xn+1]-vx1[yn+1,xn])/xstp
            EYY[yn,xn]=(vy1[yn+1,xn+1]-vy1[yn,xn+1])/(VY_Y[yn+1,xn+1]-VY_Y[yn,xn+1])

            # # SII
            # SII[yn,xn]=1/2*(SXX[yn,xn]^2 + SYY[yn,xn]^2 + SXY[yn,xn]^2 + SXY[yn,xn]^2)
            # SII[yn,xn]=SII[yn,xn]^0.5
            # # EII
            # EII[yn,xn]=1/2*(EXX[yn,xn]^2 + EYY[yn,xn]^2 + EXY[yn,xn]^2 + EXY[yn,xn]^2)
            # EII[yn,xn]=EII[yn,xn]^0.5
            # # SII
            # SII[yn,xn]=1/2*(SXX[yn,xn]^2 + SYY[yn,xn]^2 + SXY[yn,xn]^2 + SXY[yn,xn]^2)
            # SII[yn,xn]=SII[yn,xn]^0.5
    # SXY,SXX,SYY
    SXY=2*etas1*EXY
    SXX=2*etan1*EXX            
    SYY=2*etan1*EYY
    
    
    
    for xm in range(mxnum):# 
        for ym in range(mynum):
            #  xn    V(xn,yn)--------------------V(xn+1,yn)
            #           ?           **                  ?
            #           ?           ?                  ?
            #           ?          dy                  ?
            #           ?           ?                  ?
            #           ?           v                  ?
            #           ?<----dx--->o Mrho(xm,ym)       ?
            #           ?                              ?
            #           ?                              ?
            #  xn+1  V(xn,yn+1)-------------------V(xn+1,yn+1)
            
            ## Computing strain rate and pressure for markers               
            
            xn=int(MX[ym,xm]/xstp-0.5)
            yn=int(MY[ym,xm]/ystp-0.5)

            # yn=int(round(MY[ym,xm]/ystp))
            # yn=sum(MY[ym,xm]-CY[:,0]>0)-1
            # xn=double(int16(MX(ym,xm)./xstp-0.5))+1
            # yn=double(int16((MY(ym,xm)+ystp/2.0)./ystp-0.5))+1
            # Check vertical index for upper left Pr-node 
            # It must be between 1 and ynum (see picture for staggered grid)
            if xn<0:
                xn=0
            
            if xn>xnum-3:
                xn=xnum-3
            
            if yn<0:
                yn=0
            
            if yn>ynum-3:
                yn=ynum-3
            
            # Define and check normalized distances from marker to the pressure-node
            dx=(MX[ym,xm]-CX[yn,xn])/xstp
            
            if yn==ynum-2:
                dy=(MY[ym,xm]-CY[yn,xn])/(CY[yn,0]-CY[yn-1,0])
            else:
                dy=(MY[ym,xm]-CY[yn,xn])/(CY[yn+1,0]-CY[yn,0])
            # dy_sum[ym,xm]=dy
            # Calculate Marker velocity from four surrounding nodes
            MSXX[ym,xm]=0
            MSXX[ym,xm]=MSXX[ym,xm]+(1.0-dx)*(1.0-dy)*SXX[yn,xn]
            MSXX[ym,xm]=MSXX[ym,xm]+(1.0-dx)*dy*SXX[yn+1,xn]
            MSXX[ym,xm]=MSXX[ym,xm]+dx*(1.0-dy)*SXX[yn,xn+1]
            MSXX[ym,xm]=MSXX[ym,xm]+dx*dy*SXX[yn+1,xn+1] 
            
            MEXX[ym,xm]=0
            MEXX[ym,xm]=MEXX[ym,xm]+(1.0-dx)*(1.0-dy)*EXX[yn,xn]
            MEXX[ym,xm]=MEXX[ym,xm]+(1.0-dx)*dy*EXX[yn+1,xn]
            MEXX[ym,xm]=MEXX[ym,xm]+dx*(1.0-dy)*EXX[yn,xn+1]
            MEXX[ym,xm]=MEXX[ym,xm]+dx*dy*EXX[yn+1,xn+1] 
            
            MSYY[ym,xm]=0
            MSYY[ym,xm]=MSYY[ym,xm]+(1.0-dx)*(1.0-dy)*SYY[yn,xn]
            MSYY[ym,xm]=MSYY[ym,xm]+(1.0-dx)*dy*SYY[yn+1,xn]
            MSYY[ym,xm]=MSYY[ym,xm]+dx*(1.0-dy)*SYY[yn,xn+1]
            MSYY[ym,xm]=MSYY[ym,xm]+dx*dy*SYY[yn+1,xn+1] 
            
            MEYY[ym,xm]=0
            MEYY[ym,xm]=MEYY[ym,xm]+(1.0-dx)*(1.0-dy)*EYY[yn,xn]
            MEYY[ym,xm]=MEYY[ym,xm]+(1.0-dx)*dy*EYY[yn+1,xn]
            MEYY[ym,xm]=MEYY[ym,xm]+dx*(1.0-dy)*EYY[yn,xn+1]
            MEYY[ym,xm]=MEYY[ym,xm]+dx*dy*EYY[yn+1,xn+1] 
            
            MPR[ym,xm]=0
            MPR[ym,xm]=MPR[ym,xm]+(1.0-dx)*(1.0-dy)*pr1[yn,xn]
            MPR[ym,xm]=MPR[ym,xm]+(1.0-dx)*dy*pr1[yn+1,xn]
            MPR[ym,xm]=MPR[ym,xm]+dx*(1.0-dy)*pr1[yn,xn+1]
            MPR[ym,xm]=MPR[ym,xm]+dx*dy*pr1[yn+1,xn+1] 
            
            
            ## Computing strain rate and pressure for markers   
            xn=int(MX[ym,xm]/xstp)
            yn=int(MY[ym,xm]/ystp)

            # yn=int(round(MY[ym,xm]/ystp))
            # yn=sum(MY[ym,xm]-CY[:,0]>0)-1
            # xn=double(int16(MX(ym,xm)./xstp-0.5))+1
            # yn=double(int16((MY(ym,xm)+ystp/2.0)./ystp-0.5))+1
            # Check vertical index for upper left Pr-node 
            # It must be between 1 and ynum (see picture for staggered grid)
            if xn<0:
                xn=0
            
            if xn>xnum-1:
                xn=xnum-1
            
            if yn<0:
                yn=0
            
            if yn>ynum-1:
                yn=ynum-1
            
            # Define and check normalized distances from marker to the pressure-node
            dx=(MX[ym,xm]-CX[yn,xn])/xstp
            dy=(MY[ym,xm]-NY[yn,xn])/(NY[yn+1,0]-NY[yn,0])

            # Calculate Marker velocity from four surrounding nodes
            MSXY[ym,xm]=0
            MSXY[ym,xm]=MSXY[ym,xm]+(1.0-dx)*(1.0-dy)*SXY[yn,xn]
            MSXY[ym,xm]=MSXY[ym,xm]+(1.0-dx)*dy*SXY[yn+1,xn]
            MSXY[ym,xm]=MSXY[ym,xm]+dx*(1.0-dy)*SXY[yn,xn+1]
            MSXY[ym,xm]=MSXY[ym,xm]+dx*dy*SXY[yn+1,xn+1]
            
            MEXY[ym,xm]=0
            MEXY[ym,xm]=MEXY[ym,xm]+(1.0-dx)*(1.0-dy)*EXY[yn,xn]
            MEXY[ym,xm]=MEXY[ym,xm]+(1.0-dx)*dy*EXY[yn+1,xn]
            MEXY[ym,xm]=MEXY[ym,xm]+dx*(1.0-dy)*EXY[yn,xn+1]
            MEXY[ym,xm]=MEXY[ym,xm]+dx*dy*EXY[yn+1,xn+1] 
            
            
            # MEII[ym,xm]=(MSXX[ym,xm]**2+MSXY[ym,xm]**2)**0.5
            MEII[ym,xm]=0
            MEII[ym,xm]=1/2*(MEXX[ym,xm]**2 + MEYY[ym,xm]**2 + MEXY[ym,xm]**2 + MEXY[ym,xm]**2)
            MEII[ym,xm]=MEII[ym,xm]**0.5
            
            
            
##############################            
            xn=int(MX[ym,xm]/xstp-0.5)
            yn=int(MY[ym,xm]/ystp-0.5)



            # yn=int(round(MY[ym,xm]/ystp))
            # yn=sum(MY[ym,xm]-CY[:,0]>0)-1
            # xn=double(int16(MX(ym,xm)./xstp-0.5))+1
            # yn=double(int16((MY(ym,xm)+ystp/2.0)./ystp-0.5))+1
            # Check vertical index for upper left Pr-node 
            # It must be between 1 and ynum (see picture for staggered grid)
            if xn<0:
                xn=0
            
            if xn>xnum-3:
                xn=xnum-3
            
            if yn<0:
                yn=0
            
            if yn>ynum-3:
                yn=ynum-3
            
            # Define and check normalized distances from marker to the pressure-node
            dx=(MX[ym,xm]-CX[yn,xn])/xstp
            
            if yn==ynum-2:
                dy=(MY[ym,xm]-CY[yn,xn])/(CY[yn,0]-CY[yn-1,0])
            else:
                dy=(MY[ym,xm]-CY[yn,xn])/(CY[yn+1,0]-CY[yn,0])
            # dy_sum[ym,xm]=dy           
            
            ## Computing strain rate and pressure for markers   
            xn=int(MX[ym,xm]/xstp)
            yn=int(MY[ym,xm]/ystp)

            # yn=int(round(MY[ym,xm]/ystp))
            # yn=sum(MY[ym,xm]-CY[:,0]>0)-1
            # xn=double(int16(MX(ym,xm)./xstp-0.5))+1
            # yn=double(int16((MY(ym,xm)+ystp/2.0)./ystp-0.5))+1
            # Check vertical index for upper left Pr-node 
            # It must be between 1 and ynum (see picture for staggered grid)
            if xn<0:
                xn=0
            
            if xn>xnum-1:
                xn=xnum-1
            
            if yn<0:
                yn=0
            
            if yn>ynum-1:
                yn=ynum-1
            

            
            ##
            # Define indexes for upper left node in the VX-cell where the marker is
            # VX-cells are displaced upward for 1/2 of vertical gridstep
            # !!! SUBTRACT 0.5 since int16(0.5)=1
            xn=int(MX[ym,xm]/xstp)
            yn=int(MY[ym,xm]/ystp+0.5)

            # yn=int(round((MY[ym,xm]+ystp/2.0)/ystp-0.5))
            # yn=sum(MY[ym,xm]-CY[:,0]>0)-1
            # Check vertical index for upper left VX-node 
            # It must be between 1 and ynum (see picture for staggered grid)
            if xn<0:
                xn=0
            
            if xn>xnum-2:
                xn=xnum-2
            
            if yn<0:
                yn=0
            
            if yn>ynum-1:
                yn=ynum-1
            
            # Define and check normalized distances from marker to the upper left VX-node
            # dx=MX[ym,xm]/xstp-xn
            dx=(MX[ym,xm]-VX_X[yn,xn])/xstp
            dy=(MY[ym,xm]-VX_Y[yn,xn])/(VX_Y[yn+1,0]-VX_Y[yn,0])
            # dy_sum[ym,xm]=dy
            # dy=(MY[ym,xm]+ystp[yn]/2.0)/(NY[yn+1,0]-NY[yn,0])-yn
            # dy=(MY[ym,xm]+ystp[yn]/2.0)/ystp[yn]-yn

            # Calculate Marker velocity from four surrounding nodes
            MVX[ym,xm]=0
            MVX[ym,xm]=MVX[ym,xm]+(1.0-dx)*(1.0-dy)*vx1[yn,xn]
            MVX[ym,xm]=MVX[ym,xm]+(1.0-dx)*dy*vx1[yn+1,xn]
            MVX[ym,xm]=MVX[ym,xm]+dx*(1.0-dy)*vx1[yn,xn+1]
            MVX[ym,xm]=MVX[ym,xm]+dx*dy*vx1[yn+1,xn+1]

            # Define indexes for upper left node in the VY-cell where the marker is
            # VY-cells are displaced leftward for 1/2 of horizontal gridstep
            # !!! SUBTRACT 0.5 since int16(0.5)=1            
                
            # xn=int(round_up((MX[ym,xm]+xstp/2.0)/xstp-0.5))
            xn=int(MX[ym,xm]/xstp+0.5)
            yn=int(MY[ym,xm]/ystp)

            # yn=int(round(MY[ym,xm]/ystp-0.5))
            # yn=sum(MY[ym,xm]-NY[:,0]>0)-1
            # Check horizontal index for upper left VY-node 
            # It must be between 1 and xnum (see picture for staggered grid)
            if xn<0:
                xn=0
            
            if xn>xnum-1:
                xn=xnum-1
            
            if yn<0:
                yn=0
            
            if yn>ynum-2:
                yn=ynum-2
            
            # Define and check normalized distances from marker to the upper left VX-node
            # dx=(MX[ym,xm]+xstp/2.0)/xstp-xn
            
            dx=(MX[ym,xm]-VY_X[yn,xn])/xstp
            dy=(MY[ym,xm]-VY_Y[yn,xn])/(VY_Y[yn+1,0]-VY_Y[yn,0])
            # dy=MY[ym,xm]/ystp[yn]-yn
            # dy_sum[ym,xm]=dy
            # Calculate Marker velocity from four surrounding nodes
            MVY[ym,xm]=0
            MVY[ym,xm]=MVY[ym,xm]+(1.0-dx)*(1.0-dy)*vy1[yn,xn]
            MVY[ym,xm]=MVY[ym,xm]+(1.0-dx)*dy*vy1[yn+1,xn]
            MVY[ym,xm]=MVY[ym,xm]+dx*(1.0-dy)*vy1[yn,xn+1]
            MVY[ym,xm]=MVY[ym,xm]+dx*dy*vy1[yn+1,xn+1]

            # Displacing Marker according to its velocity
            
            # MV[ym,xm]=(MVX[ym,xm]**2 + MVY[ym,xm]**2)**0.5
    
    
    
    MV = (MVX ** 2+ MVY **2)**0.5
    
    filename = savename + str(ntimestep)
    
    savemat(filename + '.mat', {
            "prfirst": prfirst,'g': g,
            "etas1": etas1,
            "etan1": etan1,
            "xnum": xnum,"ynum": ynum,
            "mxnum": mxnum,"mynum": mynum,
            "xstp": xstp,"ystp": ystp,
            "RX1": RX1,"RY1": RY1,"RC1": RC1,
            "bleft": bleft,"bright": bright,"btop": btop,"bbottom": bbottom,
            "typ1": typ1, "rho1": rho1,
            "vx1": vx1, "vy1": vy1, "pr1": pr1,
            # "resx0": resx0, "resy0": resy0, "resc0": resc0,
            'MX': MX, 'MY': MY, 'MI':MI, 'MRHO':MRHO, 'META':META, 
            'MVX':MVX, 'MVY':MVY,'MV':MV,
            'NX': NX, 'NY': NY, 'CX':CX, 'CY':CY, 'VX_X':VX_X, 'VX_Y':VX_Y, 'VY_X':VY_X,'VY_Y':VY_Y,
            # 'SXX': SXX, 'SYY': SYY,'SXY': SXY,
            # 'EXX': EXX, 'EYY': EYY, 'EXY': EXY,
            # 'MSXX': MSXX, 'MSYY': MSYY, 'MSXY': MSXY, 
            # 'MEXX': MEXX, 'MEYY': MEYY, 'MEXY': MEXY,
            'MTK':MTK,
            # 'MEII':MEII,'MEII0':MEII0,'MPR': MPR,
            'file_in': file_in,
            # 'BC_Ext':BC_Ext,
            'measured_topo': measured_topo

        })

endtime = datetime.now()
print(f'Endtime:    {endtime}')
print('Total Time:  %f s' %((endtime - starttime).total_seconds()))
