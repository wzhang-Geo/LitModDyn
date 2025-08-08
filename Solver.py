#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 18:14:31 2021

@author: wzhang
"""

from __future__ import division
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix

from scipy.io import savemat
from scipy.io import loadmat

# npzfile = np.load('PYout.npz')
# # prfirst, etas1, etan1, xnum, ynum, xstp, ystp, RX1, RY1, RC1, bleft, bright, btop, bbottom
# prnorm,etas,etan,xnum,ynum,xstp,ystp,RX,RY,RC,bleft,bright,btop,bbottom \
#     = npzfile['prfirst'], npzfile['etas1'], npzfile['etan1'], \
#       npzfile['xnum'], npzfile['ynum'], npzfile['xstp'], npzfile['ystp'], \
#       npzfile['RX1'], npzfile['RY1'], npzfile['RC1'], \
#       npzfile['bleft'], npzfile['bright'], npzfile['btop'], npzfile['bbottom']

def main(prnorm,etas,etan,xnum,ynum,xstp,ystp,RX,RY,RC,bleft,bright,btop,bbottom):
    # Poisson-like equations koefficients
    xkf=1/xstp ** 2
    xkf2=2/xstp ** 2
    ykf=1/ystp ** 2
    ykf2=2/ystp ** 2
    xykf=1/(xstp*ystp)
    
    
    # Koefficient for scaling pressure
    pscale=2*etan[0,0]/(xstp+ystp)
    
    # Horizontal shift index
    ynum3=(ynum-1)*3
    
    
    # Creating matrix
    L=lil_matrix(((xnum-1)*(ynum-1)*3,(xnum-1)*(ynum-1)*3))
    R=np.zeros(((xnum-1)*(ynum-1)*3,1))
    
    # Solving of Stokes and continuity equations on nodes
    for i in range(ynum-1):
        for j in range(xnum-1):
            # Indexes for P,vx,vy
            ivx=(j * (ynum-1) + i) * 3
            ivy=ivx+1
            ipr=ivx+2
    
            # x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX
            if j < xnum - 2:
                # x-Stokes equation stensil
                #     +-------------------- -+----------------------+
                #     |                      |                      |
                #     |                      |                      |
                #     |                   vx[i-1,j]                 |
                #     |                      |                      |
                #     |                      |                      |
                #     +-----vy[i-1,j]---etas[i,j+1]---vy[i-1,j+1]---+
                #     |                      |                      |
                #     |                      |                      |
                # vx[i,j-1]  pr[i,j]      vx[i,j]     P[i,j+1]   vx[i,j+1]
                #     |     etan[i,j]        |       etan[i,j+1]    |
                #     |                      |                      |
                #     +------vy[i,j]---etas[i+1,j+1]---vy[i,j+1]----+
                #     |                      |                      |
                #     |                      |                      |
                #     |                   vx[i+1,j]                 |
                #     |                      |                      |
                #     |                      |                      |
                #     +-------------------- -+----------------------+
                # Right Part
                R[ivx, 0] = RX[i + 1, j + 1]
                # Computing Current x-Stokes coefficients
                # Central Vx node  Vx[i,j]
                L[ivx, ivx] = -xkf2 * (etan[i, j + 1] + etan[i, j]) - ykf * (etas[i + 1, j + 1] + etas[i, j + 1])
                # Left Vx node  Vx[i,j-1]
                if j > 0:
                    ivxleft = ivx - ynum3
                    L[ivx, ivxleft] = xkf2 * etan[i, j]
    
                # Right Vx node Vx[i,j+1]
                if j < xnum - 3:
                    ivxright = ivx + ynum3
                    L[ivx, ivxright] = xkf2 * etan[i, j + 1]
    
                # Top Vx node   vx[i-1,j]
                if i > 0:
                    ivxtop = ivx - 3
                    L[ivx, ivxtop] = ykf * etas[i, j + 1]
                else:
                    L[ivx, ivx] = L[ivx, ivx] + btop * ykf * etas[i, j + 1]
    
                # Bottom Vx node vx[i+1,j]
                if i < ynum - 2:
                    ivxbottom = ivx + 3
                    L[ivx, ivxbottom] = ykf * etas[i + 1, j + 1]
                else:
                    L[ivx, ivx] = L[ivx, ivx] + bbottom * ykf * etas[i + 1, j + 1]
    
                # Top Left Vy node  vy[i-1,j]
                if i > 0:
                    ivytopleft = ivx - 3 + 1
                    L[ivx, ivytopleft] = xykf * etas[i, j + 1]
    
                # Top Right Vy node  vy[i-1,j+1]
                if i > 0:
                    ivytopright = ivx - 3 + 1 + ynum3
                    L[ivx, ivytopright] = -xykf * etas[i, j + 1]
    
                # Bottom Left Vy node vy[i,j]
                if i < ynum - 2:
                    ivybottomleft = ivx + 1
                    L[ivx, ivybottomleft] = -xykf * etas[i + 1, j + 1]
    
                # Bottom Right Vy node  vy[i,j+1]
                if i < ynum - 2:
                    ivybottomright = ivx + 1 + ynum3
                    L[ivx, ivybottomright] = xykf * etas[i + 1, j + 1]
    
                # Left P node
                iprleft = ivx + 2
                L[ivx, iprleft] = pscale / xstp
                # Right P node
                iprright = ivx + 2 + ynum3
                L[ivx, iprright] = -pscale / xstp
    
            # Ghost Vx_parameter=0 used for numbering
            else:
                L[ivx, ivx] = 1
                R[ivx, 0] = 0
    
            # y-Stokes equation dSIGMAyy/dy+dSIGMAyx/dx-dP/dy=RY
            if i < ynum - 2:
                # y-Stokes equation stensil
                #     +-------------------- -+-------vy[i-1,j]------+----------------------+
                #     |                      |                      |                      |
                #     |                      |                      |                      |
                #     |                  vx[i,j-1]     P[i,j]    vx[i,j]                   |
                #     |                      |        etan[i,j]     |                      |
                #     |                      |                      |                      |
                #     +-----vy[i,j-1]---etas[i+1,j]---vy[i,j]--etas[i+1,j+1]---vy[i,j+1]---+
                #     |                      |                      |                      |
                #     |                      |                      |                      |
                #     |                  vx[i+1,j-1]  P[i+1,j]   vx[i+1,j]                 |
                #     |                      |      etan[i+1,j]     |                      |
                #     |                      |                      |                      |
                #     +----------------------+-------vy[i+1,j]------+----------------------+
                #
                # Right Part
                R[ivy, 0] = RY[i + 1, j + 1]
                # Computing Current y-Stokes coefficients
                # Central Vy node
                L[ivy, ivy] = -ykf2 * (etan[i + 1, j] + etan[i, j]) - xkf * (etas[i + 1, j + 1] + etas[i + 1, j])
                # Top Vy node
                if i > 0:
                    ivytop = ivy - 3
                    L[ivy, ivytop] = ykf2 * etan[i, j]
    
                # Bottom Vy node
                if i < ynum - 3:
                    ivybottom = ivy + 3
                    L[ivy, ivybottom] = ykf2 * etan[i + 1, j]
    
                # Left Vy node
                if j > 0:
                    ivyleft = ivy - ynum3
                    L[ivy, ivyleft] = xkf * etas[i + 1, j]
                else:
                    L[ivy, ivy] = L[ivy, ivy] + bleft * xkf * etas[i + 1, j]
    
                # Right Vy node
                if j < xnum - 2:
                    ivyright = ivy + ynum3
                    L[ivy, ivyright] = xkf * etas[i + 1, j + 1]
                else:
                    L[ivy, ivy] = L[ivy, ivy] + bright * xkf * etas[i + 1, j + 1]
    
                # Top left Vx node
                if j > 0:
                    ivxtopleft = ivy - 1 - ynum3
                    L[ivy, ivxtopleft] = xykf * etas[i + 1, j]
    
                # Bottom left Vx node
                if j > 0:
                    ivxbottomleft = ivy - 1 + 3 - ynum3
                    L[ivy, ivxbottomleft] = -xykf * etas[i + 1, j]
    
                # Top right Vx node
                if j < xnum - 2:
                    ivxtopright = ivy - 1
                    L[ivy, ivxtopright] = -xykf * etas[i + 1, j + 1]
    
                # Bottom right Vx node
                if j < xnum - 2:
                    ivxbottomright = ivy - 1 + 3
                    L[ivy, ivxbottomright] = xykf * etas[i + 1, j + 1]
    
                # Top P node
                iprtop = ivy + 1
                L[ivy, iprtop] = pscale / ystp
                # Bottom P node
                iprbottom = ivy + 1 + 3
                L[ivy, iprbottom] = -pscale / ystp
    
            # Ghost Vy_parameter=0 used for numbering
            else:
                L[ivy, ivy] = 1
                R[ivy, 0] = 0
    
    
            # Continuity equation dvx/dx+dvy/dy=RC
            if j > 0 or i > 0:
                # Continuity equation stensil
                #     +-----vy[i-1,j]--------+
                #     |                      |
                #     |                      |
                # vx[i,j-1]  pr[i,j]      vx[i,j]
                #     |                      |
                #     |                      |
                #     +------vy[i,j]---------+
                #
                # Right Part
                R[ipr, 0] = RC[i, j]
                # Computing Current Continuity coefficients
                # Left Vx node
                if j > 0:
                    ivxleft = ipr - 2 - ynum3
                    L[ipr, ivxleft] = -pscale / xstp
    
                # Right Vx node
                if j < xnum - 2:
                    ivxright = ipr - 2
                    L[ipr, ivxright] = pscale / xstp
    
                # Top Vy node
                if i > 0:
                    ivytop = ipr - 1 - 3
                    L[ipr, ivytop] = -pscale / ystp
    
                # Bottom Vy node
                if i < ynum - 2:
                    ivybottom = ipr - 1
                    L[ipr, ivybottom] = pscale / ystp
    
    
            # Pressure definition for the upper left node
            else:
                L[ipr, ipr] = 2 * pscale / (xstp + ystp)
                R[ipr, 0] = 2 * prnorm / (xstp + ystp)
    
    # Solve matrix
    S=spsolve(L.tocsr(), R)
    # savemat("PY_LRS.mat", {"L": L.todense(), "R": R, "S": S})
    # data = loadmat("matlab_LRS.mat")
    # S = data['S']
    
    
    # matfile = loadmat("PY_LRS.mat")
    # S = matfile['S'].T
    # matfile = loadmat("build/matfile_LRS.mat")
    # S = matfile['S']
    
    # Reload solution
    vx = np.zeros((ynum+1,xnum))
    vy = np.zeros((ynum,xnum+1))
    pr = np.zeros((ynum-1,xnum-1))
    resx = np.zeros((ynum+1,xnum))
    resy = np.zeros((ynum,xnum+1))
    resc = np.zeros((ynum-1,xnum-1))
    
    for i in range(ynum-1):
        for j in range(xnum - 1):
            # Indexes for P,vx,vy
            ivx = (j * (ynum - 1) + i) * 3
            ivy=ivx+1
            ipr=ivx+2
            # Reload Vx
            if j<xnum-2:
                vx[i+1,j+1]=S[ivx]
    
            # Reload Vy
            if i<ynum-2:
                vy[i+1,j+1]=S[ivy]
    
            # Reload P
            pr[i,j]=S[ipr]*pscale
    
    
    # Apply vx boundary conditions
    # Left Boundary
    vx[:,0]=0
    # Right Boundary
    vx[:,xnum-1]=0
    # Upper Boundary
    vx[0,:]=btop*vx[1,:]
    # Lower Boundary
    vx[ynum,:]=bbottom*vx[ynum-1,:]
    
    # Apply vy boundary conditions
    # Upper Boundary
    vy[0,:]=0
    # Lower Boundary
    vy[ynum-1,:]=0
    # Left Boundary
    vy[:,0]=bleft*vy[:,1]
    # Right Boundary
    vy[:,xnum]=bright*vy[:,xnum-1]
    
    # Computing residuals
    for i in range(ynum+1):
        for j in range(xnum+1):
            # x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX
            if j<xnum:
                # vx-Boundrary conditions 
                if i==0 or i==ynum or j==0 or j==xnum-1:
                    resx[i,j]=0
                # x-Stokes equation
                else:
                    # Computing Current x-Stokes residual
                    # dSIGMAxx/dx-dP/dx
                    resx[i,j]=RX[i,j]-(xkf2*(etan[i-1,j]*(vx[i,j+1]-vx[i,j])-etan[i-1,j-1]*(vx[i,j]-vx[i,j-1]))-(pr[i-1,j]-pr[i-1,j-1])/xstp)
                    # dSIGMAxy/dy
                    resx[i,j]=resx[i,j]-(etas[i,j]*((vx[i+1,j]-vx[i,j])*ykf+(vy[i,j+1]-vy[i,j])*xykf)-etas[i-1,j]*((vx[i,j]-vx[i-1,j])*ykf+(vy[i-1,j+1]-vy[i-1,j])*xykf))
                           
                
            # y-Stokes equation dSIGMAyy/dy+dSIGMAyx/dx-dP/dy=RY
            if i<ynum:
                # vy-Boundrary conditions 
                if i==0 or i==ynum-1 or j==0 or j==xnum:
                    resy[i,j]=0
                #y-Stokes equation
                else:
                    # Computing current residual
                    # dSIGMAyy/dy-dP/dy
                    resy[i,j]=RY[i,j]-(ykf2*(etan[i,j-1]*(vy[i+1,j]-vy[i,j])-etan[i-1,j-1]*(vy[i,j]-vy[i-1,j]))-(pr[i,j-1]-pr[i-1,j-1])/ystp)
                    # dSIGMAxy/dx
                    resy[i,j]=resy[i,j]-(etas[i,j]*((vy[i,j+1]-vy[i,j])*xkf+(vx[i+1,j]-vx[i,j])*xykf)-etas[i,j-1]*((vy[i,j]-vy[i,j-1])*xkf+(vx[i+1,j-1]-vx[i,j-1])*xykf))
    
                
            
                
            # Continuity equation dvx/dx+dvy/dy=RC
            if i<ynum-1 and j<xnum-1:
                # Computing current residual
                resc[i,j]=RC[i,j]-((vx[i+1,j+1]-vx[i+1,j])/xstp+(vy[i+1,j+1]-vy[i,j+1])/ystp)
    return(vx, resx, vy, resy, pr, resc)
             
                
# savemat("PY_reslut.mat", 
#         {"vx1": vx, "vy1": vy, "pr1": pr, 
#          "resx1": resx, "resy1": resy, "resc1": resc})



if __name__ == '__main__':
    main()

    # (vx1, resx1, vy1, resy1, pr1, resc1) \
    #     = main(prfirst, etas1, etan1, xnum, ynum, xstp, ystp, RX1, RY1, RC1, bleft, bright, btop, bbottom)




