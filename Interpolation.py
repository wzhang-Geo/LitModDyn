"""
Interpolation of physical parameters between markers and nodes

From (MX, MY) TO (NX, NY)

eg.
typ1 = Solver.Griddata((MX, MY), MI, (NX, NY), method='nearest')

"""
import numpy as np
from scipy.interpolate import griddata
from scipy.io import savemat

# 注意这里必须是规则网格下的方法
# size表示有效粒子的最大距离，不能超过1 (一个网格大小)
def marker2cell(MX,MY,MI,xi,size):
    if size>1:
        raise ValueError("!!! The size is too big." )
    
    X = xi[0]
    Y = xi[1]
    
    xstp=X[0,1]-X[0,0]
    ystp=Y[1,0]-Y[0,0]
        
    ynum,xnum=np.shape(X)   
    mynum,mxnum=np.shape(MX)  
    wights=np.zeros((ynum,xnum))

    ip = np.zeros((ynum,xnum))
    

        
    for xm in range(mxnum):
        for ym in range(mynum):
            
            # Define indexes for upper left node in the cell where the marker is
            xn = int(MX[ym,xm]/xstp)
            yn = int(MY[ym,xm]/ystp)
            
            dx=MX[ym,xm]/xstp-xn
            dy=MY[ym,xm]/ystp-yn
            
            if dx<=size and dy<=size:
                ip[yn,xn]=ip[yn,xn]+(1.0-dx)*(1.0-dy)*MI[ym,xm]
                wights[yn,xn]=wights[yn,xn]+(1.0-dx)*(1.0-dy)
                
            if dx<=size and dy>=1-size:
                ip[yn+1,xn]=ip[yn+1,xn]+(1.0-dx)*dy*MI[ym,xm]
                wights[yn+1,xn]=wights[yn+1,xn]+(1.0-dx)*dy
                
            if dx>=1-size and dy<=size:
                ip[yn,xn+1]=ip[yn,xn+1]+dx*(1.0-dy)*MI[ym,xm]
                wights[yn,xn+1]=wights[yn,xn+1]+dx*(1.0-dy)
                
            if dx>=1-size and dy>=1-size:
                ip[yn+1,xn+1]=ip[yn+1,xn+1]+dx*dy*MI[ym,xm]
                wights[yn+1,xn+1]=wights[yn+1,xn+1]+dx*dy
    
    for j in range(ynum):
        for i in range(xnum):
            
            ip[j,i] = ip[j,i]/wights[j,i]

    return ip

# 这个是针对网格中心插值
def marker2cell_c(MX,MY,MI,xi,size):   
    X = xi[0]
    Y = xi[1]
    
    xstp=X[0,1]-X[0,0]
    ystp=Y[1,0]-Y[0,0]
        
    ynum,xnum=np.shape(X)   
    mynum,mxnum=np.shape(MX)  
    wights=np.zeros((ynum,xnum))

    ip = np.zeros((ynum,xnum))
    
        
    for xm in range(mxnum):
        for ym in range(mynum):
            
            # Define indexes for upper left node in the cell where the marker is
            xn = int(MX[ym,xm]/xstp)
            yn = int(MY[ym,xm]/ystp)
            
            dx=MX[ym,xm]/xstp-xn
            dy=MY[ym,xm]/ystp-yn
            
            ip[yn,xn]=ip[yn,xn]+(1.0-abs(0.5-dx))*(1.0-abs(0.5-dy))*MI[ym,xm]
            wights[yn,xn]=wights[yn,xn]+(1.0-abs(0.5-dx))*(1.0-abs(0.5-dy));
    
    for j in range(ynum):
        for i in range(xnum):
            
            ip[j,i] = ip[j,i]/wights[j,i]

    return ip


def Griddata(points, values, xi, method='linear', size=1):
    MX = points[0]
    MY = points[1]
    MI = values    
   
    if method in ('nearest', 'linear', 'cubic'):
        ip = griddata((MX.flatten(), MY.flatten()), MI.flatten(), xi, method)
        return ip
    elif method in ('MIC'):
        # marker-in-cell (MIC) technique
        return marker2cell(MX,MY,MI,xi,size)
    
    elif method in ('MIC_node_in_center'):
        # marker-in-cell (MIC) technique for center of nodes
        return marker2cell_c(MX,MY,MI,xi,size)
    
    else:
        raise ValueError("Unknown interpolation method %r " % (method))


# This is examples from Gerya. (2019)
# Exercise 8.4. --- Stokes_Continuity_Markers.m
# Exercise 20.1. --- Variable_viscosity_block.m
if __name__ == "__main__":
    
    # Model size, m
    xsize=500000
    ysize=500000
       
    # Defining resolution
    xnum=55
    ynum=51
    
    # Defining gridsteps
    xstp=xsize/(xnum-1)
    ystp=ysize/(ynum-1)
      
    # Defining number of markers and steps between them in the horizontal and vertical direction
    xmx=5 #number of markers per cell in horizontal direction
    ymy=5 #number of markers per cell in vertical direction
    mxnum=(xnum-1)*xmx #total number of markers in horizontal direction
    mynum=(ynum-1)*ymy #total number of markers in vertical direction
    mxstep=xsize/mxnum #step between markers in horizontal direction   
    mystep=ysize/mynum #step between markers in vertical direction
    
    # Creating nodes, center & markers arrays
    x0 = np.arange(0 * xstp, xsize + xstp, xstp)
    y0 = np.arange(0 * ystp, ysize + ystp, ystp)
    NX,NY = np.meshgrid(x0, y0)
    
    x2 = np.arange(0.5 * xstp, xsize, xstp)
    y2 = np.arange(0.5 * ystp, ysize, ystp)
    CX,CY = np.meshgrid(x2, y2)
    del x0, y0, x2, y2
    
    # Creating markers arrays
    MX=np.zeros((mynum,mxnum)) # X coordinate
    MY=np.zeros((mynum,mxnum)) # Y coordinate
    MI=np.zeros((mynum,mxnum)) # Type
    MRHO=np.zeros((mynum,mxnum)) # Density
    META=np.zeros((mynum,mxnum)) # viscosity
    
    # Defining intial position of markers
    # Defining lithological structure of the model
    for xm in range(mxnum):
        for ym in range(mynum):
            MX[ym,xm]=xm*mxstep+mxstep/2
            MY[ym,xm]=ym*mystep+mystep/2
            MI[ym,xm]=1
            MRHO[ym,xm]=3200.0
            META[ym,xm]=1e+21
            # Density, viscosity structure definition for block
            # Relative distances for the marker inside the grid
            dx=MX[ym,xm]/xsize
            dy=MY[ym,xm]/ysize
            if( dx>=0.4 and dx<=0.6 and dy>=0.1 and dy<=0.3):
                MI[ym,xm]=2
                MRHO[ym,xm]=3300.0
                META[ym,xm]=1e+27

    
    typ1 = Griddata((MX, MY), MI, (NX, NY), method='MIC', size=0.5)
    
    rho1 = Griddata((MX, MY), MRHO, (NX, NY), method='MIC')
    
    etas1 = Griddata((MX, MY), META, (NX, NY), method='MIC', size=0.5)
    
    etan1 = Griddata((MX, MY), META, (CX, CY), method='MIC_node_in_center')
    
    # etan1 = np.log10(etan1)
    
    savemat('py.mat', {
            "typ1": typ1,
            "rho1": rho1,
            "etas1": etas1,
            "etan1": etan1

        })


