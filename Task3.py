
from random import seed, random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import fft 
import cv2
import pylab 
import time

fileName =cv2.imread('test-1.jpg',0)
dim = fileName.shape[0]
kspace = [[0 for col in range(dim)] for row in range(dim)]
Ranges={'whiteMatter':{'T1':(760,1080),'T2':(61,100)},'CSF':{'T1':(800,2000),'T2':(110,2000)},'GrayMatter':{'T1':(109,2150),'T2':(61,109)}}
TE = 0.03
TR = 0.5
xdata = []
ydata =[]
dataxPlot=[]
datayPlot=[]
print(dim)
N = 1000
z = np.linspace(0, 3, N)
df = 42.6*3
spins = []
theta = []
x = 0
y = 0
RN=np.empty((dim,dim,dim,dim))
# K = 


def RandomNumbers():
    seed(1)
    for j in range(0,dim):
        for i in range(0,dim):
            for n in range(0,dim):
                for k in range(0,dim):
                    RN[j][i][n][k] =random()*0.5 +3
    return(RN)


def mapping(value):
    T1Range =(0.8,2)
    T2Range =(0.11,2)
    val1 = ((value -T1Range[0])/(T1Range[1]-T1Range[0]))*T1Range[0]
    val2 = ((value -T2Range[0])/(T2Range[1]-T2Range[0]))*T2Range[0]    
    return(val1,val2)


def getT1AndT2():
    dim = fileName.shape[0]
    T1 = np.empty((fileName.shape[0],fileName.shape[1]))
    T2 = np.empty((fileName.shape[0],fileName.shape[1]))
    rho = np.empty((fileName.shape[0],fileName.shape[1]))
    data = fileName
    for row in range (dim):
        for col in range (dim):
            if (np.mean(data[row][col]) == 0):
                T1[row][col] = 1.09
                T2[row][col] = 0.61 
                rho[row][col] = 85              
            elif(np.mean(data[row][col])==255):
                T1[row][col] = 0.76
                T2[row][col] = 0.61
                rho[row][col] = 70            
            else:
                T1[row][col],T2[row][col] = mapping(np.mean(data[row][col]))
                rho[row][col] = 85    
    return(T1,T2,rho)

def k_space():
    T1,T2,rho = getT1AndT2()
 
    for row in range(dim):
        for col in range(dim):
            spin = (rho[row][col])*(1-np.exp(-TR/T1[row][col]))*(np.exp(-TE/T2[row][col]))
            spins.append(spin)
            for n in range(dim):
                Gy = (-np.pi) + (((2*np.pi)/dim)*n)
                for i in range(dim):
                    Gx = 0 + (((2*np.pi)/dim)*i)
                    theta.append(df*(Gx + Gy +3*np.pi/180) )
            signal = np.sum(spins)
            kspace[row][col] =signal*np.cos(np.sum(theta)) + 1j*np.sin(np.sum(theta))
    image = fft.ifftn(kspace)
    image = fft.fftshift(image)
    k = np.real(image)    
    plt.subplot(121) 
    plt.imshow(fileName, cmap = 'gray',interpolation='nearest')
    plt.title('Original image')
    plt.subplot(122) 
    plt.imshow(k, cmap = 'gray',interpolation='nearest')
    plt.title('K sapce')
    plt.show() 

def NonUniformKspace():
    T1,T2,rho = getT1AndT2()
    RN = RandomNumbers()
    for row in range(dim):
        for col in range(dim):
            spin = (rho[row][col])*(1-np.exp(-TR/T1[row][col]))*(np.exp(-TE/T2[row][col]))
            spins.append(spin)
            for n in range(dim):
                Gy = (-np.pi) + (((2*np.pi)/dim)*n)
                for i in range(dim):
                    Gx = 0 + (((2*np.pi)/dim)*i)
                    theta.append(df*(Gx + Gy + (RN[row][col][n][i]*np.pi)/180))
            signal = np.sum(spins)
            kspace[row][col] =signal*np.cos(np.sum(theta)) + 1j*np.sin(np.sum(theta))
    image = fft.ifftn(kspace)
    image = fft.fftshift(image)
    k = np.real(image)  
    plt.subplot(121) 
    plt.imshow(fileName, cmap = 'gray',interpolation='nearest')
    plt.title('Original image')
    plt.subplot(122) 
    plt.imshow(k, cmap = 'gray',interpolation='nearest')
    plt.title('K sapce non-uniformity')
    plt.show() 

def zrot(phi):   
    Rz = [[np.cos(phi) ,-np.sin(phi), 0],[np.sin(phi) ,np.cos(phi), 0],[ 0, 0 ,1]]
    return(Rz)

def freeprecess(dT ,T1 ,T2 , df):
    phi = df
    E1 = np.exp(-dT/T1) 
    E2 = np.exp(-dT/T2) 
    Afp = np.dot([[E2,0,0],[0,E2,0],[0,0,E1]],zrot(phi))     
    Bfp = [0 ,0 ,1-E1]   
    return(Afp,Bfp)

def point():
    T1,T2,rho = getT1AndT2()
    T1 = T1[x][y]*1000
    T2 = T2[x][y]*1000
    Gx = ((2*np.pi)/dim)*x
    Gy = (-np.pi) + (((2*np.pi)/dim)*y)
    df =42.4*(3 + x*Gx + y*Gy)
    blochEquation(T1,T2,df)



def Magnetization(A,B):
    M = np.empty((N,3))    
    M[0,:] =np.array([1,0,0])   
    for k in range (N-1):
        M[k+1,:] = np.dot(A,M[k,:])+B
    return(M[:,0],M[:,1],M[:,2])


def blochEquation(T1,T2,df):
    dT =1    
    A,B = freeprecess(dT,T1,T2,df)
    xdata,ydata,zdata = Magnetization(A,B)
    # print(ydata)
    plot(xdata,ydata)

def plot(xdata,ydata):
   
    ax = pylab.gca(projection='3d')   
    ax.set_xlim(min(xdata), max(xdata))
    ax.set_ylim(min(ydata),max(ydata))
    ax.set_zlim(0, 5)
    for i in range(1,500):
        dataxPlot.append(xdata[i-1])
        datayPlot.append(ydata[i-1]) 
        pylab.plot(dataxPlot,datayPlot,z[:i],color ="red",linewidth=1.5)
        pylab.draw()
        pylab.pause(1e-117)        
     
    plt.show() 

# k_space()
# point()
# blochEquation()
# getT1AndT2()
# NonUniformKspace()
# points = RandomNumbers()
# print(points[0][1][1][1])