#This program creates the potential energy of a space above the 13 atoms used in the bowl of a corner of MOF-5 using Lennard-Jones
import os
import numpy as np
import math
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt



def magnitude(r,rspace):
    #3D distance formula
    return np.sqrt((r[0]-rspace[0])**2+(r[1]-rspace[1])**2+(r[2]-rspace[2])**2)
#The standard Leonard-Jones 6-12
def LJ6_12(sigma,e,r,rspace,q):
    #charge of e-
    #Q=0.48 #quatrupole moment 
    k=320000. #e^2/a0=1Hartree=320000K
    #get magnitude of r
    R=magnitude(r,rspace)
    #these are both in the same units so its okay even if use in Angstroms
    R=2*R #1 a0=0.5 A
    sigma=2*sigma
    #Potential calculation
    if (R==0.):
        #if going to inf. make 10000 instead
        return(10000.)
    else:
        LJ = 4*e*((sigma/R)**12-(sigma/R)**6)
        #Uquad=k*((Q*q)/(R**3)) Not using quadrupole moment
        U = LJ
        if (U>=10000.):
            #max out at 10000
            return(10000.)
        else:
            return(U)

def potential(Lz,L,Nz,N,r,e,sigma,q):
    Div=L/(N-1)
    Divz=Lz/(Nz-1)
    U=np.empty([Nz,N,N])
    for z in range (0,Nz):
        #Z coordinate in 3D space
        Uz=z*(Divz)+1.5#-.5*L #(z+1)(L/N)-0.5*L #From 0 to 4
        for y in range(0,N):
            #Y coordinate in 3D space
            Uy=y*(Div)-.5*L #From -2 to 2
            for x in range(0,N):
                #X coordinate in 3D space
                Ux=x*(Div)-.5*L #From -2 to 2
                #create vector to point in 3D space
                rspace=[Ux,Uy,Uz]
                #calculate potential for point
                U[z][y][x]=LJ6_12(sigma,e,r,rspace,q)
    print "Minimum of run: ", np.min(U)
    return U


def main():
    #Zn:0-2, O1:3, O2:4-9, C1:10-12
    r=[[-1.855325180842072, 0.0, -10.563315995147232], [0.9276625904210358, -1.6067587388901914, -10.563315995147232],
        [0.9276625904210358, 1.6067587388901914, -10.563315995147232],[0.0, 0.0, -11.219359106027404],
        [-2.2071535575638297, -1.575575329839865, -9.474260182374842], [2.468065039999284, -1.1236633859835425, -9.474260182374842],
        [-0.26091148243545437, 2.6992387158234075, -9.474260182374842], [2.468065039999284, 1.1236633859835425, -9.474260182374842],
        [-0.2609114824354546, -2.6992387158234075, -9.474260182374842], [-2.2071535575638292, 1.575575329839865, -9.474260182374842],
        [-1.46152887986063, -2.53144227664784, -9.152445142328544], [2.92305775972126, 0.0, -9.152445142328544],
        [-1.4615288798606296, 2.53144227664784, -9.152445142328544]]
    sigma=[2.4616,3.118,3.118,3.431,3.431,3.431,2.571]
    #Zn,O1,O2,C1,C2,C3,H
    #epsilon
    e=[62.3993,30.19,30.19,52.84,52.84,52.84,22.14]
    #not using adjusted charges
    #q=[1.8529,-2.2568,-1.0069,1.0982,-0.1378,-0.0518,0.1489]
    #divisions of x and y
    N=41
    #length of x and y dimensions in Angstroms
    L=6.
    #divisions of z
    Nz=41
    #length of x and y dimensions in Angstroms
    Lz=6.
    #Make an empty array to as blank potential
    U=np.empty([Nz,N,N])
    #Loop through the atoms
    for i in range(len(r)):
        r[i][2]=r[i][2]+11.219359106027404 #Translate so center O is at (0,0,0)
        if i<=2:
            U+=potential(Lz,L,Nz,N,r[i],e[0],sigma[0],q[0])
        elif i==3:
            U+=potential(Lz,L,Nz,N,r[i],e[1],sigma[1],q[1])
        elif i>3 and i<=9:
            U+=potential(Lz,L,Nz,N,r[i],e[2],sigma[2],q[2])
        elif i>9 and i<=12:
            U+=potential(Lz,L,Nz,N,r[i],e[3],sigma[3],q[3])
    print np.min(U), np.where(U==np.min(U))
    
    #Convert to Hartree
    U=U*(1./320000.) #1Hartree=320000K 
    print np.min(U), np.where(U==np.min(U))
    
    #plot a slice of the potential
    fig = plt.figure()
    #set axes
    x=np.linspace(-3,3,N)#N)
    y=np.linspace(-3,3,N)#N)
    #make X,Y 2D arrays that correspond to eachother
    X,Y=np.meshgrid(x,y)
    #Plot
    ax = fig.gca(projection='3d')
    ax.plot_surface(X,Y,U[19],rstride=1, cstride=1, cmap=cm.gist_rainbow) #Go to 3.5A currently going to 3.6A
    
    
    #For printing along z axis if wanted (must comment out previous 2 lines for this to work)
    #Z=[]
    #for i in range(N):
    #    Z.append(U[i][(N/2)-1][(N/2)-1])
    #print np.min(Z), np.where(Z==np.min(Z))
    ##print Z
    #z=np.linspace(0,6,N)
    #plt.plot(z[:],Z[:])

    #Show
    plt.show()
    
    #Save as .npy file for NuSol.py
    print type(U)
    print np.shape(U)
    name = raw_input('Name the file without an extension:') #naming the file from user input
    data=U
    print ("New File being created: "+name+".npy")
    np.save(name+".npy",data)

    
main()

