#Creating The MOF-5 Structure
from paraview.simple import *
import math
#6 atoms by 48 symmetries (adding 24 atoms to finish)
M=144+24#H
N=192+48#H
#5 bonds by 24 symmetries plus 48 bonds to finish benzene rings and 48 bonds for final oxygens and 24 H bonds
K=120+24#Hbonds
L=216+48#Hbonds
#ADDED for Version 5.0
def Rotation(preRotation):
    postRotation=[0,0,0]
    #Rotation Matrix
    RotMatrix=[[-1/(math.sqrt(6)),2/(math.sqrt(6)),-1/(math.sqrt(6))],[1/(math.sqrt(2)),0,-1/(math.sqrt(2))],[1/(math.sqrt(3)),1/(math.sqrt(3)),1/(math.sqrt(3))]]
    #Identity Matrix
    #RotMatrix=[[1,0,0],[0,1,0],[0,0,1]]
    for i in range(3):
        for j in range(3):
            postRotation[i]+=RotMatrix[i][j]*preRotation[j]
    return postRotation

#Symmetry function
def sym(x,y,z,which):
    RETURNS=[]
    RETURNS=[[x,y,z],[-x,-y,z],[-x,y,-z],[x,-y,-z],[z,x,y], #1-5
    [z,-x,-y],[-z,-x,y],[-z,x,-y],[y,z,x],[-y,z,-x],        #6-10
    [y,-z,-x],[-y,-z,x],[y,x,-z],[-y,-x,-z],[y,-x,z],       #11-15      Base Symmetries that produce most of the unit cell without overlap
    [-y,x,z],[x,z,-y],[-x,z,y],[-x,-z,-y],[x,-z,y],         #16-20
    [z,y,-x],[z,-y,x],[-z,y,x],[-z,-y,-x],                  #21-24
#The symmetries in addition to crystal maker's explicit symmetries
    [z,y,x],[-z,x,y],[y,x,z],[x,z,y],[x,-y,z],              #1-5
    [y,-x,-z],[x,-z,-y],[y,z,-x],[-z,y,-x],[z,x,-y],        #6-10
    [y,-z,x],[-z,-y,x],[-x,z,-y],[-y,-z,-x],[-y,-x,z],      #11-15      Symmetries that finish the unit cell
    [-x,-y,-z],[-x,-z,y],[-y,z,x],[z,-x,y],[x,y,-z],        #16-20
    [-y,x,-z],[-x,y,z],[z,-y,-x],[-z,-x,-y]]                #21-24
    
    preRotation=RETURNS[which]
    postRotation=Rotation(preRotation)
    return postRotation[0],postRotation[1],postRotation[2]
    
    
#create empty list of bonds
Bonds=[]
for i in range (L):
    Bonds.append(Line())
#Define Display Properties
BondsDP=[]
for i in range(L):
    BondsDP.append(GetDisplayProperties(Bonds[i]))
    
#create empty list of atoms
unit = []
#make spheres radius LARGER TODO .015 .02
for i in range(N):
    unit.append(Sphere())

        
#Define Display Properties
unitDP=[]
for i in range(N):
    unitDP.append(GetDisplayProperties(unit[i]))
     
#Array of Coordinates (6 atoms by 3 coordinates O1-Zn-O2-C1-C2-C3-H)
unitCoord = [[6.4775, 6.4775, 6.4775], [7.6136, 5.3413, 5.3413], [7.2721, 5.6830, 3.4548], [6.4775, 6.4775, 2.8975], [6.4775, 6.4775, 1.3820], [7.3224, 5.6325, 0.6905], [7.9981, 4.9568, 1.2522]]
#U=number of atoms in basic group
U=7

#Application of Translation vector
#Tx, Ty and Tz denote a translation given as the components of a vector
Tx=0.
Ty=0.
Tz=0.

#Radius of all atoms-visulatization of atomic radius is not to scale
RADIUS=.3
#Alter parameters for each atom
for j in range (24): #Do the main picture with the first 24 symmetries
    for i in range(U):
        x=unitCoord[i][0]
        y=unitCoord[i][1]
        z=unitCoord[i][2]
        # Call one symmetry method
        x,y,z=sym(x,y,z,j)
        unit[i+U*j].Radius = RADIUS
        unit[i+U*j].PhiResolution = 15
        unit[i+U*j].ThetaResolution = 15
        unit[i+U*j].Center = (x+Tx,y+Ty,z+Tz)
    #make Zn light blue
        if (i==1):    
            unitDP[i+U*j].DiffuseColor=[.6,.8,.9]
    #make O1,O2 Red
        elif (i==0 or i==2):
            unitDP[i+U*j].DiffuseColor=[.8,.1,.1]
    #make C1,C2,C3 Black
        elif (i>2 and i<6):
            unitDP[i+U*j].DiffuseColor=[0.3,0.3,0.3]
    #don't alter H color
    for k in range(U-1):
        #make start point for bond
        x1=unitCoord[k][0]
        y1=unitCoord[k][1]
        z1=unitCoord[k][2]
        #make endpoint for bond
        x2=unitCoord[k+1][0]
        y2=unitCoord[k+1][1]
        z2=unitCoord[k+1][2]
        #Call same symmetry for both points
        x1,y1,z1=sym(x1,y1,z1,j)
        x2,y2,z2=sym(x2,y2,z2,j)
        Bonds[k+(U-1)*j].Point1=(x1+Tx,y1+Ty,z1+Tz)
        Bonds[k+(U-1)*j].Point2=(x2+Tx,y2+Ty,z2+Tz)
        BondsDP[k+(U-1)*j].DiffuseColor=[0.,0.,0.] #MAKE IT BLACK!

#Finish Benzene Rings and Final Oxygens
for a in range(24):
    #Oxygens (Spheres 192)
    xOxygen=unitCoord[2][0]
    yOxygen=unitCoord[2][1]
    zOxygen=unitCoord[2][2]
    # Call one symmetry method
    xOxygen,yOxygen,zOxygen=sym(xOxygen,yOxygen,zOxygen,a+24)
    unit[M+a].Radius = RADIUS
    unit[M+a].PhiResolution = 15
    unit[M+a].ThetaResolution = 15
    unit[M+a].Center = (xOxygen+Tx,yOxygen+Ty,zOxygen+Tz)
    unitDP[M+a].DiffuseColor=[.8,.1,.1]
    
    #OXYGEN BONDS
    #Bond to Oxygen 2-3
    x1=unitCoord[1][0]
    y1=unitCoord[1][1]
    z1=unitCoord[1][2]
    #make endpoint for bond
    x2=unitCoord[2][0]
    y2=unitCoord[2][1]
    z2=unitCoord[2][2]
    #Call same symmetry for both points
    x1,y1,z1=sym(x1,y1,z1,a+24)
    x2,y2,z2=sym(x2,y2,z2,a+24)
    Bonds[K+a].Point1=(x1+Tx,y1+Ty,z1+Tz)
    Bonds[K+a].Point2=(x2+Tx,y2+Ty,z2+Tz)
    BondsDP[K+a].DiffuseColor=[0.,0.,0.] #MAKE IT BLACK!
    #Bond from Oxygen 3-4
    x1=unitCoord[2][0]
    y1=unitCoord[2][1]
    z1=unitCoord[2][2]
    #make endpoint for bond
    x2=unitCoord[3][0]
    y2=unitCoord[3][1]
    z2=unitCoord[3][2]
    #Call same symmetry for both points
    x1,y1,z1=sym(x1,y1,z1,a+24)
    x2,y2,z2=sym(x2,y2,z2,a+24)
    #First 
    Bonds[K+a+24].Point1=(x1+Tx,y1+Ty,z1+Tz) 
    Bonds[K+a+24].Point2=(x2+Tx,y2+Ty,z2+Tz)
    BondsDP[K+a+24].DiffuseColor=[0.,0.,0.] #MAKE IT BLACK!
    
    #Carbons
    xCarbon=unitCoord[5][0]
    yCarbon=unitCoord[5][1]
    zCarbon=unitCoord[5][2]
    # Call one symmetry method
    xCarbon,yCarbon,zCarbon=sym(xCarbon,yCarbon,zCarbon,a+24)
    unit[M+a+24].Radius = RADIUS
    unit[M+a+24].PhiResolution = 15
    unit[M+a+24].ThetaResolution = 15
    unit[M+a+24].Center = (xCarbon+Tx,yCarbon+Ty,zCarbon+Tz)    
    unitDP[M+a+24].DiffuseColor=[0.3,0.3,0.3]
    
    #Bond to Carbon Oxygen 5-6
    x1=unitCoord[4][0]
    y1=unitCoord[4][1]
    z1=unitCoord[4][2]
    #make endpoint for bond
    x2=unitCoord[5][0]
    y2=unitCoord[5][1]
    z2=unitCoord[5][2]
    #Call same symmetry for both points
    x1,y1,z1=sym(x1,y1,z1,a+24)
    x2,y2,z2=sym(x2,y2,z2,a+24)
    Bonds[K+a+48].Point1=(x1+Tx,y1+Ty,z1+Tz)
    Bonds[K+a+48].Point2=(x2+Tx,y2+Ty,z2+Tz)
    BondsDP[K+a+48].DiffuseColor=[0.,0.,0.] #MAKE IT BLACK!
    
    #Final H
    Hx=unitCoord[6][0]
    Hy=unitCoord[6][1]
    Hz=unitCoord[6][2]
    #Call one symmetry method
    Hx,Hy,Hz=sym(Hx,Hy,Hz,a+24)
    #print(M+a+48)
    unit[M+a+48].Radius = RADIUS
    unit[M+a+48].PhiResolution = 15
    unit[M+a+48].ThetaResolution = 15
    unit[M+a+48].Center = (Hx+Tx,Hy+Ty,Hz+Tz)
    
    #Bond to Carbon Hydrogen 6-7
    x1=unitCoord[5][0]
    y1=unitCoord[5][1]
    z1=unitCoord[5][2]
    #make endpoint for bond
    x2=unitCoord[6][0]
    y2=unitCoord[6][1]
    z2=unitCoord[6][2]
    #Call same symmetry for both points
    x1,y1,z1=sym(x1,y1,z1,a+24)
    x2,y2,z2=sym(x2,y2,z2,a+24)
    Bonds[K+a+72].Point1=(x1+Tx,y1+Ty,z1+Tz)#3*24=72
    Bonds[K+a+72].Point2=(x2+Tx,y2+Ty,z2+Tz)
    BondsDP[K+a+72].DiffuseColor=[0.,0.,0.] #MAKE IT BLACK!
    
    
#Final Benzene Ring Bonds
for j in range (24):
    x1=unitCoord[5][0]
    y1=unitCoord[5][1]
    z1=unitCoord[5][2]
    #make endpoint for bond
    x2=unitCoord[5][0]
    y2=unitCoord[5][1]
    z2=-unitCoord[5][2]
    #Call same symmetry for both points
    x1,y1,z1=sym(x1,y1,z1,j)
    x2,y2,z2=sym(x2,y2,z2,j)
    Bonds[L-j-1].Point1=(x1+Tx,y1+Ty,z1+Tz)
    Bonds[L-j-1].Point2=(x2+Tx,y2+Ty,z2+Tz)
    BondsDP[L-j-1].DiffuseColor=[0.,0.,0.] #MAKE IT BLACK!

#Checks if coordinate is within range we are considering
#center = middle of circle, accept = radius of circle, coord = considered coordinate
def inRangeSphere(center,accept,coord):
    if center=="None":
        return True
    else:
        #if distance is witin acceptable range then return True
        distance=magnitude(coord,center)
        if distance<=accept:
            return True
        else:
            return False
    
def magnitude(r,rspace):
    #3D distance formula
    return math.sqrt((r[0]-rspace[0])**2+(r[1]-rspace[1])**2+(r[2]-rspace[2])**2)

#Decide a portion of the structure to show
#If you want to show the whole thing put center="None"
center=[0,0,-11.2193591060274] #Currently set to have center at the lowest atom on the z-axis which is a center oxygen in and ion cluster
#Acceptance radius
accept=3.61
#Show Spheres
for i in range(N):
    s=unit[i]
    if inRangeSphere(center,accept,s.Center):
        Show(s)
    else:
        Hide(s)
#Show Bonds
for i in range (L):
    Line=Bonds[i]
    if inRangeSphere(center,accept,Line.Point1) and inRangeSphere(center,accept,Line.Point2):
        Show(Line)
    else:
        Hide(Line)



Render()


