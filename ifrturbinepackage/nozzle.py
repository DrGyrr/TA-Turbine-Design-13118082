import CoolProp
from CoolProp.CoolProp import PropsSI as Props
import numpy as np
from ifrturbinepackage.rotor import *

def VariableNozz(): # INPUT => NR,r4,Alpha4,b4,Ct4,rho4,mflow
  global aOc,Betha2,Nn,r3,Ct3,Cm3,r2Or3,r2,s3,o3,Alpha3,t2Oc,t3Oc,tmaxOc,dOc
  aOc = 0.3               #Varopt
  Betha2=np.radians(21)   #Varopt
  Nn = NR + 1             #Varopt
  r3=r4*(1+(2*b4*np.sin(np.radians(Alpha4)))/r4)
  Ct3=Ct4*r4/r3
  Cm3=mflow/(2*np.pi*r3*b4*rho4)
  r2Or3=1.2       #Varopt
  r2=r2Or3*r3
  s3=2*np.pi*r3/Nn
  o3=s3*Cm3/(np.sqrt(Cm3**2+Ct3**2))
  Alpha3=np.arcsin(o3/s3)
  t2Oc = 0.03     #Varopt
  t3Oc = 0.015    #Varopt
  tmaxOc = 0.06   #Varopt
  dOc = 0.4       #Varopt
# OUTPUT => return(Betha2,Nn,r3,Ct3,Cm3,r2,s3,o3,Alpha3,t2Oc,t3Oc,tmaxOc,dOc)
 
def NozzCalc1(): # INPUT => ai,bi,ci
  #Variabel Desain
  global Y03,Tetha,Y2
  #Kalkulasi Y03
  X2=np.arctan(4*bi/(4*ai-ci))
  X3=np.arctan(4*bi/(3*ci-4*ai))
  Tetha=X2+X3
  Y2=Betha2+X2
  Y03=r2/r3*np.cos(Y2)
# OUTPUT => return (Y03)
       
def NozzCalc2(): # INPUT => ai,bi,ci
  global a,b,c,Errora,Errorb,Errorc
  Y3=np.arccos(Y03)
  c=(r2-r3)/np.sin(Y3)
  a=c*aOc
  b=((np.sqrt(1+(4*np.tan(Tetha)**2)*(ai/ci-(ai/ci)**3-3/16)))-1)/(4*np.tan(Tetha))*ci
  Errora=(ai-a)
  Errorb=(bi-b)
  Errorc=(ci-c)
# OUTPUT => return(a,b,c,Errora,Errorb,Errorc)  

def VariantSearchNozz():
    SListB=[]
    SListC=[]
    SListA=[]
    ErrorList=[]
    # Y03List=[]
    SListSize=1000
    for k in range(0,SListSize):
        Y03=2
        while Y03>1:
            bi=np.abs(np.random.uniform(0.1,0.8)*o3)
            ci=np.abs(np.random.uniform(0,5,1)*(r3-r2))  
            ai=ci*aOc
            Y03=NozzCalc1(ai,bi,ci)[0]
            if Y03<=1:
                break
        ListNozzCalc1=NozzCalc2(ai,bi,ci)
        Errora=ListNozzCalc1[3]
        Errorb=ListNozzCalc1[4]
        Errorc=ListNozzCalc1[5]
        Errorx=np.abs(Errora)+np.abs(Errorb)+np.abs(Errorc)
        SListB.append(bi)
        SListC.append(ci)
        SListA.append(ai)
        ErrorList.append(Errorx)
        # Y03.append(Y03)
        
        # print(SListA)
        # print()
        # print(SListB)
        # print()
        # print(SListC)
        # print()
        # print(ErrorList)
    return(SListA,SListB,SListC,ErrorList)

def CheckValNozz(): # INPUT => anew,bnew,cnew
  X2=np.arctan(4*bnew/(4*anew-cnew))
  X3=np.arctan(4*bnew/(3*cnew-4*anew))
  Tetha=X2+X3
  Y2=Betha2+X2
  Y03=r2/r3*np.cos(Y2)
  Y3=np.arccos(Y03)
  c=(r2-r3)/np.sin(Y3)
  a=c*aOc
  b=((np.sqrt(1+(4*np.tan(Tetha)**2)*(anew/cnew-(anew/cnew)**3-3/16)))-1)/(4*np.tan(Tetha))*cnew
  Errora=(anew-a)
  Errorb=(bnew-b)
  Errorc=(cnew-c)
  Errorxnoz=abs(Errora)+abs(Errorb)+abs(Errorc)
# OUTPUT => return (anew,bnew,cnew,Errorxnoz,X2,X3,Y2,Y3)


#Menggambar 3D Nozzle
def nozzleconts(): # INPUT => a,b,c,t2Oc,t3Oc,dOc
  t2=c*t2Oc
  t3=c*t3Oc
  tmax=c*tmaxOc
  d=c*dOc

  #Camberline Design
  NozzGrid=50
  Xc=[]
  Yc=[]
  Xi=[]
  Xc.append(0)
  Yc.append(0)
  Xi.append(0)
  tref=[]
  epsi=[]
  enozz=[]
  tnozz=[]
  Xa=[]
  Ya=[]
  Xb=[]
  Yb=[]
  Xa.append(0)
  Ya.append(0)
  Xb.append(0)
  Yb.append(0)
  for i in range(1,NozzGrid):
    Xc.append(Xc[i-1]+c/NozzGrid)
    Yc.append(Xc[i]*(c-Xc[i])/(((c-2*a)**2)/(4*b**2)+(c-2*a)/b*Xc[i]-(c**2-4*a*c)/(4*b)))
    Xi.append(np.arctan(Yc[i]-Yc[i-1])/(Xc[i]-Xc[i-1]))

#T Calculation
  for i in range(1,len(Xc)):
    tref.append(t2+(t3-t2)*(Xc[i]/d))
    if Xc[i]<=d:
      epsi.append(Xc[i]/d)
    else:
      epsi.append((c-Xc[i])/(c-d))

  for i in range(1,len(epsi)):
    enozz.append(np.sqrt((0.4*d/c)*(0.95*(1-Xc[i]/c)*(1-epsi[i]+0.05))))
  for i in range(1,len(enozz)):
    tnozz.append(tref[i]+(tmax-tref)*(epsi[i]**enozz[i]))

#Coordinate Calculation
  for i in range(1,len(tnozz)):
    Xa.append(Xc[i]+0.5*tnozz[i]*np.sin(Xi[i]))
    Ya.append(Yc[i]+0.5*tnozz[i]*np.cos(Xi[i]))
    Xb.append(Xc[i]-0.5*tnozz[i]*np.sin(Xi[i]))
    Yb.append(Yc[i]-0.5*tnozz[i]*np.cos(Xi[i]))

  print(Xc)
  print()
  print(Xa)
#OUTPUT => return(Xa,Xb,Ya,Yb)

def proceedN(savenumber):
    VariableNozz()
    # so on and so on. //// still under construction
    
