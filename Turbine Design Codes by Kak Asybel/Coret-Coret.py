from math import degrees
from pydoc import doc
import numpy as np
import CoolProp
from CoolProp.CoolProp import PropsSI as Pr
from Compute_Rotor import*

NR,r4,Alpha4,b4,Ct4,rho4,mflow=22,0.0291,12.58,3.25064/1000,70.7,32.0081,0.3
aOc = 0.3   #Varopt
Betha2=np.radians(21)   #Varopt
Nn = NR + 1       #Varopt
r3=r4*(1+(2*b4*np.sin(np.radians(Alpha4)))/r4)
Ct3=Ct4*r4/r3
Cm3=mflow/(2*np.pi*r3*b4*rho4)
r2Or3=1.2     #Varopt
r2=r2Or3*r3
s3=2*np.pi*r3/Nn
o3=s3*Cm3/(np.sqrt(Cm3**2+Ct3**2))
Alpha3=np.arcsin(o3/s3)
t2Oc = 0.03   #Varopt
t3Oc = 0.015    #Varopt
tmaxOc = 0.06   #Varopt
dOc = 0.4   #Varopt

bi=0.0003454
ci=0.025497
ai= ci*aOc


#Kalkulasi Y03
X2=np.arctan(4*bi/(4*ai-ci))
X3=np.arctan(4*bi/(3*ci-4*ai))
Tetha=X2+X3
Y2=Betha2+X2
Y03=r2/r3*np.cos(Y2)

Y3=np.arccos(Y03)
c=(r2-r3)/np.sin(Y3)
a=ci*aOc
b=((np.sqrt(1+(4*np.tan(Tetha)**2)*(ai/ci-(ai/ci)**3-3/16)))-1)/(4*np.tan(Tetha))*ci
Errora=(ai-a)
Errorb=(bi-b)
Errorc=(ci-c)

print(ai,bi,ci)
print()
print(a,b,c)
print()
print(Errora,Errorb,Errorc)
print()
