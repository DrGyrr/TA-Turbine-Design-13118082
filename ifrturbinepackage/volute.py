import CoolProp
from CoolProp.CoolProp import PropsSI as Props
import numpy as np
from ifrturbinepackage.rotor import *
from ifrturbinepackage.nozzle import *

def ComputeV(): # INPUT => fluid,r2,r3,Cm3,Betha2,H_1,T_1,P_1,mflow
    Cm2=r2/r3*Cm3
    Ct2=Cm2/np.tan(np.radians(Betha2))
    AR=1
    Vaneless=0.002
    vB=0.001
    vA=AR*vB
    r1=r2+vA+Vaneless
    C1=Ct2*r2/r1
    h1=H_1-1/2*C1**2
    rmax=r1+vB
    A1=(3/4*np.pi+1)*vA*vB
    rho1=Props('D','T',T_1,'P',P_1*1e6,fluid)
    A1i=mflow/(rho1*C1)
    vBi=A1i/(3/4*np.pi+1)/vA
    ErrorA1=abs(A1-A1i)
    # print(A1N,A1iN,C1,ErrorA1,ErrorPercent)
    while ErrorA1 > 10e-7:
        vB=vB+(vBi-vB)/2
        vA=AR*vB
        r1=r2+vA+Vaneless
        C1=Ct2*r2/r1
        h1=H_1-1/2*C1**2
        rmax=r1+vB
        A1=(3/4*np.pi+1)*vA*vB
        A1i=mflow/rho1*C1
        vBi=A1i/(3/4*np.pi+1)/vA
        ErrorA1=np.abs(A1-A1i)
        #print(ErrorA1/A1*100)
        if ErrorA1 <=10e-7:
            print(ErrorA1/A1*100)
            break
# OUTPUT => return(vA,vB,r1)