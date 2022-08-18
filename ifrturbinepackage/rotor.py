import numpy as np
import CoolProp
from CoolProp.CoolProp import PropsSI as Props
import pandas as pd
import matplotlib.pyplot as plt
from ifrturbinepackage.definitions import *
from ifrturbinepackage.inputs import *
import scipy
from scipy import optimize

#ComputeR1 is deprecated
def ComputeR1(tenflow_coeff,tenwork_coeff,k):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4ss,T05ss,T05,T5ss,T5
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1C5,k2C5,errorC5
    global TotalLoss,LossInc,LossPass,LossTip,LossWind,LossTE,LossExit,rho4m,S5,O5
    global Effts,Efftt,Efftspred,Reaction,vNu

    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10

    whichcycle (k)       # The cycle to be computed
    
    Cp4 = Props('C','T',T_1,'P',P_1,fluid)
    Cv4 = Props('O','T',T_1,'P',P_1,fluid)
    gamma = Cp4/Cv4
    Rx = 8.31446261815324   #J/K.mol

    #General Properties inlet outlet turbin (Total)
    H_1     = Props('H','T',T_1,'P',P_1,fluid)     #J/kg
    s01     = Props('S','T',T_1,'P',P_1,fluid)     #J/kg.K 
    T_5     = Props('T','P',P_5,'S',s01,fluid)  # =>asumsi nozzle isenthalpy DAN Isentropic
    H_5     = Props('H','T',T_5,'P',P_5,fluid)  # meski pada kenyataannya isenthalpic nozzle tidak isentropic
    DeltaH  = H_1-H_5            #Ideal === Isentropic Total Enthalpy change 

    C0s     = np.sqrt(2*DeltaH)         #Spouting Velocity

    #Segitiga Kecepatan Inlet, m/s, radians
    U4      = np.sqrt(DeltaH/work_coeff)
    Cm4     = U4*flow_coeff
    Ct4     = DeltaH/U4                 # => DeltaH = U4*Ct4-U5*Ct5 ; Alpha5=0 => Ct5=0
    C4      = np.sqrt(Cm4**2+Ct4**2)
    Alpha4  = np.arctan(Ct4/Cm4)
    W4      = np.sqrt(Cm4**2+(U4-Ct4)**2)
    Beta4   = np.arctan((U4-Ct4)/Cm4)

    #Perhitungan Properties ideal lain (Total)
    p01     = P_1           #inlet volute [1], Total
    T01     = T_1
    h01     = H_1
    p1      = p01
    T1      = T_1
    h01     = H_1
    rho1   = Props('D','P',p1,'T',T1,fluid)
    
    p05ss   = P_5
    T05ss   = T_5
    h05ss   = H_5
    s05ss   = Props('S','H',h05ss,'P',p05ss,fluid)
    rho05ss = Props('D','P',p05ss,'T',T05ss,fluid)

    Q5      = mflow/rho05ss
    ns      = np.radians(rpm*6)*np.sqrt(Q5)/DeltaH**0.75
    Efftspred    = 0.81-1.07*(ns-0.55)**2-0.5*(ns-0.55)**3       #predicted total-to-static efficiency

    
    h02s    = H_1           #inlet nozzle [2], Total
    s02s    = s01            #ideal volute === approx. as isentropic
    p02s    = p01
    T02s    = T01
    h03s    = h02s           #outlet nozzle [3], Total
    s03s    = s02s            #ideal nozzle === approx. as isentropic (in Total)
    p03s    = p02s
    T03s    = T02s
    h04s    = h03s           #inlet rotor [4], Total
    s04s    = s03s           #outlet nozzle === inlet rotor
    p04s    = p03s
    T04s    = T03s
    h04     = h04s          # Nozzle isenthalpic but not isentropic
    p04     = p01-rho1*DeltaH*(1-Efftspred)/4
    T04     = Props('T','P',p04,'H',h04s,fluid)
    s04     = Props('S','P',p04,'T',T04,fluid)
    
    h05ss   = H_5
    p05ss   = P_5
    T05ss   = T_5
    s05ss   = Props('S','H',h05ss,'P',p05ss,fluid)

    #Perhitungan Properties ideal lain (Static)
    h4s     = h04s-1/2*C4**2
    p4s     = Props('P','H',h4s,'S',s04s,fluid)
    p4      = Props('P','H',h4s,'S',s04,fluid)
    T4s     = Props('T','H',h4s,'S',s04s,fluid)
    rho04s  = Props('D','P',p04s,'T',T04s,fluid)
    rho4s   = Props('D','P',p4s,'T',T4s,fluid)
    rho4sm  = 2*(p04s-p4s)/C4**2
    h4      = h04-1/2*C4**2
    p4      = Props('P','H',h4,'S',s04,fluid)
    T4      = Props('T','H',h4,'S',s04,fluid)
    rho04   = Props('D','P',p04,'T',T04,fluid)
    rho4    = Props('D','P',p4,'H',h04,fluid)
    rho4m   = 2*(p04-p4)/C4**2
    a01     = Props('A','P',p01,'T',T01,fluid)
    a4s     = Props('A','P',p4s,'T',T4s,fluid)
    a4      = Props('A','P',p4,'T',T4,fluid)
    Ma4s    = C4/a4s
    Ma4     = C4/a4



    #Perhitungan Geometri 
    r4  = U4/np.radians(rpm*6)
    b4  = mflow/(2*np.pi*r4*rho4s*Cm4)
    D4  = r4*2
    Zr  = Zratio*r4
    rh5 = Rrh5r4*r4
    b5  = Rb5b4*b4  
    rs5 = rh5+b5
    r5  = (rs5+rh5)/2

    Re4s    = rho4s*C4*b4/Props('V','P',p4s,'T',T4s,fluid)
    Re4     = rho4*C4*b4/Props('V','P',p4,'T',T4,fluid)

    #Segitiga Kecepatan Outlet
    U5      = r5*np.radians(rpm*6)
    Cm5_0    = 0
    rho5ss_0= rho05ss
    Cm5ii    = Cm5_0
    rho5ssii= rho5ss_0
    Cm5didconverge1 = False
    Cm5didconverge2 = False
    k1Cm5    = 0
    k2Cm5    = 0
    while Cm5didconverge1 == False:
        k1Cm5     = k1Cm5+1             # => iteration amount
        Cm5i      = Cm5ii
        rho5ssi = rho5ssii
        Cm5ii    = mflow/(rho5ssi*2*np.pi*r5*b5)
        h5ss    = h05ss-1/2*Cm5ii**2    # => it is predetermined that Ct5=0 since Alpha5===0
        rho5ssii  = Props('D','H',h5ss,'S',s05ss,fluid)
        errorCm5= mflow/(rho5ssii*Cm5ii*2*np.pi*b5*r5)-1
        if np.abs(errorCm5) <= 10**-10:
            Cm5didconverge1 = True
            Cm5didconverge2 = True
            Cm5  = Cm5ii
            rho5ss=rho5ssii
            break
        if (rho5ssi*Cm5i-rho5ssii*Cm5ii)*(Cm5i-Cm5ii)<0:
            Cm5      = Cm5ii
            rho5ss  = rho5ssii
            break
    while Cm5didconverge2 == False:
        k2Cm5    = k2Cm5 +1         # => iteration amount
        Cm5      = mflow/(rho5ss*2*np.pi*r5*b5)
        h5ss    = h05ss-1/2*Cm5**2
        rho5ss  = Props('D','H',h5ss,'S',s05ss,fluid)
        errorCm5= mflow/(rho5ss*Cm5*2*np.pi*b5*r5)-1
        if errorCm5 <= 10**-10:
            Cm5didconverge2 = True
            break
    h5ss    = h05ss-1/2*Cm5**2
    Alpha5  = 0
    Ct5     = 0
    C5      = np.sqrt(Ct5**2+Cm5**2)
    Beta5   = np.arctan((U5-Ct5)/Cm5)
    W5      = Cm5/np.cos(Beta5)

    S5      = 2*np.pi*r5/NR
    O5      = S5*Cm5/W5

    p5ss    = Props('P','H',h5ss,'S',s05ss,fluid)
    a5ss    = Props('A','H',h5ss,'P',p5ss,fluid)
    Ma5ss   = C5/a5ss


    # \\\\\\\ <<---------<<----||----->>------------>> ////////
    ## Losses Coefficient ##

    #Rotor Incidence Losses 
    Beta4opt2 = np.arctan((-1.98/NR)/(1-1.98/NR)*np.tan(Alpha4))
    Beta4opt= np.arctan(np.tan(Alpha4)*(work_coeff-1+2/NR)/work_coeff)  #(Chen)
    LossInc0 = 0.5*(W4**2)*(np.sin(np.abs(np.abs(Beta4)-np.abs(Beta4opt))))**2  #m2/s2
    LossInc     = 0.5*(W4**2)*(np.sin(Beta4)-np.sin(Beta4opt))**2
       
    #Blade loading efficiency (Chen)
    vNu = U4/np.sqrt(2*Cp4*T01*(1-(p5ss/p01)**((gamma-1)/gamma))) #blade/isentropic jet speed ratio
    Effreductbladeloading = flow_coeff**2*vNu**2
    #Rotor Passage Losses ([Uusitalo] from Moustapha PLM3)
    LH = np.pi/4*((Zr-b4/2)+(r4-rh5-b5/2))                                                              #m
    DH = 0.5*((4*np.pi*r4*b4/(2*np.pi*r4+Zr*rh5))+((2*np.pi*(rs5**2-rh5**2)/(np.pi*(rs5-rh5))+Zr*b5)))  #m
    Y5 = np.arctan(0.5*(np.tan(Beta4)+np.tan(Beta5)))
    C = Zr/np.cos(Y5)
    if (r4-rs5)/b5>=0.2:
        KpCETI = 0.11
    else:
        KpCETI = 0.22
    LossPass = KpCETI*(LH/DH+0.68*((1-(r5/r4)**2)*np.cos(Beta5)/(b5/C))*((W4**2+W5**2)/2))
    
    #Rotor Clearance Losses
    Ca = (1-(rs5/r4))/(Cm4*b4)
    Cr = (rs5/r4)*((Zr-b4)/(Cm5*r5*b5))
    Ka = 0.4
    Kr = 0.75
    Kar = -0.3
    Ea = 0.0003
    Er = 0.0003
    if Ea*Er*Ca*Cr>=0:
        LossTip = (U4**3*NR/(8*np.pi))*(Ka*Ea*Ca+Kr*Er*Cr+Kar*np.sqrt(Ea*Er*Ca*Cr))
    else:
        LossTip = (U4**3*NR/(8*np.pi))*(Ka*Ea*Ca+Kr*Er*Cr)
    #Windage Losses # disk friction losses (fiaschi, 2015: 4.2.5)
    Eb = 0.0003
    Kf = 3.7*(Eb/r4)**0.1/Re4s**0.5
    LossWind = Kf*((rho4s+rho5ss)/2)*U4**3*r4**2/(2*mflow*W5**2)

    #Trailing Edge Losses
    tb4 = 0.04*r4
    tb5 = 0.02*r4
    LossTE = rho5ss*W5**2/2*(NR*tb5/(np.pi*(rh5+rs5)*np.cos(Beta5)))**2
    
    #Exit Losses
    LossExit = 0.5*C5**2    # => mungkin untuk ubah dari total jadi static. abaikan dulu

    # => Sum Enthalpy Losses
    TotalLoss = LossInc + LossPass + LossTip + LossWind + LossTE
    # \\\\\\\ <<---------<<----||----->>------------>> ////////


    #Perhitungan Properties considering losses
    h05 = h05ss+ (LossInc0+LossPass+LossTip+LossWind+LossTE)             #nozzle masih diasumsikan isentropic dan isenthalpic
    h5  = h5ss+ (LossInc+LossPass+LossTip+LossWind+LossTE  )   
    p05 = p05ss
    T05 = Props('T','H',h05,'P',p05,fluid)
    p5  = p5ss
    T5  = Props('T','H',h5,'P',p5,fluid)

    #Effisiensi 
    Reaction    = (h4s-h5)/(h01-h05ss)
    Efftt       = ((h01-h05)/(h01-h05ss)-Effreductbladeloading)*100
    Effts       = ((h01-h05)/(h01-h5ss)-Effreductbladeloading)*100

def ComputeR2(tenflow_coeff,tenwork_coeff,k,l,m):

    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4ss,T05ss,T05,T5ss,T5
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm4didconverge1,Cm4didconverge2,k1Cm4,k2Cm4,errorCm4
    global TotalLoss,LossInc,LossPass,LossTip,LossWind,LossTE,LossExit,rho4m,S5,O5
    global Effts,Efftt,Efftspred,Reaction,choked4

    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10

    cycledict=whichcycle (k)       # The cycle to be computed
    locals().update(cycledict)
    gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    locals().update(gparamdict)        # rpm at m will be used
    rpm =whatrpm(m)          # rpm at m will be used

    P_1 = P_1*10**6
    P_5 = P_5*10**6
    Cp4 = Props('C','T',T_1,'P',P_1,fluid)
    Cv4 = Props('O','T',T_1,'P',P_1,fluid)
    gamma = Cp4/Cv4
    Rx = 8.31446261815324   #J/K.mol

    #General Properties inlet outlet turbin (Total)
    H_1     = Props('H','T',T_1,'P',P_1,fluid)     #J/kg
    s01     = Props('S','T',T_1,'P',P_1,fluid)     #J/kg.K 
    T_5     = Props('T','P',P_5,'S',s01,fluid)  # =>asumsi nozzle isenthalpy DAN Isentropic
    H_5     = Props('H','T',T_5,'P',P_5,fluid)  # meski pada kenyataannya isenthalpic nozzle tidak isentropic
    DeltaH  = H_1-H_5            #Ideal === Isentropic Total Enthalpy change 

    C0s     = np.sqrt(2*DeltaH)         #Spouting Velocity

    #Perhitungan Properties ideal lain (Total)
    p01     = P_1           #inlet volute [1], Total
    T01     = T_1
    h01     = H_1
    p1      = p01           # inlet turbine, V~0 
    T1      = T_1
    h01     = H_1
    rho1   = Props('D','P',p1,'T',T1,fluid)
    h02s    = H_1           #inlet nozzle [2], Total
    s02s    = s01            #ideal volute === approx. as isentropic
    p02s    = p01
    T02s    = T01
    h03s    = h02s           #outlet nozzle [3], Total
    s03s    = s02s            #ideal nozzle === approx. as isentropic (in Total)
    p03s    = p02s
    T03s    = T02s
    h04s    = h03s           #inlet rotor [4], Total
    s04s    = s03s           #outlet nozzle === inlet rotor
    p04s    = p03s
    T04s    = T03s
    h04     = h04s          # Nozzle isenthalpic but not isentropic
    p05ss   = P_5
    T05ss   = T_5
    h05ss   = H_5
    s05ss   = s04s

    #Segitiga Kecepatan , m/s, radians
    U4      = np.sqrt(DeltaH/work_coeff)
    Ct4     = DeltaH/U4
    Cm5     = flow_coeff*U4
    h5ss    = h05ss-1/2*Cm5**2 # => Ct5=0, jadi C=Cm
    rho5ss  = Props('D','H',h5ss,'S',s01,fluid)
    rho05ss = Props('D','H',h05ss,'S',s01,fluid)
    # mflow   = rho5ss*2*np.pi*b5*r5*Cm5
    rho04s  = Props('D','H',h04,'S',s01,fluid)

    Cm4_0    = 0
    rho4s_0= rho04s        # => initial value for iteration
    Cm4ii    = Cm4_0
    rho4sii= rho4s_0
    Cm4didconverge1 = False
    Cm4didconverge2 = False
    choked4     = False
    k1Cm4    = 0
    k2Cm4    = 0
    while Cm4didconverge1 == False:
        k1Cm4       = k1Cm4+1             # => iteration amount
        Cm4i        = Cm4ii
        rho4si      = rho4sii
        Cm4ii       = Rb5b4*Rr5r4*(rho5ss/rho4si)*Cm5
        h4s         = h04s-1/2*(Cm4ii**2+Ct4**2)
        rho4sii     = Props('D','H',h4s,'S',s04s,fluid)
        errorCm4    = np.abs(Rb5b4*Rr5r4*(rho5ss/rho4sii)*(Cm5/Cm4ii)-1)
        # errorCm4    = mflow/(rho5ssii*Cm4ii*2*np.pi*b4*r4)-1
        if errorCm4 <= 10**-10:
            Cm4didconverge1 = True
            Cm4didconverge2 = True
            Cm4     = Cm4ii
            rho4s   = rho4sii
            break
        if (rho4si*Cm4i-rho4sii*Cm4ii)*(Cm4i-Cm4ii)<0:
            Cm4      = Cm4ii
            rho4s  = rho4sii
            break
    while Cm4didconverge2 == False:
        k2Cm4     = k2Cm4 +1         # => iteration amount
        Cm4       = Rb5b4*Rr5r4*(rho5ss/rho4s)*Cm5
        h4s       = h04ss-1/2*(Cm4ii**2+Ct4**2)
        rho4s     = Props('D','H',h4s,'S',s04s,fluid)
        if np.abs(1-Cm4/Props('A','H',h4s,'S',s04s,fluid)) < 5*1e-3:
            choked4 = True
            break
        errorCm4  = np.abs(Rb5b4*Rr5r4*(rho5ss/rho4s)*(Cm5/Cm4)-1)
        if errorCm4 <= 10**-10:
            Cm4didconverge2 = True
            break
    h4s     = h04s-1/2*(Cm4**2+Ct4**2)
    C4      = np.sqrt(Cm4**2+Ct4**2)
    W4      = np.sqrt(Cm4**2+(U4-Ct4)**2)
    Alpha4  = np.arccos(Cm4/C4)
    Beta4   = np.arccos(Cm4/W4)
    Alpha5  = 0
    Ct5     = Cm5*np.tan(Alpha5)
    U5      = U4*Rr5r4
    W5      = np.sqrt(Cm5**2+(U5-Ct5)**2)
    C5      = np.sqrt(Ct5**2+Cm5**2)
    Beta5   = np.arccos(Cm5/W5)

    #Perhitungan geometri
    r4      = U4/np.radians(rpm*6)
    r5      = Rr5r4*r4
    b4      = Rb4r4*r4
    b5      = Rb5b4*b4
    rs5     = (2*r5+b5)/2
    rh5     = rs5-b5
    if rh5 < 0.0015:
        print("For flow coeff =",flow_coeff,"and work coeff=",work_coeff,"rh5 too small(<1.5mm, adjust gparams")
        return
    Zr      = RZrr4*r4

    mflow   = 2*np.pi()*b5*r5*rho5ss*Cm5

    Q5      = mflow/rho05ss
    ns      = np.radians(rpm*6)*np.sqrt(Q5)/DeltaH**0.75
    Efftspred    = 0.81-1.07*(ns-0.55)**2-0.5*(ns-0.55)**3       #predicted total-to-static efficiency

    

    p04     = p01-rho1*DeltaH*(1-Efftspred)/4   #predicting loss due to nozzle
    T04     = Props('T','P',p04,'H',h04s,fluid)
    s04     = Props('S','P',p04,'T',T04,fluid)
    

    #Perhitungan Properties ideal lain (Static)
    h4s     = h04s-1/2*C4**2
    p4s     = Props('P','H',h4s,'S',s04s,fluid)
    p4      = Props('P','H',h4s,'S',s04,fluid)
    T4s     = Props('T','H',h4s,'S',s04s,fluid)
    rho04s  = Props('D','P',p04s,'T',T04s,fluid)
    rho4s   = Props('D','P',p4s,'T',T4s,fluid)
    rho4sm  = 2*(p04s-p4s)/C4**2
    h4      = h04-1/2*C4**2
    p4      = Props('P','H',h4,'S',s04,fluid)
    T4      = Props('T','H',h4,'S',s04,fluid)
    rho04   = Props('D','P',p04,'T',T04,fluid)
    rho4    = Props('D','P',p4,'H',h04,fluid)
    rho4m   = 2*(p04-p4)/C4**2
    a01     = Props('A','P',p01,'T',T01,fluid)
    a4s     = Props('A','P',p4s,'T',T4s,fluid)
    a4      = Props('A','P',p4,'T',T4,fluid)
    Ma4s    = C4/a4s
    Ma4     = C4/a4


    Re4s    = rho4s*Cm4*b4/Props('V','P',p4s,'T',T4s,fluid)
    # Re4     = rho4*C4*b4/Props('V','P',p4,'T',T4,fluid)



    S5      = 2*np.pi*r5/NR
    O5      = S5*Cm5/W5

    p5ss    = Props('P','H',h5ss,'S',s05ss,fluid)
    a5ss    = Props('A','H',h5ss,'P',p5ss,fluid)
    Ma5ss   = C5/a5ss


    # \\\\\\\ <<---------<<----||----->>------------>> ////////
    ## Losses Coefficient ##

    #Rotor Incidence Losses (Chen)
    Beta4opt2 = np.arctan((-1.98/NR)/(1-1.98/NR)*np.tan(Alpha4)) #
    Beta4opt  = np.arctan(np.tan(Alpha4)*(work_coeff-1+2/NR)/work_coeff)    #(Chen)
    LossInc0  = 0.5*(W4**2)*(np.sin(np.abs(np.abs(Beta4)-np.abs(Beta4opt))))**2  #m2/s2
    LossInc   = 0.5*(W4**2)*(np.sin(Beta4)-np.sin(Beta4opt))**2

    #Blade loading efficiency (Chen) 
    vNu = U4/np.sqrt(2*Cp4*T01*(1-(p5ss/p01)**((gamma-1)/gamma))) #blade/isentropic jet speed ratio
    Effreductbladeloading = flow_coeff**2*vNu**2

    #Rotor Passage Losses([Uusitalo] from Moustapha. PLM3)
    LH = np.pi/4*((Zr-b4/2)+(r4-rh5-b5/2))                                                              #m
    DH = 0.5*((4*np.pi*r4*b4/(2*np.pi*r4+Zr*rh5))+((2*np.pi*(rs5**2-rh5**2)/(np.pi*(rs5-rh5))+Zr*b5)))  #m
    Y5 = np.arctan(0.5*(np.tan(Beta4)+np.tan(Beta5)))
    C = Zr/np.cos(Y5)
    if (r4-rs5)/b5>=0.2:
        KpCETI = 0.11
    else:
        KpCETI = 0.22
    LossPass = KpCETI*(LH/DH+0.68*((1-(r5/r4)**2)*np.cos(Beta5)/(b5/C))*((W4**2+W5**2)/2))
    
    #Rotor Clearance Losses
    Ca = (1-(rs5/r4))/(Cm4*b4)
    Cr = (rs5/r4)*((Zr-b4)/(Cm5*r5*b5))
    Ka = 0.4
    Kr = 0.75
    Kar = -0.3
    Ea = 0.0003
    Er = 0.0003
    if Ea*Er*Ca*Cr>=0:
        LossTip = (U4**3*NR/(8*np.pi))*(Ka*Ea*Ca+Kr*Er*Cr+Kar*np.sqrt(Ea*Er*Ca*Cr))
    else:
        LossTip = (U4**3*NR/(8*np.pi))*(Ka*Ea*Ca+Kr*Er*Cr)

    #Windage Losses
    Eb = 0.0003
    Kf = 3.7*(Eb/r4)**0.1/Re4s**0.5
    LossWind = Kf*((rho4s+rho5ss)/2)*U4**3*r4**2/(2*mflow*W5**2)

    #Trailing Edge Losses
    tb4 = 0.04*r4
    tb5 = 0.02*r4
    LossTE = rho5ss*W5**2/2*(NR*tb5/(np.pi*(rh5+rs5)*np.cos(Beta5)))**2
    
    #Exit Losses
    LossExit = 0.5*C5**2    # => mungkin untuk ubah dari total jadi static. abaikan dulu

    # => Sum Enthalpy Losses
    TotalLoss = LossInc + LossPass + LossTip + LossWind + LossTE
    # \\\\\\\ <<---------<<----||----->>------------>> ////////


    #Perhitungan Properties considering losses
    h05 = h05ss+ (LossInc0+LossPass+LossTip+LossWind+LossTE)           #bentar masih salah lossnya       #nozzle masih diasumsikan isentropic dan isenthalpic
    h5  = h5ss+ (LossInc+LossPass+LossTip+LossWind+LossTE) 
    p05 = p05ss
    T05 = Props('T','H',h05,'P',p05,fluid)
    p5  = p5ss
    T5  = Props('T','H',h5,'P',p5,fluid)

    #Effisiensi 
    Reaction    = (h4s-h5)/(h01-h05ss)
    Efftt       = ((h01-h05)/(h01-h05ss)-Effreductbladeloading)*100
    Effts       = ((h01-h05)/(h01-h5ss)-Effreductbladeloading)*100

    geomdict    = dict()
    thermodict  = dict()
    veltridict  = dict()
    effdict     = dict()
    lossdict    = dict()

    for i in ('r4','r5','rs5','rh5','b4','b5','Zr','NR','tb4','tb5'):
        geomdict[i]     = globals()[i]
    for i in ('T_1','T_5','P_1','P_5','p04s','p04','p4s','p4','p5ss','p5','p05ss','p05','T05ss','T05','T5ss','T5','rho4s','rho5ss'):
        thermodict[i]   = globals()[i]
    for i in ('C4','Ct4','Cm4','W4','U4','Alpha4','Beta4','C5','Ct5','Cm5','W5','U5','Alpha5','Beta5','Beta4opt','Beta4opt2','Cm5didconverge1','Cm5didconverge2','k1Cm5','k2Cm5'):
        veltridict[i]   = globals()[i]
    for i in ('Reaction','Effts','Efftt'):
        effdict[i]      = globals()[i]
    for i in ('LossInc0','LossInc','LossPass','LossTip','LossWind','LossTE','Effreductbladeloading','LossExit'):
        lossdict[i]     = globals()[i]
    for i in ('Beta4','b4','r4','Zr','rs5','rh5'):
        proceeddict[i]  = globals()[i]
    outputdict  = {
        'geometry'  : geomdict,
        'thermo'    : thermodict,
        'velocity'  : veltridict,
        'efficiency': effdict,
        'losses'    : lossdict,
        'proceed'   : proceeddict
        }
    
    return outputdict

def ComputeR3(tenflow_coeff,tenwork_coeff,k,l,m):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu

    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10

    cycledict=whichcycle (k)       # The cycle to be computed
    globals().update(cycledict)
    gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    globals().update(gparamdict)
    rpm =whatrpm(m)          # rpm at m will be used
    P_1 = P_1*10**6
    P_5 = P_5*10**6

    Cp4 = Props('C','T',T_1,'P',P_1,fluid)
    Cv4 = Props('O','T',T_1,'P',P_1,fluid)
    gamma = Cp4/Cv4
    Rx = 8.31446261815324   #J/K.mol

    #General Properties inlet outlet turbin (Total)
    H_1     = Props('H','T',T_1,'P',P_1,fluid)     #J/kg
    s01     = Props('S','T',T_1,'P',P_1,fluid)     #J/kg.K 
    T_5     = Props('T','P',P_5,'S',s01,fluid)  # =>asumsi nozzle isenthalpy DAN Isentropic
    H_5     = Props('H','T',T_5,'P',P_5,fluid)  # meski pada kenyataannya isenthalpic nozzle tidak isentropic
    DeltaH  = H_1-H_5            #Ideal === Isentropic Total Enthalpy change 

    C0s     = np.sqrt(2*DeltaH)         #Spouting Velocity

    #Perhitungan Properties ideal lain (Total)
    p01     = P_1           #inlet volute [1], Total
    T01     = T_1
    h01     = H_1
    p1      = p01           # inlet turbine, V~0 
    T1      = T_1
    h01     = H_1
    rho1   = Props('D','P',p1,'T',T1,fluid)
    h02s    = H_1           #inlet nozzle [2], Total
    s02s    = s01            #ideal volute === approx. as isentropic
    p02s    = p01
    T02s    = T01
    h03s    = h02s           #outlet nozzle [3], Total
    s03s    = s02s            #ideal nozzle === approx. as isentropic (in Total)
    p03s    = p02s
    T03s    = T02s
    h04s    = h03s           #inlet rotor [4], Total
    s04s    = s03s           #outlet nozzle === inlet rotor
    p04s    = p03s
    T04s    = T03s
    h04     = h04s          # Nozzle isenthalpic but not isentropic
    p05ss   = P_5
    T05ss   = T_5
    h05ss   = H_5
    s05ss   = s04s

    #Segitiga Kecepatan Inlet, m/s, radians
    U4      = np.sqrt(DeltaH/work_coeff)
    Cm4     = U4*flow_coeff
    Ct4     = DeltaH/U4                 # => DeltaH = U4*Ct4-U5*Ct5 ; Alpha5=0 => Ct5=0
    C4      = np.sqrt(Cm4**2+Ct4**2)
    Alpha4  = np.arctan(Ct4/Cm4)
    W4      = np.sqrt(Cm4**2+(U4-Ct4)**2)
    Beta4   = np.arctan((U4-Ct4)/Cm4)

    h4s    = h04s-1/2*Cm4**2
    rho4s   = Props('D','H',h4s,'S',s04s,fluid)
    rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)
    
    
    Ct5 = 0 # => it is predetermined that Alpha5=0
    Alpha5 = 0
    Cm5_0    = 0
    rho5ss_0= rho05ss        # => initial value for iteration
    Cm5ii    = Cm5_0
    rho5ssii= rho5ss_0
    Cm5didconverge1 = False
    Cm5didconverge2 = False
    choked5     = False
    k1Cm5    = 0
    k2Cm5    = 0
    while Cm5didconverge1 == False:
        k1Cm5       = k1Cm5+1             # => iteration amount
        Cm5i        = Cm5ii
        rho5ssi      = rho5ssii
        Cm5ii       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ssi)*Cm4
        # Cm4ii       = mflow/(2*np.pi()*b5*)
        h5ss         = h05ss-1/2*(Cm5ii**2+Ct5**2)
        rho5ssii     = Props('D','H',h5ss,'S',s05ss,fluid)
        errorCm5    = np.abs((Rb5b4*Rr5r4*(rho5ssii/rho4s)*(Cm5ii/Cm4))-1)
        # errorCm4    = mflow/(rho5ssii*Cm4ii*2*np.pi*b4*r4)-1
        if errorCm5 <= 5*1e-3:
            Cm5didconverge1 = True
            Cm5didconverge2 = True
            Cm5     = Cm5ii
            rho5ss   = rho5ssii
            break
        if (rho5ssi*Cm5i-rho5ssii*Cm5ii)*(Cm5i-Cm5ii)<0:
            Cm5      = Cm5ii
            rho5ss  = rho5ssii
            break
        if k1Cm5>200:
            print(f"loop1 iterates too long ({k1Cm5}) at {flow_coeff,work_coeff} with errorCm5 = {errorCm5}")
            break
    while Cm5didconverge2 == False:
        k2Cm5     = k2Cm5 +1         # => iteration amount
        Cm5       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ss)*Cm4
        h5ss       = h05ss-1/2*(Cm5**2+Ct5**2)
        rho5ss     = Props('D','H',h5ss,'S',s05ss,fluid)
        if np.abs(1-Cm5/Props('A','H',h5ss,'S',s05ss,fluid)) < 5*1e-3:
            choked5 = True
            break
        errorCm5  = np.abs((Rb5b4*Rr5r4*(rho5ss/rho4s)*(Cm5/Cm4))-1)
        if errorCm5 <= 5*1e-5:
            Cm5didconverge2 = True
            break
        if k2Cm5>200:
            print(f"loop2 iterates too long ({k1Cm5},{k2Cm5}) at {flow_coeff,work_coeff} with errorCm5 = {errorCm5}")
            break
    h5ss    = h05ss-1/2*(Cm5**2+Ct5**2)
    C5      = np.sqrt(Cm5**2+Ct5**2)
    U5      = U4*Rr5r4
    W5      = np.sqrt(Cm5**2+(U5-Ct5)**2)
    Beta5   = np.arccos(Cm5/W5)
    
    #Perhitungan geometri
    r4      = U4/np.radians(rpm*6)
    r5      = Rr5r4*r4
    b4      = Rb4r4*r4
    b5      = Rb5b4*b4
    rs5     = (2*r5+b5)/2
    rh5     = rs5-b5
    if rh5 < 0.0015:
        print(f"For flow coeff ={flow_coeff} and work coeff={work_coeff} rh5 too small(<1.5mm), adjust gparams")
        return
    Zr      = RZrr4*r4

    mflow   = 2*np.pi*b5*r5*rho5ss*Cm5

    Q5      = mflow/rho05ss
    ns      = np.radians(rpm*6)*np.sqrt(Q5)/DeltaH**0.75
    Efftspred    = 0.81-1.07*(ns-0.55)**2-0.5*(ns-0.55)**3       #predicted total-to-static efficiency

    p04     = p01-rho1*DeltaH*(1-Efftspred)/4
    T04     = Props('T','P',p04,'H',h04s,fluid)
    s04     = Props('S','P',p04,'T',T04,fluid)

    #Perhitungan Properties ideal lain (Static)
    h4s     = h04s-1/2*C4**2
    p4s     = Props('P','H',h4s,'S',s04s,fluid)
    p4      = Props('P','H',h4s,'S',s04,fluid)
    T4s     = Props('T','H',h4s,'S',s04s,fluid)
    rho04s  = Props('D','P',p04s,'T',T04s,fluid)
    rho4s   = Props('D','P',p4s,'T',T4s,fluid)
    rho4sm  = 2*(p04s-p4s)/C4**2
    h4      = h04-1/2*C4**2
    p4      = Props('P','H',h4,'S',s04,fluid)
    T4      = Props('T','H',h4,'S',s04,fluid)
    rho04   = Props('D','P',p04,'T',T04,fluid)
    rho4    = Props('D','P',p4,'H',h04,fluid)
    rho4m   = 2*(p04-p4)/C4**2
    a01     = Props('A','P',p01,'T',T01,fluid)
    a4s     = Props('A','P',p4s,'T',T4s,fluid)
    a4      = Props('A','P',p4,'T',T4,fluid)
    Ma4s    = C4/a4s
    Ma4     = C4/a4
    T5ss    = Props('T','H',h5ss,'S',s05ss,fluid)
    p5ss    = Props('P','H',h5ss,'S',s05ss,fluid)


    Re4s    = rho4s*C4*b4/Props('V','P',p4s,'T',T4s,fluid)
    Re4     = rho4*C4*b4/Props('V','P',p4,'T',T4,fluid)



    S5      = 2*np.pi*r5/NR
    O5      = S5*Cm5/W5

    p5ss    = Props('P','H',h5ss,'S',s05ss,fluid)
    a5ss    = Props('A','H',h5ss,'P',p5ss,fluid)
    Ma5ss   = C5/a5ss


    # \\\\\\\ <<---------<<----||----->>------------>> ////////
    ## Losses Coefficient ##

    #Rotor Incidence Losses 
    Beta4opt2= np.arctan((-1.98/NR)/(1-1.98/NR)*np.tan(Alpha4))
    Beta4opt = np.arctan(np.tan(Alpha4)*(work_coeff-1+2/NR)/work_coeff)  #(Chen)
    LossInc0 = 0.5*(W4**2)*(np.sin(np.abs(np.abs(Beta4)-np.abs(Beta4opt))))**2  #m2/s2
    LossInc  = 0.5*(W4**2)*(np.sin(Beta4)-np.sin(Beta4opt))**2
       
    #Blade loading efficiency (Chen)
    vNu = U4/np.sqrt(2*Cp4*T01*(1-(p5ss/p01)**((gamma-1)/gamma))) #blade/isentropic jet speed ratio
    Effreductbladeloading = flow_coeff**2*vNu**2

    #Rotor Passage Losses ([Uusitalo] from Moustapha PLM3)
    LH = np.pi/4*((Zr-b4/2)+(r4-rh5-b5/2))                                                              #m
    DH = 0.5*((4*np.pi*r4*b4/(2*np.pi*r4+Zr*rh5))+((2*np.pi*(rs5**2-rh5**2)/(np.pi*(rs5-rh5))+Zr*b5)))  #m
    Y5 = np.arctan(0.5*(np.tan(Beta4)+np.tan(Beta5)))
    C = Zr/np.cos(Y5)
    if (r4-rs5)/b5>=0.2:
        KpCETI = 0.11
    else:
        KpCETI = 0.22
    LossPass = KpCETI*(LH/DH+0.68*((1-(r5/r4)**2)*np.cos(Beta5)/(b5/C))*((W4**2+W5**2)/2))
    
    #Rotor Clearance Losses
    Ca = (1-(rs5/r4))/(Cm4*b4)
    Cr = (rs5/r4)*((Zr-b4)/(Cm5*r5*b5))
    Ka = 0.4
    Kr = 0.75
    Kar = -0.3
    Ea = 0.0003
    Er = 0.0003
    if Ea*Er*Ca*Cr>=0:
        LossTip = (U4**3*NR/(8*np.pi))*(Ka*Ea*Ca+Kr*Er*Cr+Kar*np.sqrt(Ea*Er*Ca*Cr))
    else:
        LossTip = (U4**3*NR/(8*np.pi))*(Ka*Ea*Ca+Kr*Er*Cr)
    #Windage Losses # disk friction losses (fiaschi, 2015: 4.2.5)
    Eb = 0.0003
    Kf = 3.7*(Eb/r4)**0.1/Re4s**0.5
    LossWind = Kf*((rho4s+rho5ss)/2)*U4**3*r4**2/(2*mflow*W5**2)

    #Trailing Edge Losses
    # if 0.04*r4>0.001:
    tb4 = 0.04*r4
    # else:
    # tb4 = 0.001
    # if 0.02*r4>0.001:
    tb5 = 0.02*r4
    # else:
    #     tb5 = 0.001
    LossTE = rho5ss*W5**2/2*(NR*tb5/(np.pi*(rh5+rs5)*np.cos(Beta5)))**2
    
    #Exit Losses
    LossExit = 0.5*C5**2    # => mungkin untuk ubah dari total jadi static. abaikan dulu

    # => Sum Enthalpy Losses
    TotalLoss = LossInc + LossPass + LossTip + LossWind + LossTE
    # \\\\\\\ <<---------<<----||----->>------------>> ////////


    #Perhitungan Properties considering losses
    h05 = h05ss+ (LossInc0+LossPass+LossTip+LossWind+LossTE)             #nozzle masih diasumsikan isentropic dan isenthalpic
    h5  = h5ss+ (LossInc+LossPass+LossTip+LossWind+LossTE  )   
    p05 = p05ss
    p5  = Props('P','H',h5,'S',s05ss,fluid)
    T05 = Props('T','H',h05,'P',p05,fluid)
    T5  = Props('T','H',h5,'P',p5,fluid)

    #Effisiensi 
    Reaction    = (h4s-h5)/(h01-h05ss)
    Efftt       = ((h01-h05)/(h01-h05ss)-0)*100
    Effts       = ((h01-h05)/(h01-h5ss)-Effreductbladeloading)*100

    geomdict    = dict()
    thermodict  = dict()
    veltridict  = dict()
    effdict     = dict()
    lossdict    = dict()
    proceeddict = dict()

    for i in ('r4','r5','rs5','rh5','b4','b5','Zr','NR','tb4','tb5'):
        geomdict[i]     = globals()[i]
    for i in ('T_1','T_5','P_1','P_5','p04s','p04','p4s','p4','p5ss','p5','p05ss','p05','T05ss','T05','T5ss','T5','rho4s','rho5ss','mflow','h4s','h04s'):
        thermodict[i]   = globals()[i]
    for i in ('C4','Ct4','Cm4','W4','U4','Alpha4','Beta4','C5','Ct5','Cm5','W5','U5','Alpha5','Beta5','Beta4opt','Beta4opt2','Cm5didconverge1','Cm5didconverge2','k1Cm5','k2Cm5'):
        veltridict[i]   = globals()[i]
    for i in ('Reaction','Effts','Efftt'):
        effdict[i]      = globals()[i]
    for i in ('LossInc0','LossInc','LossPass','LossTip','LossWind','LossTE','Effreductbladeloading','LossExit'):
        lossdict[i]     = globals()[i]
    for i in ('Beta4','b4','r4','Zr','rs5','rh5'):
        proceeddict[i]  = globals()[i]
    outputdict  = {
        'geometry'  : geomdict,
        'thermo'    : thermodict,
        'velocity'  : veltridict,
        'efficiency': effdict,
        'losses'    : lossdict,
        'proceed'   : proceeddict
        }

    return outputdict

# func to find effts value from fitted func of 
# n-set of cycle,gparams,rpm
# for specific 10xflow and 10xwork coeff point
def fiteffts(tenflow_coeff,tenwork_coeff,n): 
    p = whichfitfun(n)
    effts=0
    for i in range(0,6):
        for j in range (0,6):
            effts=effts+p[i][j]*tenflow_coeff**i*tenwork_coeff**j
    return effts

# def efftsforn(tenflow_coeff,tenwork_coeff):
#     p 
#     fiteffts=0
#     for i in range(0,6):
#         for j in range (0,6):
#             fiteffts=fiteffts+p[i][j]*tenflow_coeff**i*tenwork_coeff**j
# func to find optimum effts value of the fitted func
def mfiteffts(x):
    tenflow_coeff   = x[0]
    tenwork_coeff   = x[1]
    n               = x[2]
    return fiteffts(tenflow_coeff,tenwork_coeff,n)*-1
def constraint1(x):
    return x[0]-0
def constraint2(x):
    return x[1]-0
def constraint3(x):
    return x[2]-0
def optfiteffts(n):

    tenflowb=(0,12)  # Constraint: 0<=10*flowcoeff<=9
    tenworkb=(0,25)  # Constraint: 0<=10*workcoeff<=20
    nb=(n,n)
    
    
    constr1 = {'type': 'ineq', 'fun':constraint1 }
    constr2 = {'type': 'ineq', 'fun':constraint2 }
    constr3 = {'type': 'eq', 'fun': constraint3 }
    constr  = [constr1,constr2,constr3]
    bnds=(tenflowb,tenworkb,nb)
    initval=[0.5,3,n]

    opteffts=optimize.minimize(mfiteffts,initval,method='BFGS',bounds=bnds,constraints=constr)
    return opteffts
    


#Create 2D Quasi
def QuasiNorm(indict,n,Z5):    # INPUT => b4,r4,Zr,rs5,rh5,n,Z5; n dan Z5 harus dinyatakan
    
    global Z,r,Zh,rh,m,Ash,Bsh,Csh,Dsh,Esh,Fsh
    #local Beta4,b4,r4
    #translist=['Beta4','b4','r4','Zr','rs5','rh5']#,'n','Z5']

    locals().update(indict)

    # for i in range(len(translist)):
    #     globals()[translist[i]]=inlist[i]
    
    Betha4 = Beta4#np.arctan((Ct4-U4)/Cm4) # tan kan depan/samping
    QuasiSec=50
    # dZ=Zr/QuasiSec
    # dzb4=abs(Zr-b4)

    # C2=(r4-rs5)/((dZ-b4)**n)
    # Zrb4=abs(b4-Zr)
    # r5=(rs5+rh5)/2
    # b5=(rs5-rh5)
    dzb4=abs(Zr-b4)
    dZ=dzb4/QuasiSec
    C2=(r4-rs5)/((Zr-b4)**n)
    r5=(rs5+rh5)/2
    b5=(rs5-rh5)
    #Shroud Sections
    # Z=[]
    # Z.append(0)
    # m=[]
    # m.append(0)
    # for i in range(0,QuasiSec-1):
    #     Z.append(Z[i]+dZ)
    # Z=np.array(Z)
    # dzb4=dZ-b4
    # for i in [Z]:
    #     Epsi = i/dzb4
    #     r=rs5+(r4-rs5)*Epsi**n
    # for i in [Z]:
    #     fz=np.sqrt(1+(C2*n*(i-(i-1))**(n-1))**2)
    # m.append((dZ/3)*(fz[0]+4*fz[1]) )
    # for i in range (2, len(fz)):
    #     if (i % 2) == 0:
    #         m.append(m[i-1]+4*(fz[i])/3)
    #     else:
    #         m.append(m[i-1]+2*(fz[i])/3)                            #perbaiki value fz #check OK

    Z=[]
    Z.append(0)
    m=[]
    m.append(0)
    for i in range(0,QuasiSec):
        Z.append(Z[i]+dZ)

    Epsi=[]
    r=[]
    fz=[]
    for i in range(0,len(Z)):
        Epsi.append(Z[i]/dzb4)
        r.append(rs5+(r4-rs5)*Epsi[i]**n)
    for i in range(0,len(Z)):
        fz.append(np.sqrt(1+(C2*n*(Z[i]-(Z5))**(n-1))**2))

    m.append((dZ/3)*(fz[0]+4*fz[1]))
    for i in range (2, len(fz)):
        if (i % 2) == 0:
            m.append(m[i-1]+4*(fz[i])*dZ/3) 
        else:
            m.append(m[i-1]+2*(fz[i])*dZ/3)


    #Hub Sections
    # Zh=[]
    # Zh.append(0)
    # mh=[]
    # mh.append(0)
    # rh=[]
    # fzh=[]
    # for i in range(1,QuasiSec):
    #     Zh.append(Zh[i-1]+dZ)
    #     Ephi = Zh[i]/(dZ-b4)
    # for i in range(0,QuasiSec):
    #     rh.append(rh5+Zr-np.sqrt((Zr**2)-(Zh[i]**2)))
    #     fzh.append(np.sqrt(1+(C2*n*(Zh[i]-Z5)**(n-1))**2))
    # mh.append((dZ/3)*(fzh[0]+4)/fzh[0])                      #Perbaiki
    # for i in range(1,len(fzh)):
    #     mh.append(mh[i]+2*(i)/3)                                 #Perbaiki value mh 

    Zh=[]
    Zh.append(0)
    mh=[]
    mh.append(0)
    rh=[]
    dZh=Zr/50
    for i in range(0,QuasiSec-5):
        Zh.append(Zh[i]+dZh)
    for i in range(0,4):
        Zh.append(Zh[QuasiSec-5])

    for i in range(0,len(Zh)-5):
        rh.append(rh5+Zr-np.sqrt((Zr**2)-(Zh[i]**2)))
    for i in range(QuasiSec-5,QuasiSec):
        rh.append((r[QuasiSec]-rh[QuasiSec-6])/5+rh[i-1])

    for i in range(1,len(Z)):
        mh.append(i/QuasiSec*90/360*np.pi*Zr*2)


    Betha5  = np.arccos(Cm5/W5)
    Betha5s = np.arctan(np.tan(Betha5)/rs5)
    Betha5h = np.arctan(np.tan(Betha5)/rh5)
    ms      = m[QuasiSec-1]
    Ash     = 1/np.tan(Betha5s)/rs5
    Bsh     = (1/np.tan(Betha4)-1/np.tan(Betha5s))/(ms**2)
    Csh     = -Bsh/(2*ms)
    Tetha4  = ms/2*(1/np.tan(Betha4)/(r4+1/np.tan(Betha5))/rs5)
    mhh     = mh[50]
    Dsh     = 1/np.tan(Betha5h)/rh5
    Esh     = 3*Tetha4/(mhh**2)-(1/mhh)*(2*1/np.tan(Betha5h))
    Fsh     = 1/mhh**2*(1/np.tan(Betha5h)/rh5+1/np.tan(Betha4)/r4)-2*Tetha4/mhh**3
    #outlist = [Z,r,Zh,rh,m,Ash,Bsh,Csh,Dsh,Esh,Fsh]
    outputdict = dict()
    for i in ('Z','r','Zh','rh','m','Ash','Bsh','Csh','Dsh','Esh','Fsh'):
        outputdict[i]=globals()[i]

    return outputdict
    # OUTPUT => Z,r,Zh,rh,m,Ash,Bsh,Csh,Dsh,Esh,Fsh

#Perhitungan sudut   
def Zrregs(indict): # INPUT => Z,r,Zh,rh,m,Ash,Bsh,Csh,Dsh,Esh,Fsh
    global phis,tethas,Bethacs,phih,tethah,Bethach

    locals().update(indict)
    #Shroud
    

    # dms=[]
    # dms.append(0)
    # for i in range(1,len(Z)):
    #     dms=np.sqrt((Z[i]-Z[i-1])**2+(r[i]-r[i-1])**2)
    # msi=[]
    # msi.append(0)
    # for i in range(1,len(Z)):
    #     msi.append(msi[i-1]+dms)
    # polynom = 6             #Coeficient harus 6 (Belum otomatis)
    # mymodel = np.poly1d(np.polyfit(Z,r,polynom))
    # #need to regress the mymodel
    # #Define the Function of phi
    # coeffphis=[]
    # phis=[]
    # tethas=[]
    # Bethacs=[]
    # for coeff in mymodel:
    #     coeffphis.append(coeff)
    # for i in range(0,len(msi)):
    #     phis.append(6*coeffphis[5]*msi[i]**5+5*coeffphis[4]*msi[i]**4+4*coeffphis[3]*msi[i]**3+3*coeffphis[2]*msi[i]**2+2*coeffphis[1]*msi[i]**1+coeffphis[0])
    #     tethas.append(msi[i]*(Ash+Bsh**3+Csh**4))
    #     Bethacs.append(1/np.arctan(r[i]*(Ash+3*Bsh**2+4*Csh**3)))

    dms=[]
    dms.append(0)
    for i in range(1,len(Z)):
        dms.append(np.sqrt((Z[i]-Z[i-1])**2+(r[i]-r[i-1])**2))
    msi=[]
    msi.append(0)
    for i in range(1,len(Z)):
        msi.append(msi[i-1]+dms[i])

    phis=[]
    tethas=[]
    Bethacs=[]
    phis.append(0)
    tethas.append(0)

    for i in range(1,len(msi)):
        if ((r[i]-r[i-1])/(msi[i]-msi[i-1]))>1:
            phis.append(np.arcsin(1))
        else:
            phis.append(np.arcsin((r[i]-r[i-1])/(msi[i]-msi[i-1])))
        tethas.append(Ash*msi[i]+Bsh*msi[i]**3+Csh*msi[i]**4)
    for i in range(0,len(msi)):
        Bethacs.append(np.arctan(1/(r[i]*(Ash+3*Bsh*msi[i]**2+4*Csh*msi[i]**3))))

    #Hub
    # dmh=[]
    # dmh.append(0)
    # for i in range(1,len(Zh)):
    #     dmh.append(np.sqrt((Zh[i]-Zh[i-1])**2+(rh[i]-rh[i-1])**2))
    # mhi=[]
    # mhi.append(0)
    # for i in range(1,len(dmh)):
    #     mhi.append(mhi[i-1]+dmh[i])
    # polynomh = 6
    # mymodelh = np.poly1d(np.polyfit(Zh,rh,polynomh))
    # coeffphih=[]
    # phih=[]
    # tethah=[]
    # Bethach=[]
    # for coeffh in mymodelh:
    #     coeffphih.append(coeffh)
    # for i in range(0,len(mhi)):
    #     phih.append(6*coeffphih[5]*mhi[i]**5+5*coeffphih[4]*mhi[i]**4+4*coeffphih[3]*mhi[i]**3+3*coeffphih[2]*mhi[i]**2+2*coeffphih[1]*mhi[i]**1+coeffphih[0])
    #     tethah.append(1/np.tan(rh[i]*(mhi[i]*(Dsh+Esh**2+Fsh**3))))
    #     Bethach.append(1/np.arctan(rh[i]*(Dsh+3*Esh**2+4*Fsh**3)))

    dmh=[]
    dmh.append(0)
    for i in range(1,len(Zh)):
        dmh.append(np.sqrt((Zh[i]-Zh[i-1])**2+(rh[i]-rh[i-1])**2))
    mhi=[]
    mhi.append(0)
    for i in range(1,len(Zh)):
        mhi.append(mhi[i-1]+dmh[i])

    phih=[]
    tethah=[]
    Bethach=[]
    phih.append(0)
    tethah.append(0)

    for i in range(1,len(mhi)):
        if (rh[i]-rh[i-1])/(mhi[i]-mhi[i-1])>1:
            phih.append(np.arcsin(1))
        else:
            phih.append(np.arcsin((rh[i]-rh[i-1])/(mhi[i]-mhi[i-1])))
        tethah.append(mhi[i]*Dsh+Esh*mhi[i]**2+Fsh*mhi[i]**3)
    for i in range(0,len(mhi)):
        Bethach.append(np.arctan(1/(rh[i]*(Dsh+3*Esh*mhi[i]**2+4*Fsh*mhi[i]**3))))

    outputdict = dict()
    for i in ('phis','tethas','Bethacs','phih','tethah','Bethach'):
        outputdict[i]=globals()[i]
    return outputdict
    # OUTPUT => return(phis,tethas,Bethacs,phih,tethah,Bethach)

#Meridional Coordinate
def meridional(indict): # INPUT => phis,tethas,Bethacs,phih,tethah,Bethach
    global X,Y,Z,Xh,Yh,Zh
    locals().update(indict)

    X=[]
    Y=[]
    Xh=[]
    Yh=[]
    #Shroud
    for i in range(0,len(tethas)):
        X.append(r[i]*np.sin(tethas[i]))
        Y.append(r[i]*np.cos(tethas[i]))
    Z=Z
    #Hub
    for i in range(0,len(tethah)):
        Xh.append(rh[i]*np.sin(tethah[i]))
        Yh.append(rh[i]*np.cos(tethah[i]))
    Zh=Zh
    outputdict=dict()
    for i in ('X','Y','Z','Xh','Yh','Zh'):
        outputdict[i]=globals()[i]
    return outputdict
    # OUTPUT => return(X,Y,Z,Xh,Yh,Zh)

#Vector Component
def VectorComp(indict): # INPUT => X,Y,Z,Xh,Yh,Zh
    global Txs,Tys,Tzs,Txh,Tyh,Tzh
    locals().update(indict)
    #Length Calculation
    L=[]
    Bx=[]
    By=[]
    Bz=[]
    Sxs=[]
    Sys=[]
    Szs=[]
    Sxh=[]
    Syh=[]
    Szh=[]
    Txs=[]
    Tys=[]
    Tzs=[]
    Txh=[]
    Tyh=[]
    Tzh=[]
    for i in range(0,len(X)-1):
        L.append(np.sqrt((X[i]-Xh[i])**2+(Y[i]-Yh[i])**2+(Z[i]-Zh[i])**2))

    #Vector B Hub and Shroud
    
    for i in range(0,len(X)-1):
        Bx.append((X[i]-Xh[i])/L[i])
        By.append((Y[i]-Yh[i])/L[i])
        Bz.append((Z[i]-Zh[i])/L[i])

    #Vector S Shroud
        Sxs.append(np.sin(phis[i])*np.sin(tethas[i])*np.sin(Bethacs[i])+np.cos(tethas[i])*np.cos(Bethacs[i]))
        Sys.append(np.cos(phis[i])*np.sin(tethas[i])*np.sin(Bethacs[i])-np.sin(tethas[i])*np.cos(Bethacs[i]))
        Szs.append(np.sin(tethas[i])*np.sin(Bethacs[i]))

    #Vector S Hub
        Sxh.append(np.sin(phih[i])*np.sin(tethah[i])*np.sin(Bethach[i])+np.cos(tethah[i])*np.cos(Bethach[i]))
        Syh.append(np.cos(phih[i])*np.sin(tethah[i])*np.sin(Bethach[i])-np.sin(tethah[i])*np.cos(Bethach[i]))
        Szh.append(np.sin(tethah[i])*np.sin(Bethach[i]))

    #Vector T Shroud
        Txs.append(Szs[i]*By[i]-Sys[i]*Bz[i])
        Tys.append(Sxs[i]*Bz[i]-Szs[i]*Bx[i])
        Tzs.append(Sys[i]*Bx[i]-Sxs[i]*By[i])
    #Vector T Hub
        Txh.append(Szh[i]*By[i]-Syh[i]*Bz[i])
        Tyh.append(Sxh[i]*Bz[i]-Szh[i]*Bx[i])
        Tzh.append(Syh[i]*Bx[i]-Sxh[i]*By[i])

    outputdict=dict()
    for i in ('Txs','Tys','Tzs','Txh','Tyh','Tzh'):
        outputdict[i]=globals()[i]
    return outputdict
    # OUTPUT => Txs,Tys,Tzs,Txh,Tyh,Tzh,tb4,tb5


#3D Coordinate

def Coord3D(indict,tb4,tb5): # INPUT =>Txs,Tys,Tzs,Txh,Tyh,Tzh,tb4,tb5
    # global XcorS,YcorS,ZcorS,XcorH,YcorH,ZcorH
    global XcorSP,YcorSP,ZcorSP,XcorSN,YcorSN,ZcorSN,XcorHP,YcorHP,ZcorHP,XcorHN,YcorHN,ZcorHN
    locals().update(indict)
    #Tebal Sudu
    tb =[]
    tb.append(tb4)
    for i in range(1,len(Txs)):
        tb.append(tb[i-1]+(abs(tb4-tb5)/len(Txs)))

    #Shroud 
    # XcorS=[]
    # YcorS=[]
    # ZcorS=[]
    # for i in range(0,len(Txs)):
    #     XcorS.append(X[i]+0.5*Txs[i]*tb[i])
    #     XcorS.append(X[i]-0.5*Txs[i]*tb[i])
    #     YcorS.append(Y[i]+0.5*Tys[i]*tb[i])
    #     YcorS.append(Y[i]-0.5*Tys[i]*tb[i])
    #     ZcorS.append(Z[i]+0.5*Tzs[i]*tb[i])
    #     ZcorS.append(Z[i]-0.5*Tzs[i]*tb[i])

    XcorSP=[]
    YcorSP=[]
    ZcorSP=[]
    XcorSN=[]
    YcorSN=[]
    ZcorSN=[]
    for i in range(0,len(Txs)):
        XcorSP.append(X[i]+0.5*Txs[i]*tb[i])
        YcorSP.append(Y[i]+0.5*Tys[i]*tb[i])
        ZcorSP.append(Z[i]+0.5*Tzs[i]*tb[i])
        
    for i in range(0,len(Txs)):
        XcorSN.append(X[i]-0.5*Txs[i]*tb[i])
        YcorSN.append(Y[i]-0.5*Tys[i]*tb[i])
        ZcorSN.append(Z[i]-0.5*Tzs[i]*tb[i])

    #Hub
    # XcorH=[]
    # YcorH=[]
    # ZcorH=[]
    # for i in range(0,len(Txs)):
    #     XcorH.append(Xh[i]+0.5*Txh[i]*tb[i])
    #     XcorH.append(Xh[i]-0.5*Txh[i]*tb[i])
    #     YcorH.append(Yh[i]+0.5*Tyh[i]*tb[i])
    #     YcorH.append(Yh[i]-0.5*Tyh[i]*tb[i])
    #     ZcorH.append(Zh[i]+0.5*Tzh[i]*tb[i])
    #     ZcorH.append(Zh[i]-0.5*Tzh[i]*tb[i])

    #Hub
    XcorHP=[]
    YcorHP=[]
    ZcorHP=[]
    XcorHN=[]
    YcorHN=[]
    ZcorHN=[]
    for i in range(0,len(Txs)):
        XcorHP.append(Xh[i]+0.5*Txh[i]*tb[i])
        YcorHP.append(Yh[i]+0.5*Tyh[i]*tb[i])
        ZcorHP.append(Zh[i]+0.5*Tzh[i]*tb[i])
 
    for i in range(0,len(Txs)):
        XcorHN.append(Xh[i]-0.5*Txh[i]*tb[i])
        YcorHN.append(Yh[i]-0.5*Tyh[i]*tb[i])
        ZcorHN.append(Zh[i]-0.5*Tzh[i]*tb[i])    

    outputdict=dict()
    # for i in ('XcorS','YcorS','ZcorS','XcorH','YcorH','ZcorH'):
    for i in ('XcorSP','YcorSP','ZcorSP','XcorSN','YcorSN','ZcorSN','XcorHP','YcorHP','ZcorHP','XcorHN','YcorHN','ZcorHN'):
        outputdict[i]=globals()[i]
    return outputdict
    # OUTPUT => return (XcorS,YcorS,ZcorS,XcorH,YcorH,ZcorH)

# def proceedR(savenumber):
#     n=4
#     Z5=0
#     QuasiNorm()
#     Zrregs()
#     meridional()
#     VectorComp()
#     Coord3D()
#     dict = {'xshroud':XcorS,'yshroud':YcorS,'zshroud':ZcorS,'xhub':XcorH,'yhub':YcorH,'zhub':ZcorH}
#     df3drotor=pd.DataFrame(dict)
#     savetoas=os.path.join(ROOT_DIR,"outputs",f"rotor{savenumber}coordinate.csv")
#     df3drotor.to_csv(savetoas,index=False)



