import numpy as np
import CoolProp
from CoolProp.CoolProp import PropsSI as Pr

def Compute(flow_coeff,work_coeff,fluid):
    global H_1,H_5,T_1,P_1
    global NR,r4,Alpha4,b4,Ct4,rho4,mflow
    #Dimensionless Number Input

    #Alpha4     = 80.42                                     #Fluid inlet angle
    #flow_coeff =   0.2633879                               #Flow Coefficient (DeltHid/U^2)
    #work_coeff =  1.11609232                               #Work Coefficient (Cm/U)
    rath5 =  0.4
    rats5 =  0.64
    rpm =  50000
    Zratio =  0.36                                         #Axial Length Ratio vs Diameter
    NR =  16 #random.uniform(4,20)                         #Blade Number
    #Thermodynamic Condition
    P_1 = 1.2                                                #MPa
    P_5 = 0.6                                                #MPa
    T_1 = 373.15                                             #K
    Cp4 = Pr('C','T',T_1,'P',P_1*1e6,fluid)                #- 
    Cv4 = Pr('O','T',T_1,'P',P_1*1e6,fluid)                #-   
    gamma = Cp4/Cv4                                       #Heat Capacity Ratio, Affecting fluid Properties and Mach Number, function of T,P, and fluid
    Rx = 8.31446261815324                                 #J/K⋅mol

    #Targeted Power and mflow
    mflow = 0.5                   #kg/s

    #General Properties IDEAL
    H_1 = Pr('H','T',T_1,'P',P_1*1e6,fluid)/1000                                  #Entalphy Inlet (kJ/kg)
    sa1 = Pr('S','H',H_1*1000,'P',P_1*1e6,fluid)                                  #Entalphy Outlet (kJ/kg)
    T_5 = Pr('T','P',P_5*1e6,'S',sa1,fluid)  
    H_5 = Pr('H','T',T_5,'P',P_5*1e6,fluid)/1000                                                              
    Delth = H_1-H_5
    rpm_rad = rpm/(60/(2*np.pi))  #rad/s
    rho4 = Pr('D','T',T_1,'P',P_1*1e6,fluid)      #kg/m3
    rho5 = Pr('D','T',T_5,'P',P_5*1e6,fluid)      #kg/m3

    #Segitiga Kecepatan Inlet
    C0s     = np.sqrt(2*Delth*1000)
    #Beta4=0 & Alpha5=0 => U^2 = DeltaH = (C0s^2)/2. Yang dibawah ini aneh, karena harusnya work coeff selalu 1
    U4      = 0.707*C0s                                 
    Cm4     = U4*flow_coeff                             #m/s 
    Ct4     = U4*work_coeff                             #m/s at Radial Turbine, Ct5 = 0, thus
    C4      = np.sqrt(Cm4**2+Ct4**2)                    #m/s
    W4      = np.sqrt(Cm4**2+(U4-Ct4)**2)               #m/s
    Alpha4  = np.arctan(Ct4/Cm4)*180/np.pi              #Degree
    Beta4   = np.arctan(Cm4/(U4-Ct4))*180/np.pi         #Degree

    #Perhitungan Geometri
    r4 = U4/rpm_rad                                 #m/s
    D4 = 2*r4                                       #m
    Zr  = Zratio*r4                                 #m
    rh5 = rath5*r4                                  #m
    rs5 = rats5*r4                                  #m
    r5  = np.sqrt(rh5**2+rs5**2)    # harusnya dijumlah bagi 2                #m

    #Segitiga Kecepatan Outlet
    Cm5 = Cm4
    Ct5 = 0
    C5 = np.sqrt(Cm5**2+Ct5**2)
    U5 = rpm_rad*r5
    W5 = np.sqrt(Cm5**2+U5**2)
    Beta5 = np.arccos(U5/W5)*180/np.pi #kok pakai arccos

    b5 = rs5-rh5
    b4 = mflow/(2*np.pi*r4*rho4*Cm4)
    deltho=(Ct4*U4-Ct5*U5)/1000
    #kalkulasi efisiensi TS   
    h5s = H_5+0.5*(C5**2)/1000
    H_4  = H_1                                                                     #Adiabatic Assumption
    h4s = H_4+0.5*(C4**2)/1000
    rho_1 = Pr('D','T',T_1,'P',P_1*1e6,fluid)                                     #Density di inlet turbin (kg/m^3)
    rho_5 = Pr('D','T',T_5,'P',P_5*1e6,fluid) 
    Pt_4 = (P_1*1000-(rho_1*(H_1-H_5)/1000*(1-Rx))/4)/1000                        #Tekanan Total Inlet Rotor

    #Perhitungan Properti Lain
    a1= np.sqrt(gamma*Rx*T_1)  #Local Speed of Sound
    Ma1 = U4/a1               #Mach Number
    Re4 = Pr('V','T',T_1,'P',P_1*1e6,fluid)/(Pr('D','T',T_1,'P',P_1*1e6,fluid)*U4*D4) 

    #Losses Coefficient
    #Losses Turbin  Case studies on performance prediction of radial inflow turbine under multiple rotor structures and organic fluids Wenyu 2021

    Re = C4*b4/Pr('viscosity','H',H_5*1000,'P',P_1*1e6,fluid)

    #Rotor Incidence
    Beta4opt = np.arctan((-1.98/NR)/NR/(1-1.98/NR)*np.tan(np.radians(Alpha4)))*180/np.pi
    LossInc = 1e-3*0.5*(W4**2)*(np.sin(np.radians(np.abs(Beta4-Beta4opt))))**2                               #m^2/s^2

    #Rotor Passage
    LH = np.pi/4*((Zr-b4/2)+(r4-rh5-b5/2))                                                              #m
    DH = 0.5*((4*np.pi*r4*b4/(2*np.pi*r4+Zr*rh5))+((2*np.pi*(rs5**2-rh5**2)/(np.pi*(rs5-rh5))+Zr*b5)))  #m
    Y5 = np.arctan(0.5*(np.tan(Beta4)+np.tan(Beta5)))
    C = Zr/np.cos(np.radians(Y5))
    if (r4-rs5)/b5>=0.2:
        KpCETI = 0.11
    else:
        KpCETI = 0.22
        
    LossPass = 1e-3*KpCETI*(LH/DH+0.68*((1-(r5/r4)**2)*np.cos(Beta5*np.pi/180)/(b5/C))*((W4**2+W5**2)/2)) 

    #Rotor Clearance Losses
    Ca = (1-(rs5/r4))/(Cm4*b4)
    Cr = (rs5/r4)*((Zr-b4)/(Cm5*r5*b5))
    Ka = 0.4
    Kr = 0.75
    Kar = -0.3
    Ea = 0.0003
    Er = 0.0003
    LossTip = 1e-3*(U4**3*NR/(8*np.pi))*(Ka*Ea*Ca+Kr*Er*Cr+Kar*np.sqrt(Ea*Er*Ca*Cr))

    #Windage Losses
    Eb = 0.0003
    Kf = 3.7*(Eb/r4)**0.1/Re4**0.5
    LossWind = 1e-3*Kf*((rho4+rho5)/2)*U4**3*r4**2/(2*mflow*W5**2)

    #Trailing Edge Losses
    tb4 = 0.04*r4
    tb5 = 0.02*r4
    LossTE = 1e-3*rho5*W5**2/2*(NR*tb5/(np.pi*(rh5+rs5)*np.cos(Beta5*np.pi/180)))**2

    #Exit Losses Losses
    LossExit = 1e-3*0.5*C5**2

    #Sum Losses
    Total_Loss = LossInc+LossTE+LossTip

    #Efficiency Calclations
    #Efficiency_TS = (h4s-h5s)/(h4s-h5s)
    #Losses_Percentage = Total_Loss/(H_4-h5s)
    #Efficiency = Efficiency_TS-Losses_Percentage
    Reaction = (H_4-h5s)/(h4s-h5s)
    s5 = 2*np.pi*r5/NR
    o5 = Cm5*s5/W5

    ha1=Pr('H','T',T_1,'P',P_1*1e6,fluid)/1000#+0.5*(C4**2)/1000
    ha3=Pr('H','T',T_5,'P',P_5*1e6,fluid)/1000#+0.5*(C5**2)/1000
    saa1=Pr('S','H',ha1*1000,'P',P_1*1e6,fluid)
    ha3ss=Pr('H','S',saa1,'P',P_5*1e6,fluid)/1000-0.5*(C5**2)/1000
    Effx = (ha1-ha3)/(ha1-ha3ss)- Total_Loss/(ha1-ha3ss)
    #t1  = Pr('T','P',P_1*1e6,'H',ha1*1000,fluid)
    #t3  = Pr('T','P',P_5*1e6,'H',ha3*1000,fluid)
    #t3s = Pr('T','P',P_5*1e6,'H',ha3ss*1000,fluid)

    return Effx,Reaction,NR,r4,Alpha4,b4,Ct4,rho4,mflow

    print(r4)
    Compute(0.2,1.1,'R245fa')
    print(Total_Loss)
    