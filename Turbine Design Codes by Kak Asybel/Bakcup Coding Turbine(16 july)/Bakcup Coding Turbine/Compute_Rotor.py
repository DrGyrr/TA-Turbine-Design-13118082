import numpy as np
import CoolProp
from CoolProp.CoolProp import PropsSI as Pr

#Compute Initial Design
def Compute(flow_coeff,work_coeff,fluid):
    global H_1,H_5,T_1,P_1
    global NR,r4,Alpha4,b4,Ct4,rho4,mflow,Cm4,Cm5,C4,C5,Ct5,Ct5,W4,W5,U4,U5,tb4,tb5
    rath5 =  0.4
    rats5 =  0.64
    rpm =  50000
    Zratio =  0.36                                         
    NR =  16 
    P_1 = 1.2                                              #MPa
    P_5 = 0.6                                              #MPa
    T_1 = 373.15                                           #K
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
    rpm_rad = rpm/(60/(2*np.pi))                                                  #rad/s
    rho4 = Pr('D','T',T_1,'P',P_1*1e6,fluid)                                      #kg/m3
    rho5 = Pr('D','T',T_5,'P',P_5*1e6,fluid)                                      #kg/m3

    #Segitiga Kecepatan Inlet
    C0s     = np.sqrt(2*Delth*1000)
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
    r5  = np.sqrt(rh5**2+rs5**2)                    #m

    #Segitiga Kecepatan Outlet
    Cm5 = Cm4
    Ct5 = 0
    C5 = np.sqrt(Cm5**2+Ct5**2)
    U5 = rpm_rad*r5
    W5 = np.sqrt(Cm5**2+U5**2)
    Beta5 = np.arccos(U5/W5)*180/np.pi

    b5 = rs5-rh5
    b4 = mflow/(2*np.pi*r4*rho4*Cm4)
    deltho=(Ct4*U4-Ct5*U5)/1000
    #kalkulasi efisiensi TS   
    h5s = H_5+0.5*(C5**2)/1000
    H_4  = H_1                                                                    #Adiabatic Assumption
    h4s = H_4+0.5*(C4**2)/1000
    rho_1 = Pr('D','T',T_1,'P',P_1*1e6,fluid)                                     #Density di inlet turbin (kg/m^3)
    rho_5 = Pr('D','T',T_5,'P',P_5*1e6,fluid) 
    Pt_4 = (P_1*1000-(rho_1*(H_1-H_5)/1000*(1-Rx))/4)/1000                        #Tekanan Total Inlet Rotor

    #Perhitungan Properti Lain
    a1= np.sqrt(gamma*Rx*T_1)       #Local Speed of Sound
    Ma1 = U4/a1                     #Mach Number
    Re4 = Pr('V','T',T_1,'P',P_1*1e6,fluid)/(Pr('D','T',T_1,'P',P_1*1e6,fluid)*U4*D4) 

    #Losses Coefficient
    #Losses Turbin  Case studies on performance prediction of radial inflow turbine under multiple rotor structures and organic fluids Wenyu 2021

    Re = C4*b4/Pr('viscosity','H',H_5*1000,'P',P_1*1e6,fluid)

    #Rotor Incidence
    Beta4opt = np.arctan((-1.98/NR)/NR/(1-1.98/NR)*np.tan(np.radians(Alpha4)))*180/np.pi
    LossInc = 1e-3*0.5*(W4**2)*(np.sin(np.radians(np.abs(Beta4-Beta4opt))))**2                               #m^2/s^2

    #Rotor Clearance Losses
    Ca = (1-(rs5/r4))/(Cm4*b4)
    Cr = (rs5/r4)*((Zr-b4)/(Cm5*r5*b5))
    Ka = 0.4
    Kr = 0.75
    Kar = -0.3
    Ea = 0.0003
    Er = 0.0003
    LossTip = 1e-3*(U4**3*NR/(8*np.pi))*(Ka*Ea*Ca+Kr*Er*Cr+Kar*np.sqrt(Ea*Er*Ca*Cr))

    #Trailing Edge Losses
    tb4 = 0.04*r4
    tb5 = 0.02*r4
    LossTE = 1e-3*rho5*W5**2/2*(NR*tb5/(np.pi*(rh5+rs5)*np.cos(Beta5*np.pi/180)))**2

    #Sum Losses
    Total_Loss = LossInc+LossTE+LossTip

    Reaction = (H_4-h5s)/(h4s-h5s)
    s5 = 2*np.pi*r5/NR
    o5 = Cm5*s5/W5

    ha1=Pr('H','T',T_1,'P',P_1*1e6,fluid)/1000#+0.5*(C4**2)/1000
    ha3=Pr('H','T',T_5,'P',P_5*1e6,fluid)/1000#+0.5*(C5**2)/1000
    saa1=Pr('S','H',ha1*1000,'P',P_1*1e6,fluid)
    ha3ss=Pr('H','S',saa1,'P',P_5*1e6,fluid)/1000-0.5*(C5**2)/1000
    Effx = (ha1-ha3)/(ha1-ha3ss)- Total_Loss/(ha1-ha3ss)

    return Effx,Reaction,NR,r4,Alpha4,b4,Ct4,rho4,mflow,H_1,T_1,P_1,Zr,rs5,rh5,tb4,tb5

#Create 2D Quasi
def QuasiNorm(b4,r4,Zr,rs5,rh5,n,Z5):
    global Z,r,Zh,rh,m,Ash,Bsh,Csh,Dsh,Esh,Fsh
    Betha4 = np.arctan(Cm4/(Ct4-U4))
    QuasiSec=50
    dZ=Zr/QuasiSec
    C2=(r4-rs5)/((dZ-b4)**n)
    Zrb4=abs(b4-Zr)
    r5=(rs5+rh5)/2
    b5=(rs5-rh5)
    #Shroud Sections
    Z=[]
    Z.append(0)
    m=[]
    m.append(0)
    for i in range(0,QuasiSec-1):
        Z.append(Z[i]+dZ)
    Z=np.array(Z)
    dzb4=dZ-b4
    for i in [Z]:
        Epsi = i/dzb4
        r=rs5+(r4-rs5)*Epsi**n
    for i in [Z]:
        fz=np.sqrt(1+(C2*n*(i-(i-1))**(n-1))**2)
    m.append((dZ/3)*(fz[0]+4*fz[1]) )
    for i in range (2, len(fz)):
        if (i % 2) == 0:
            m.append(m[i-1]+4*(fz[i])/3)
        else:
            m.append(m[i-1]+2*(fz[i])/3)                            #perbaiki value fz #check OK

    #Hub Sections
    Zh=[]
    Zh.append(0)
    mh=[]
    mh.append(0)
    rh=[]
    fzh=[]
    for i in range(1,QuasiSec):
        Zh.append(Zh[i-1]+dZ)
        Ephi = Zh[i]/(dZ-b4)
    for i in range(0,QuasiSec):
        rh.append(rh5+Zr-np.sqrt((Zr**2)-(Zh[i]**2)))
        fzh.append(np.sqrt(1+(C2*n*(Zh[i]-Z5)**(n-1))**2))
    mh.append((dZ/3)*(fzh[0]+4)/fzh[0])                      #Perbaiki
    for i in range(1,len(fzh)):
        mh.append(mh[i]+2*(i)/3)                                 #Perbaiki value mh 


    # Zm=(Z-Zh)/2
    # rm=(r-rh)/2
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
    return(Z,r,Zh,rh,m,Ash,Bsh,Csh,Dsh,Esh,Fsh)

#Perhitungan Sudut   
def Zrregs(Z,r,Zh,rh,m,Ash,Bsh,Csh,Dsh,Esh,Fsh):
    global phis,tethas,Bethacs,phih,tethah,Bethach
    #Shroud
    dms=[]
    dms.append(0)
    for i in range(1,len(Z)):
        dms=np.sqrt((Z[i]-Z[i-1])**2+(r[i]-r[i-1])**2)
    msi=[]
    msi.append(0)
    for i in range(1,len(Z)):
        msi.append(msi[i-1]+dms)
    polynom = 6             #Coeficient harus 6 (Belum otomatis)
    mymodel = np.poly1d(np.polyfit(Z,r,polynom))
    #need to regress the mymodel
    #Define the Function of phi
    coeffphis=[]
    phis=[]
    tethas=[]
    Bethacs=[]
    for coeff in mymodel:
        coeffphis.append(coeff)
    for i in range(0,len(msi)):
        phis.append(6*coeffphis[5]*msi[i]**5+5*coeffphis[4]*msi[i]**4+4*coeffphis[3]*msi[i]**3+3*coeffphis[2]*msi[i]**2+2*coeffphis[1]*msi[i]**1+coeffphis[0])
        tethas.append(msi[i]*(Ash+Bsh**3+Csh**4))
        Bethacs.append(1/np.arctan(r[i]*(Ash+3*Bsh**2+4*Csh**3)))

    #Hub
    dmh=[]
    dmh.append(0)
    for i in range(1,len(Zh)):
        dmh.append(np.sqrt((Zh[i]-Zh[i-1])**2+(rh[i]-rh[i-1])**2))
    mhi=[]
    mhi.append(0)
    for i in range(1,len(dmh)):
        mhi.append(mhi[i-1]+dmh[i])
    polynomh = 6
    mymodelh = np.poly1d(np.polyfit(Zh,rh,polynomh))
    coeffphih=[]
    phih=[]
    tethah=[]
    Bethach=[]
    for coeffh in mymodelh:
        coeffphih.append(coeffh)
    for i in range(0,len(mhi)):
        phih.append(6*coeffphih[5]*mhi[i]**5+5*coeffphih[4]*mhi[i]**4+4*coeffphih[3]*mhi[i]**3+3*coeffphih[2]*mhi[i]**2+2*coeffphih[1]*mhi[i]**1+coeffphih[0])
        tethah.append(1/np.tan(rh[i]*(mhi[i]*(Dsh+Esh**2+Fsh**3))))
        Bethach.append(1/np.arctan(rh[i]*(Dsh+3*Esh**2+4*Fsh**3)))

    #Return Results
    return(phis,tethas,Bethacs,phih,tethah,Bethach)

#Meridional Coordinate
def meridional(phis,tethas,Bethacs,phih,tethah,Bethach):
    global X,Y,Z,Xh,Yh,Zh
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
    return(X,Y,Z,Xh,Yh,Zh)

#Vector Component
def VectorComp(X,Y,Z,Xh,Yh,Zh):
    global Txs,Tys,Tzs,Txh,Tyh,Tzh
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
    for i in range(0,len(X)):
        L.append(np.sqrt((X[i]-Xh[i])**2+(Y[i]-Yh[i])**2+(Z[i]-Zh[i])**2))

    #Vector B Hub and Shroud
    
    for i in range(0,len(X)):
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

    return(Txs,Tys,Tzs,Txh,Tyh,Tzh,tb4,tb5)

#3D Coordinate

def Coord3D(Txs,Tys,Tzs,Txh,Tyh,Tzh,tb4,tb5):
    #Tebal Sudu
    tb =[]
    tb.append(tb4)
    for i in range(1,len(Txs)):
        tb.append(tb[i-1]+(abs(tb4-tb5)/len(Txs)))

    #Shroud 
    XcorS=[]
    YcorS=[]
    ZcorS=[]
    for i in range(0,len(Txs)):
        XcorS.append(X[i]+0.5*Txs[i]*tb[i])
        XcorS.append(X[i]-0.5*Txs[i]*tb[i])
        YcorS.append(Y[i]+0.5*Tys[i]*tb[i])
        YcorS.append(Y[i]-0.5*Tys[i]*tb[i])
        ZcorS.append(Z[i]+0.5*Tzs[i]*tb[i])
        ZcorS.append(Z[i]-0.5*Tzs[i]*tb[i])

    #Hub
    XcorH=[]
    YcorH=[]
    ZcorH=[]
    for i in range(0,len(Txs)):
        XcorH.append(Xh[i]+0.5*Txh[i]*tb[i])
        XcorH.append(Xh[i]-0.5*Txh[i]*tb[i])
        YcorH.append(Yh[i]+0.5*Tyh[i]*tb[i])
        YcorH.append(Yh[i]-0.5*Tyh[i]*tb[i])
        ZcorH.append(Zh[i]+0.5*Tzh[i]*tb[i])
        ZcorH.append(Zh[i]-0.5*Tzh[i]*tb[i])
    return (XcorS,YcorS,ZcorS,XcorH,YcorH,ZcorH)