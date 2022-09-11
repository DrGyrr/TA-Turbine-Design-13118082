
from itertools import combinations_with_replacement
from ifrturbinepackage.definitions import *
from ifrturbinepackage.inputs import *
from ifrturbinepackage.rotor import *


k = 10#cycle

def ComputeR3t(x):
    try:
        global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
        global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
        global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
        global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
        global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
        global Effts,Efftt,Efftspred,Reaction,vNu
        tenflow_coeff = x[0]
        tenwork_coeff = x[1]
        Rr5r4         = x[2]
        Rb5b4         = x[3]
        Rb4r4         = x[4]
        RZrr4         = x[5]
        NR            = x[6]
        rpm           = x[7]
        flow_coeff=tenflow_coeff/10
        work_coeff=tenwork_coeff/10

        cycledict=whichcycle (k)       # The cycle to be computed
        globals().update(cycledict)
        # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
        # globals().update(gparamdict)
        # rpm =whatrpm(m)          # rpm at m will be used

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
        # if rh5 < 0.0015:
        #     print(f"For flow coeff ={flow_coeff} and work coeff={work_coeff} rh5 too small(<1.5mm), adjust gparams")
        #     return
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
        if 0.04*r4>0.001:
            tb4 = 0.04*r4
        else:
            tb4 = 0.001
        if 0.02*r4>0.001:
            tb5 = 0.02*r4
        else:
            tb5 = 0.001
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

        return Effts*-1
    except ValueError:
        return 0

# tenflow_coeff = x[0]
# tenwork_coeff = x[1]
# Rr5r4         = x[2]
# Rb5b4         = x[3]
# Rb4r4         = x[4]
# RZrr4         = x[5]
# NR            = x[6]
# rpm           = x[7]

def hconstr(x):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    tenflow_coeff = x[0]
    tenwork_coeff = x[1]
    Rr5r4         = x[2]
    Rb5b4         = x[3]
    Rb4r4         = x[4]
    RZrr4         = x[5]
    NR            = x[6]
    rpm           = x[7]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    cycledict=whichcycle (k)       # The cycle to be computed
    globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
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
    return (h04s-h4s)/h04s-0.2

def sonicconstr(x):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    tenflow_coeff = x[0]
    tenwork_coeff = x[1]
    Rr5r4         = x[2]
    Rb5b4         = x[3]
    Rb4r4         = x[4]
    RZrr4         = x[5]
    NR            = x[6]
    rpm           = x[7]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    cycledict=whichcycle (k)       # The cycle to be computed
    globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used

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
    a4      = Props('A','H',h4s,'S',s04s,fluid)
    return (a4-Cm4)/a4-0.025

def rh5constr(x):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    tenflow_coeff = x[0]
    tenwork_coeff = x[1]
    Rr5r4         = x[2]
    Rb5b4         = x[3]
    Rb4r4         = x[4]
    RZrr4         = x[5]
    NR            = x[6]
    rpm           = x[7]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    cycledict=whichcycle (k)       # The cycle to be computed
    globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
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

    r4      = U4/np.radians(rpm*6)
    r5      = Rr5r4*r4
    b4      = Rb4r4*r4
    b5      = Rb5b4*b4
    rs5     = (2*r5+b5)/2
    rh5     = rs5-b5
    return rh5-0.003

def rs5constr(x):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    tenflow_coeff = x[0]
    tenwork_coeff = x[1]
    Rr5r4         = x[2]
    Rb5b4         = x[3]
    Rb4r4         = x[4]
    RZrr4         = x[5]
    NR            = x[6]
    rpm           = x[7]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    cycledict=whichcycle (k)       # The cycle to be computed
    globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
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

    r4      = U4/np.radians(rpm*6)
    r5      = Rr5r4*r4
    b4      = Rb4r4*r4
    b5      = Rb5b4*b4
    rs5     = (2*r5+b5)/2
    rh5     = rs5-b5
    return (r4-rs5)/r4-0.3

# def mflowconstr(x):
#     global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
#     global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
#     global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
#     global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
#     global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
#     global Effts,Efftt,Efftspred,Reaction,vNu
#     tenflow_coeff = x[0]
#     tenwork_coeff = x[1]
#     Rr5r4         = x[2]
#     Rb5b4         = x[3]
#     Rb4r4         = x[4]
#     RZrr4         = x[5]
#     NR            = x[6]
#     rpm           = x[7]
#     flow_coeff=tenflow_coeff/10
#     work_coeff=tenwork_coeff/10
#     cycledict=whichcycle (k)       # The cycle to be computed
#     globals().update(cycledict)
#     # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
#     # globals().update(gparamdict)
#     # rpm =whatrpm(m)          # rpm at m will be used
#     Cp4 = Props('C','T',T_1,'P',P_1,fluid)
#     Cv4 = Props('O','T',T_1,'P',P_1,fluid)
#     gamma = Cp4/Cv4
#     Rx = 8.31446261815324   #J/K.mol
#     #General Properties inlet outlet turbin (Total)
#     H_1     = Props('H','T',T_1,'P',P_1,fluid)     #J/kg
#     s01     = Props('S','T',T_1,'P',P_1,fluid)     #J/kg.K 
#     T_5     = Props('T','P',P_5,'S',s01,fluid)  # =>asumsi nozzle isenthalpy DAN Isentropic
#     H_5     = Props('H','T',T_5,'P',P_5,fluid)  # meski pada kenyataannya isenthalpic nozzle tidak isentropic
#     DeltaH  = H_1-H_5            #Ideal === Isentropic Total Enthalpy change 
#     C0s     = np.sqrt(2*DeltaH)         #Spouting Velocity
#     #Perhitungan Properties ideal lain (Total)
#     p01     = P_1           #inlet volute [1], Total
#     T01     = T_1
#     h01     = H_1
#     p1      = p01           # inlet turbine, V~0 
#     T1      = T_1
#     h01     = H_1
#     rho1   = Props('D','P',p1,'T',T1,fluid)
#     h02s    = H_1           #inlet nozzle [2], Total
#     s02s    = s01            #ideal volute === approx. as isentropic
#     p02s    = p01
#     T02s    = T01
#     h03s    = h02s           #outlet nozzle [3], Total
#     s03s    = s02s            #ideal nozzle === approx. as isentropic (in Total)
#     p03s    = p02s
#     T03s    = T02s
#     h04s    = h03s           #inlet rotor [4], Total
#     s04s    = s03s           #outlet nozzle === inlet rotor
#     p04s    = p03s
#     T04s    = T03s
#     h04     = h04s          # Nozzle isenthalpic but not isentropic
#     p05ss   = P_5
#     T05ss   = T_5
#     h05ss   = H_5
#     s05ss   = s04s
#     #Segitiga Kecepatan Inlet, m/s, radians
#     U4      = np.sqrt(DeltaH/work_coeff)
#     Cm4     = U4*flow_coeff
#     Ct4     = DeltaH/U4                 # => DeltaH = U4*Ct4-U5*Ct5 ; Alpha5=0 => Ct5=0
#     C4      = np.sqrt(Cm4**2+Ct4**2)
#     Alpha4  = np.arctan(Ct4/Cm4)
#     W4      = np.sqrt(Cm4**2+(U4-Ct4)**2)
#     Beta4   = np.arctan((U4-Ct4)/Cm4)
#     h4s    = h04s-1/2*Cm4**2
#     rho4s   = Props('D','H',h4s,'S',s04s,fluid)
#     rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)


#     #Perhitungan geometri
#     r4      = U4/np.radians(rpm*6)
#     r5      = Rr5r4*r4
#     b4      = Rb4r4*r4
#     b5      = Rb5b4*b4
#     rs5     = (2*r5+b5)/2
#     rh5     = rs5-b5
#     Zr      = RZrr4*r4
#     mflow   = 2*np.pi*b4*r4*rho4s*Cm4
#     return mflow-0.4

def workconstr(x):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    tenflow_coeff = x[0]
    tenwork_coeff = x[1]
    Rr5r4         = x[2]
    Rb5b4         = x[3]
    Rb4r4         = x[4]
    RZrr4         = x[5]
    NR            = x[6]
    rpm           = x[7]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    cycledict=whichcycle (k)       # The cycle to be computed
    globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
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
    # if rh5 < 0.0015:
    #     print(f"For flow coeff ={flow_coeff} and work coeff={work_coeff} rh5 too small(<1.5mm), adjust gparams")
    #     return
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
    T05 = Props('T','H',h05,'P',p05,fluid)
    p5  = p5ss
    T5  = Props('T','H',h5,'P',p5,fluid)
    return (h01-h05)*mflow-5*10**3
def r4constr(x):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    tenflow_coeff = x[0]
    tenwork_coeff = x[1]
    Rr5r4         = x[2]
    Rb5b4         = x[3]
    Rb4r4         = x[4]
    RZrr4         = x[5]
    NR            = x[6]
    rpm           = x[7]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    cycledict=whichcycle (k)       # The cycle to be computed
    globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used

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

    r4      = U4/np.radians(rpm*6)
    r5      = Rr5r4*r4
    b4      = Rb4r4*r4
    b5      = Rb5b4*b4
    rs5     = (2*r5+b5)/2
    rh5     = rs5-b5
    return 0.1-r4

#  # def rhoconstr(x):

#     global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
#     global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
#     global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
#     global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
#     global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
#     global Effts,Efftt,Efftspred,Reaction,vNu
#     tenflow_coeff = x[0]
#     tenwork_coeff = x[1]
#     Rr5r4         = x[2]
#     Rb5b4         = x[3]
#     Rb4r4         = x[4]
#     RZrr4         = x[5]
#     NR            = x[6]
#     rpm           = x[7]
#     flow_coeff=tenflow_coeff/10
#     work_coeff=tenwork_coeff/10
#     cycledict=whichcycle (k)       # The cycle to be computed
#     globals().update(cycledict)
#     # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
#     # globals().update(gparamdict)
#     # rpm =whatrpm(m)          # rpm at m will be used
#     Cp4 = Props('C','T',T_1,'P',P_1,fluid)
#     Cv4 = Props('O','T',T_1,'P',P_1,fluid)
#     gamma = Cp4/Cv4
#     Rx = 8.31446261815324   #J/K.mol
#     #General Properties inlet outlet turbin (Total)
#     H_1     = Props('H','T',T_1,'P',P_1,fluid)     #J/kg
#     s01     = Props('S','T',T_1,'P',P_1,fluid)     #J/kg.K 
#     T_5     = Props('T','P',P_5,'S',s01,fluid)  # =>asumsi nozzle isenthalpy DAN Isentropic
#     H_5     = Props('H','T',T_5,'P',P_5,fluid)  # meski pada kenyataannya isenthalpic nozzle tidak isentropic
#     DeltaH  = H_1-H_5            #Ideal === Isentropic Total Enthalpy change 
#     C0s     = np.sqrt(2*DeltaH)         #Spouting Velocity
#     #Perhitungan Properties ideal lain (Total)
#     p01     = P_1           #inlet volute [1], Total
#     T01     = T_1
#     h01     = H_1
#     p1      = p01           # inlet turbine, V~0 
#     T1      = T_1
#     h01     = H_1
#     rho1   = Props('D','P',p1,'T',T1,fluid)
#     h02s    = H_1           #inlet nozzle [2], Total
#     s02s    = s01            #ideal volute === approx. as isentropic
#     p02s    = p01
#     T02s    = T01
#     h03s    = h02s           #outlet nozzle [3], Total
#     s03s    = s02s            #ideal nozzle === approx. as isentropic (in Total)
#     p03s    = p02s
#     T03s    = T02s
#     h04s    = h03s           #inlet rotor [4], Total
#     s04s    = s03s           #outlet nozzle === inlet rotor
#     p04s    = p03s
#     T04s    = T03s
#     h04     = h04s          # Nozzle isenthalpic but not isentropic
#     p05ss   = P_5
#     T05ss   = T_5
#     h05ss   = H_5
#     s05ss   = s04s
#     #Segitiga Kecepatan Inlet, m/s, radians
#     U4      = np.sqrt(DeltaH/work_coeff)
#     Cm4     = U4*flow_coeff
#     Ct4     = DeltaH/U4                 # => DeltaH = U4*Ct4-U5*Ct5 ; Alpha5=0 => Ct5=0
#     C4      = np.sqrt(Cm4**2+Ct4**2)
#     Alpha4  = np.arctan(Ct4/Cm4)
#     W4      = np.sqrt(Cm4**2+(U4-Ct4)**2)
#     Beta4   = np.arctan((U4-Ct4)/Cm4)
#     h4s    = h04s-1/2*Cm4**2
#     rho4s   = Props('D','H',h4s,'S',s04s,fluid)
#     rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)

def cm5loop1constr(x):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    tenflow_coeff = x[0]
    tenwork_coeff = x[1]
    Rr5r4         = x[2]
    Rb5b4         = x[3]
    Rb4r4         = x[4]
    RZrr4         = x[5]
    NR            = x[6]
    rpm           = x[7]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10

    cycledict=whichcycle (k)       # The cycle to be computed
    globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
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

    # plimit = Props('P_min',fluid)
    T5sslimit = 273 + 30
    h5sslimit=Props('H','T',T5sslimit,'S',s05ss,fluid)
    # rho5sslimit=Props('D','T',T5sslimit,'S',s05ss,fluid)

 
    while Cm5didconverge1 == False:
        k1Cm5       = k1Cm5+1             # => iteration amount
        Cm5i        = Cm5ii
        rho5ssi      = rho5ssii
        Cm5ii       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ssi)*Cm4
        # Cm4ii       = mflow/(2*np.pi()*b5*)
        h5ss         = h05ss-1/2*(Cm5ii**2+Ct5**2)
        if h5sslimit>h5ss:
            break
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

    return h5ss-h5sslimit

def cm5loop2constr(x):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    tenflow_coeff = x[0]
    tenwork_coeff = x[1]
    Rr5r4         = x[2]
    Rb5b4         = x[3]
    Rb4r4         = x[4]
    RZrr4         = x[5]
    NR            = x[6]
    rpm           = x[7]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    cycledict=whichcycle (k)       # The cycle to be computed
    globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
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

    T5sslimit = 273 + 30
    h5sslimit=Props('H','T',T5sslimit,'S',s05ss,fluid)

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
        if h5sslimit>h5ss:
            break
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
    return h5ss-h5sslimit

def twophaseconstr(x):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    tenflow_coeff = x[0]
    tenwork_coeff = x[1]
    Rr5r4         = x[2]
    Rb5b4         = x[3]
    Rb4r4         = x[4]
    RZrr4         = x[5]
    NR            = x[6]
    rpm           = x[7]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    cycledict=whichcycle (k)       # The cycle to be computed
    globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
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
    satfors=getsatpvfors(P_5,s05ss,fluid)
    satpv=satfors['P']+satfors['acc']
    hsatv = Props('H','P',satpv,'S',s05ss,fluid)


    return (h5ss-hsatv)/hsatv-0.05

tenflowb    =(1,15)
tenworkb    =(10,25)
Rr5r4b      =(0.15,0.4)
Rb5b4b      =(1,2.5)
Rb4r4b      =(0.08,0.2)
RZrr4b      =(0.5,2)
NRb         =(8,20)
rpmb        =(15000,50000)
bnds        =(tenflowb,tenworkb,Rr5r4b,Rb5b4b,Rb4r4b,RZrr4b,NRb,rpmb)
constr1      ={'type': 'ineq', 'fun':constraint1}
constr2      ={'type': 'ineq', 'fun':constraint1}
constr3      ={'type': 'ineq', 'fun':constraint1}
constr4      ={'type': 'ineq', 'fun':constraint1}
constr5      ={'type': 'ineq', 'fun':constraint1}
constr6      ={'type': 'ineq', 'fun':constraint1}
constr7      ={'type': 'ineq', 'fun':constraint1}
constr8      ={'type': 'ineq', 'fun':constraint1}
constrh     ={'type': 'ineq','fun':hconstr}
constrsonic ={'type': 'ineq','fun':sonicconstr}
constrrh5 ={'type': 'ineq','fun':rh5constr}
constrrs5 ={'type': 'ineq','fun':rs5constr}
constrwork ={'type':'ineq','fun':workconstr}
constrr4={'type':'ineq','fun':r4constr}
# constrmflow={'type': 'ineq','fun':mflowconstr}
constrcm5loop1={'type': 'ineq','fun':cm5loop1constr}
constrcm5loop2={'type': 'ineq','fun':cm5loop2constr}
constrtwophase={'type':'ineq','fun':twophaseconstr}

constr       = [constr1,constr2,constr3,constr4,constr5,constr6,constr7,constr8,constrsonic,constrrh5,constrrs5,constrh,constrwork,constrr4,constrcm5loop1,constrcm5loop2,constrtwophase]
initval      = [1.2,10,0.35,2,0.15,1,10,50000]

opteffts_tesuto = optimize.minimize(ComputeR3t,initval,method='SLSQP',bounds=bnds,constraints=constr)
print(opteffts_tesuto)