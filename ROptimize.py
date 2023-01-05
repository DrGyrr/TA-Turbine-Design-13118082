from time import sleep
from ifrturbinepackage.definitions  import *
from ifrturbinepackage.inputs       import *
from ifrturbinepackage.rotor        import *



def ComputeR4varg(x:list) -> 'Effts':

    try:
        
        global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
        global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
        global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
        global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
        global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
        global Effts,Efftt,Efftspred,Reaction,vNu

        Rr5r4   = x[0]
        Rb5b4   = x[1]
        Rb4r4   = x[2]
        RZrr4   = x[3]
        NR      = x[4]

        flow_coeff=tenflow_coeff/10
        work_coeff=tenwork_coeff/10

        # cycledict=whichcycle (cyclenum)       # The cycle to be computed
        # globals().update(cycledict)
        
        # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
        # globals().update(gparamdict)
        # rpm =whatrpm(m)          # rpm at m will be used
        # P_1 = P_1*10**6   sudah diubah jadi pa di fungsi whichcycle
        # P_5 = P_5*10**6

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

        h4s    = h04s-1/2*C4**2
        rho4s   = Props('D','H',h4s,'S',s04s,fluid)
        rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)

        initr4guess=0.045 # => initial guess r4 = 4.5cm
        r40     = initr4guess
        r4      = 0.04  # => set asal untuk memulai loop
        while np.abs(r40-r4)/r40 > 0.005: # residual/error harus lebih kecil dari 0.5%
            r40     = r4
            if 0.04*r4>0.001:
                tb4 = 0.04*r4
            else:
                tb4 = 0.001
            Bk4     = (tb4*NR*0.05)/(2*np.pi*r40*np.cos(Beta4))
            r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s*(1-Bk4))) # mflow sebagai input

        angvel  = U4/r4
        rpm     = angvel*(60/(2*np.pi))
        r5      = Rr5r4*r4
        b4      = Rb4r4*r4
        b5      = Rb5b4*b4
        rs5     = (2*r5+b5)/2
        rh5     = rs5-b5

        if 0.02*r4>0.001:
            tb5 = 0.02*r4
        else:
            tb5 = 0.001

        Ct5 = 0 # => it is predetermined that Alpha5=0
        Alpha5 = 0
        Cm5_0    = 10
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

            C5      = np.sqrt(Cm5i**2+Ct5**2)
            U5      = U4*Rr5r4
            W5      = np.sqrt(Cm5i**2+(U5-Ct5)**2)
            Beta5   = np.arccos(Cm5i/W5)

            Bk5     = (tb5*0.05*NR)/(2*np.pi*r5*np.cos(Beta5))
            Cm5ii       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ssi)*(1-Bk4)/(1-Bk5) *Cm4
            # Cm4ii       = mflow/(2*np.pi()*b5*)
            h5ss         = h05ss-1/2*(Cm5ii**2+Ct5**2)
            rho5ssii     = Props('D','H',h5ss,'S',s05ss,fluid)
            errorCm5    = np.abs((Rb5b4*Rr5r4*(rho5ssii/rho4s)*(Cm5ii/Cm4)*((1-Bk5)/(1-Bk4)))-1)
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

            C5      = np.sqrt(Cm5**2+Ct5**2)
            U5      = U4*Rr5r4
            W5       = np.sqrt(Cm5**2+(U5-Ct5)**2)
            Beta5   = np.arccos(Cm5/W5)

            Bk5     = (tb5*0.05*NR)/(2*np.pi*r5*np.cos(Beta5))

            Cm5       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ss)*(1-Bk4)/(1-Bk5) *Cm4 # -_- -_- -_-
            h5ss       = h05ss-1/2*(Cm5**2+Ct5**2)
            rho5ss     = Props('D','H',h5ss,'S',s05ss,fluid)
            if np.abs(1-Cm5/Props('A','H',h5ss,'S',s05ss,fluid)) < 5*1e-3:
                choked5 = True
                break
            errorCm5  = np.abs((Rb5b4*Rr5r4*(rho5ss/rho4s)*(Cm5/Cm4)*((1-Bk5)/(1-Bk4)))-1)
            if errorCm5 <= 5*1e-5:
                Cm5didconverge2 = True
                break
            if k2Cm5>200:
                print(f"loop2 iterates too long ({k1Cm5},{k2Cm5}) at {flow_coeff,work_coeff} with errorCm5 = {errorCm5}")
                break
        h5ss    = h05ss-1/2*(Cm5**2+Ct5**2)
        C5      = np.sqrt(Cm5**2+Ct5**2)
        U5      = U4*Rr5r4
        W5       = np.sqrt(Cm5**2+(U5-Ct5)**2)
        Beta5   = np.arccos(Cm5/W5)



        #Perhitungan geometri
        # r4      = U4/np.radians(rpm*6)
        

        # if rh5 < 0.0015:
        #     print(f"For flow coeff ={flow_coeff} and work coeff={work_coeff} rh5 too small(<1.5mm), adjust gparams")
        
        Zr      = RZrr4*r4

        # mflow   = 2*np.pi*b5*r5*rho5ss*Cm5

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
        tesdict     = dict()

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
        for i in ('Beta4','Beta5','b4','r4','Zr','rs5','rh5'):
            proceeddict[i]  = globals()[i]
        tesdict['rpm'] = rpm

        outputdict  = {
            'geometry'  : geomdict,
            'thermo'    : thermodict,
            'velocity'  : veltridict,
            'efficiency': effdict,
            'losses'    : lossdict,
            'proceed'   : proceeddict,
            'tes'       : tesdict
            }

        return Effts*-1
    except ValueError:
        return 0


def ComputeR4vargd(x:list) -> dict :

    try:
        
        global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
        global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
        global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
        global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
        global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
        global Effts,Efftt,Efftspred,Reaction,vNu

        Rr5r4   = x[0]
        Rb5b4   = x[1]
        Rb4r4   = x[2]
        RZrr4   = x[3]
        NR      = x[4]

        flow_coeff=tenflow_coeff/10
        work_coeff=tenwork_coeff/10

        cycledict=whichcycle (cyclenum)       # The cycle to be computed
        globals().update(cycledict)
        # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
        # globals().update(gparamdict)
        # rpm =whatrpm(m)          # rpm at m will be used
        # P_1 = P_1*10**6   sudah diubah jadi pa di fungsi whichcycle
        # P_5 = P_5*10**6

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

        h4s    = h04s-1/2*C4**2
        rho4s   = Props('D','H',h4s,'S',s04s,fluid)
        rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)

        initr4guess=0.045 # => initial guess r4 = 4.5cm
        r40     = initr4guess
        r4      = 0.04  # => set asal untuk memulai loop
        while np.abs(r40-r4)/r40 > 0.01: # residual/error harus lebih kecil dari 1%
            r40     = r4
            if 0.04*r4>0.001:
                tb4 = 0.04*r4
            else:
                tb4 = 0.001
            Bk4     = (tb4*0.05*NR)/(2*np.pi*r40*np.cos(Beta4))
            r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s*(1-Bk4))) # mflow sebagai input

        angvel  = U4/r4
        rpm     = angvel*(60/(2*np.pi))
        r5      = Rr5r4*r4
        b4      = Rb4r4*r4
        b5      = Rb5b4*b4
        rs5     = (2*r5+b5)/2
        rh5     = rs5-b5

        if 0.02*r4>0.001:
            tb5 = 0.02*r4
        else:
            tb5 = 0.001

        Ct5 = 0 # => it is predetermined that Alpha5=0
        Alpha5 = 0
        Cm5_0    = 10
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

            C5      = np.sqrt(Cm5i**2+Ct5**2)
            U5      = U4*Rr5r4
            W5      = np.sqrt(Cm5i**2+(U5-Ct5)**2)
            Beta5   = np.arccos(Cm5i/W5)

            Bk5     = (tb5*0.05*NR)/(2*np.pi*r5*np.cos(Beta5))
            Cm5ii       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ssi)*(1-Bk4)/(1-Bk5) *Cm4
            # Cm4ii       = mflow/(2*np.pi()*b5*)
            h5ss         = h05ss-1/2*(Cm5ii**2+Ct5**2)
            rho5ssii     = Props('D','H',h5ss,'S',s05ss,fluid)
            errorCm5    = np.abs((Rb5b4*Rr5r4*(rho5ssii/rho4s)*(Cm5ii/Cm4)*((1-Bk5)/(1-Bk4)))-1)
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
            
            C5      = np.sqrt(Cm5**2+Ct5**2)
            U5      = U4*Rr5r4
            W5       = np.sqrt(Cm5**2+(U5-Ct5)**2)
            Beta5   = np.arccos(Cm5/W5)

            Bk5     = (tb5*0.05*NR)/(2*np.pi*r5*np.cos(Beta5))

            Cm5       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ss)*(1-Bk4)/(1-Bk5) *Cm4 # -_- -_- -_-
            h5ss       = h05ss-1/2*(Cm5**2+Ct5**2)
            rho5ss     = Props('D','H',h5ss,'S',s05ss,fluid)
            if np.abs(1-Cm5/Props('A','H',h5ss,'S',s05ss,fluid)) < 5*1e-3:
                choked5 = True
                break
            errorCm5  = np.abs((Rb5b4*Rr5r4*(rho5ss/rho4s)*(Cm5/Cm4)*((1-Bk5)/(1-Bk4)))-1)
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
        # r4      = U4/np.radians(rpm*6)
        # r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s)) # mflow sebagai input
        # angvel  = U4/r4
        # rpm     = angvel*(60/(2*np.pi))
        # r5      = Rr5r4*r4
        # b4      = Rb4r4*r4
        # b5      = Rb5b4*b4
        # rs5     = (2*r5+b5)/2
        # rh5     = rs5-b5
        # if rh5 < 0.0015:
        #     print(f"For flow coeff ={flow_coeff} and work coeff={work_coeff} rh5 too small(<1.5mm), adjust gparams")
        
        Zr      = RZrr4*r4

        # mflow   = 2*np.pi*b5*r5*rho5ss*Cm5

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
        tesdict     = dict()

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
        for i in ('Beta4','Beta5','b4','r4','Zr','rs5','rh5'):
            proceeddict[i]  = globals()[i]
        tesdict['rpm'] = rpm

        outputdict  = {
            'geometry'  : geomdict,
            'thermo'    : thermodict,
            'velocity'  : veltridict,
            'efficiency': effdict,
            'losses'    : lossdict,
            'proceed'   : proceeddict,
            'tes'       : tesdict
            }

        return outputdict
    except ValueError:
        return 0


def rh5constrvarg(x:list):

    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    Rr5r4   = x[0]
    Rb5b4   = x[1]
    Rb4r4   = x[2]
    RZrr4   = x[3]
    NR      = x[4]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    # cycledict=whichcycle (cyclenum)       # The cycle to be computed
    # globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
    # P_1 = P_1*10**6 sudah diubah jadi pa di fungsi whichcycle
    # P_5 = P_5*10**6
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

    h4s    = h04s-1/2*C4**2
    rho4s   = Props('D','H',h4s,'S',s04s,fluid)
    rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)

    initr4guess=0.045 # => initial guess r4 = 4.5cm
    r40     = initr4guess
    r4      = 0.04  # => set asal untuk memulai loop
    while np.abs(r40-r4)/r40 > 0.01: # residual/error harus lebih kecil dari 1%
        r40     = r4
        if 0.04*r4>0.001:
            tb4 = 0.04*r4
        else:
            tb4 = 0.001
        Bk4     = (tb4*0.05*NR)/(2*np.pi*r40*np.cos(Beta4))
        r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s*(1-Bk4))) # mflow sebagai input
    angvel  = U4/r4
    rpm     = angvel*(60/(2*np.pi))
    r5      = Rr5r4*r4
    b4      = Rb4r4*r4
    b5      = Rb5b4*b4
    rs5     = (2*r5+b5)/2
    rh5     = rs5-b5
    if 0.02*r4>0.001:
        tb5 = 0.02*r4
    else:
        tb5 = 0.001


    return rh5-0.0075

def rs5constrvarg(x:list):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    Rr5r4   = x[0]
    Rb5b4   = x[1]
    Rb4r4   = x[2]
    RZrr4   = x[3]
    NR      = x[4]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    # cycledict=whichcycle (cyclenum)       # The cycle to be computed
    # globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
    # P_1 = P_1*10**6   sudah diubah jadi pa di fungsi whichcycle
    # P_5 = P_5*10**6
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

    h4s    = h04s-1/2*C4**2
    rho4s   = Props('D','H',h4s,'S',s04s,fluid)
    rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)

    initr4guess=0.045 # => initial guess r4 = 4.5cm
    r40     = initr4guess
    r4      = 0.04  # => set asal untuk memulai loop
    while np.abs(r40-r4)/r40 > 0.01: # residual/error harus lebih kecil dari 1%
        r40     = r4
        if 0.04*r4>0.001:
            tb4 = 0.04*r4
        else:
            tb4 = 0.001
        Bk4     = (tb4*0.05*NR)/(2*np.pi*r40*np.cos(Beta4))
        r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s*(1-Bk4))) # mflow sebagai input
    angvel  = U4/r4
    rpm     = angvel*(60/(2*np.pi))
    r5      = Rr5r4*r4
    b4      = Rb4r4*r4
    b5      = Rb5b4*b4
    rs5     = (2*r5+b5)/2
    rh5     = rs5-b5
    if 0.02*r4>0.001:
        tb5 = 0.02*r4
    else:
        tb5 = 0.001

    return (r4-rs5)/r4-0.2

def GetGParams(coeffs,mflows): #tenflow_coeff,tenwork_coeff
    global tenflow_coeff,tenwork_coeff,mflow

    tenflow_coeff = coeffs [0]
    tenwork_coeff = coeffs [1]
    mflow    = mflows
    
    Rr5r4b   = (0.15,0.4)
    Rb5b4b   = (1,10)
    Rb4r4b   = (0.09,0.3)
    RZrr4b   = (0.75,1.4)
    NRb      = (10,17)

    bnds        = (Rr5r4b,Rb5b4b,Rb4r4b,RZrr4b,NRb)
    constrrh5   = {'type': 'ineq', 'fun':rh5constrvarg}
    constrrs5   = {'type': 'ineq', 'fun':rs5constrvarg}

    constr      = [constrrh5,constrrs5]
    initval     = [0.15, 1.6, 0.3, 0.8, 12]

    opteffts_geo    = optimize.minimize(ComputeR4varg,initval,method='SLSQP',bounds=bnds,constraints=constr)

    outputdict = {'success':opteffts_geo.success,'Effts':opteffts_geo.fun,'gparams':opteffts_geo.x}
    return outputdict

def ComputeR4vart(x:list) -> 'Effts':

    try:
        
        global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
        global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
        global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
        global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
        global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
        global Effts,Efftt,Efftspred,Reaction,vNu

        tenflow_coeff   = x[0]
        tenwork_coeff   = x[1]


        flow_coeff=tenflow_coeff/10
        work_coeff=tenwork_coeff/10

        # cycledict=whichcycle (cyclenum)       # The cycle to be computed
        # globals().update(cycledict)
        # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
        # globals().update(gparamdict)
        # rpm =whatrpm(m)          # rpm at m will be used
        # P_1 = P_1*10**6   sudah diubah jadi pa di fungsi whichcycle
        # P_5 = P_5*10**6

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

        h4s    = h04s-1/2*C4**2
        rho4s   = Props('D','H',h4s,'S',s04s,fluid)
        rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)

        initr4guess=0.045 # => initial guess r4 = 4.5cm
        r40     = initr4guess
        r4      = 0.04  # => set asal untuk memulai loop
        while np.abs(r40-r4)/r40 > 0.01: # residual/error harus lebih kecil dari 1%
            r40     = r4
            if 0.04*r4>0.001:
                tb4 = 0.04*r4
            else:
                tb4 = 0.001
            Bk4     = (tb4*0.05*NR)/(2*np.pi*r40*np.cos(Beta4))
            r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s*(1-Bk4))) # mflow sebagai input

        angvel  = U4/r4
        rpm     = angvel*(60/(2*np.pi))
        r5      = Rr5r4*r4
        b4      = Rb4r4*r4
        b5      = Rb5b4*b4
        rs5     = (2*r5+b5)/2
        rh5     = rs5-b5

        if 0.02*r4>0.001:
            tb5 = 0.02*r4
        else:
            tb5 = 0.001


        Ct5 = 0 # => it is predetermined that Alpha5=0
        Alpha5 = 0
        Cm5_0    = 10
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
            
            C5      = np.sqrt(Cm5i**2+Ct5**2)
            U5      = U4*Rr5r4
            W5      = np.sqrt(Cm5i**2+(U5-Ct5)**2)
            Beta5   = np.arccos(Cm5i/W5)

            Bk5     = (tb5*0.05*NR)/(2*np.pi*r5*np.cos(Beta5))
            Cm5ii       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ssi)*(1-Bk4)/(1-Bk5) *Cm4
            # Cm4ii       = mflow/(2*np.pi()*b5*)
            h5ss         = h05ss-1/2*(Cm5ii**2+Ct5**2)
            rho5ssii     = Props('D','H',h5ss,'S',s05ss,fluid)
            errorCm5    = np.abs((Rb5b4*Rr5r4*(rho5ssii/rho4s)*(Cm5ii/Cm4)*((1-Bk5)/(1-Bk4)))-1)
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
            
            C5      = np.sqrt(Cm5**2+Ct5**2)
            U5      = U4*Rr5r4
            W5       = np.sqrt(Cm5**2+(U5-Ct5)**2)
            Beta5   = np.arccos(Cm5/W5)

            Bk5     = (tb5*0.05*NR)/(2*np.pi*r5*np.cos(Beta5))

            Cm5       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ss)*(1-Bk4)/(1-Bk5) *Cm4 # -_- -_- -_-
            h5ss       = h05ss-1/2*(Cm5**2+Ct5**2)
            rho5ss     = Props('D','H',h5ss,'S',s05ss,fluid)
            if np.abs(1-Cm5/Props('A','H',h5ss,'S',s05ss,fluid)) < 5*1e-3:
                choked5 = True
                break
            errorCm5  = np.abs((Rb5b4*Rr5r4*(rho5ss/rho4s)*(Cm5/Cm4)*((1-Bk5)/(1-Bk4)))-1)
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
        # r4      = U4/np.radians(rpm*6)
        # r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s)) # mflow sebagai input
        # angvel  = U4/r4
        # rpm     = angvel*(60/(2*np.pi))
        # r5      = Rr5r4*r4
        # b4      = Rb4r4*r4
        # b5      = Rb5b4*b4
        # rs5     = (2*r5+b5)/2
        # rh5     = rs5-b5
        # if rh5 < 0.0015:
        #     print(f"For flow coeff ={flow_coeff} and work coeff={work_coeff} rh5 too small(<1.5mm), adjust gparams")
        
        Zr      = RZrr4*r4

        # mflow   = 2*np.pi*b5*r5*rho5ss*Cm5

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
        tesdict     = dict()

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
        for i in ('Beta4','Beta5','b4','r4','Zr','rs5','rh5'):
            proceeddict[i]  = globals()[i]
        tesdict['rpm'] = rpm

        outputdict  = {
            'geometry'  : geomdict,
            'thermo'    : thermodict,
            'velocity'  : veltridict,
            'efficiency': effdict,
            'losses'    : lossdict,
            'proceed'   : proceeddict,
            'tes'       : tesdict
            }

        return Effts*-1
    except ValueError:
        return 0


def rh5constrvart(x:list):

    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu

    tenflow_coeff   = x[0]
    tenwork_coeff   = x[1]

    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    # cycledict=whichcycle (cyclenum)       # The cycle to be computed
    # globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
    # P_1 = P_1*10**6 sudah diubah jadi pa di fungsi whichcycle
    # P_5 = P_5*10**6
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

    h4s    = h04s-1/2*C4**2
    rho4s   = Props('D','H',h4s,'S',s04s,fluid)
    rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)
    
    initr4guess=0.045 # => initial guess r4 = 4.5cm
    r40     = initr4guess
    r4      = 0.04  # => set asal untuk memulai loop
    while np.abs(r40-r4)/r40 > 0.01: # residual/error harus lebih kecil dari 1%
        r40     = r4
        if 0.04*r4>0.001:
            tb4 = 0.04*r4
        else:
            tb4 = 0.001
        Bk4     = (tb4*0.05*NR)/(2*np.pi*r40*np.cos(Beta4))
        r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s*(1-Bk4))) # mflow sebagai input
    angvel  = U4/r4
    rpm     = angvel*(60/(2*np.pi))
    r5      = Rr5r4*r4
    b4      = Rb4r4*r4
    b5      = Rb5b4*b4
    rs5     = (2*r5+b5)/2
    rh5     = rs5-b5
    if 0.02*r4>0.001:
        tb5 = 0.02*r4
    else:
        tb5 = 0.001
    
    return rh5-0.0075

def rs5constrvart(x:list):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu

    tenflow_coeff   = x[0]
    tenwork_coeff   = x[1]

    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    # cycledict=whichcycle (cyclenum)       # The cycle to be computed
    # globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
    # P_1 = P_1*10**6 sudah diubah jadi pa di fungsi whichcycle
    # P_5 = P_5*10**6
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

    h4s    = h04s-1/2*C4**2
    rho4s   = Props('D','H',h4s,'S',s04s,fluid)
    rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)

    initr4guess=0.045 # => initial guess r4 = 4.5cm
    r40     = initr4guess
    r4      = 0.04  # => set asal untuk memulai loop
    while np.abs(r40-r4)/r40 > 0.01: # residual/error harus lebih kecil dari 1%
        r40     = r4
        if 0.04*r4>0.001:
            tb4 = 0.04*r4
        else:
            tb4 = 0.001
        Bk4     = (tb4*0.05*NR)/(2*np.pi*r40*np.cos(Beta4))
        r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s*(1-Bk4))) # mflow sebagai input
    angvel  = U4/r4
    rpm     = angvel*(60/(2*np.pi))
    r5      = Rr5r4*r4
    b4      = Rb4r4*r4
    b5      = Rb5b4*b4
    rs5     = (2*r5+b5)/2
    rh5     = rs5-b5
    if 0.02*r4>0.001:
        tb5 = 0.02*r4
    else:
        tb5 = 0.001

    return (r4-rs5)/r4-0.2

def sonicconstrvart(x:list):
    
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    tenflow_coeff   = x[0]
    tenwork_coeff   = x[1]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    # cycledict=whichcycle (cyclenum)       # The cycle to be computed
    # globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
    # P_1 = P_1*10**6 sudah diubah jadi pa di fungsi whichcycle
    # P_5 = P_5*10**6
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
    h4s    = h04s-1/2*C4**2
    rho4s   = Props('D','H',h4s,'S',s04s,fluid)
    a4      = Props('A','H',h4s,'S',s04s,fluid)

    return (a4-C4)/a4-0.025

def Beta4constrvart(x:list):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu
    tenflow_coeff   = x[0]
    tenwork_coeff   = x[1]
    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    # cycledict=whichcycle (cyclenum)       # The cycle to be computed
    # globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
    # P_1 = P_1*10**6 sudah diubah jadi pa di fungsi whichcycle
    # P_5 = P_5*10**6
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
    return 70-np.abs(np.degrees(Beta4)) #untuk batasi Beta4, ikuti format: MaxAngle[deg] - abs(deg(Beta4))

def GetCoeffs(gparams,mflows) : # gparams 
    global  Rr5r4,Rb5b4,Rb4r4,RZrr4,NR,\
            mflow
    mflow   = mflows
    Rr5r4   = gparams[0]
    Rb5b4   = gparams[1]
    Rb4r4   = gparams[2]
    RZrr4   = gparams[3]
    NR      = gparams[4]

    tenflow_coeffb      = (0.6,7)
    tenwork_coeffb      = (7,30)
    bnds                = (tenflow_coeffb,tenwork_coeffb)

    constrrh5   = {'type': 'ineq', 'fun':rh5constrvart}
    constrrs5   = {'type': 'ineq', 'fun':rs5constrvart}
    constrsonic = {'type': 'ineq', 'fun':sonicconstrvart}
    constrbeta4 = {'type': 'ineq', 'fun':Beta4constrvart}
    constr      = [constrrh5,constrrs5,constrsonic,constrbeta4]
    initval     = [1,15]

    opteffts_coeffs = optimize.minimize(ComputeR4vart,initval,method='SLSQP',bounds=bnds,constraints=constr)

    outputdict  = {'success':opteffts_coeffs.success,'Effts':opteffts_coeffs.fun,'coeffs':opteffts_coeffs.x}
    return outputdict


def ComputeR4varm(x) -> 'Effts': 

    try:
        
        global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
        global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
        global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
        global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
        global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
        global Effts,Efftt,Efftspred,Reaction,vNu

        mflow   = x
        


        flow_coeff=tenflow_coeff/10
        work_coeff=tenwork_coeff/10

        # cycledict=whichcycle (cyclenum)       # The cycle to be computed
        # globals().update(cycledict)
        # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
        # globals().update(gparamdict)
        # rpm =whatrpm(m)          # rpm at m will be used
        # P_1 = P_1*10**6   sudah diubah jadi pa di fungsi whichcycle
        # P_5 = P_5*10**6

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

        h4s    = h04s-1/2*C4**2
        rho4s   = Props('D','H',h4s,'S',s04s,fluid)
        rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)

        initr4guess=0.045 # => initial guess r4 = 4.5cm
        r40     = initr4guess
        r4      = 0.04  # => set asal untuk memulai loop
        while np.abs(r40-r4)/r40 > 0.01: # residual/error harus lebih kecil dari 1%
            r40     = r4
            if 0.04*r4>0.001:
                tb4 = 0.04*r4
            else:
                tb4 = 0.001
            Bk4     = (tb4*0.05*NR)/(2*np.pi*r40*np.cos(Beta4))
            r4      = np.sqrt(mflow/ (2*np.pi*Rb4r4*Cm4*rho4s*(1-Bk4)) ) # mflow sebagai input

        angvel  = U4/r4
        rpm     = angvel*(60/(2*np.pi))
        r5      = Rr5r4*r4
        b4      = Rb4r4*r4
        b5      = Rb5b4*b4
        rs5     = (2*r5+b5)/2
        rh5     = rs5-b5

        if 0.02*r4>0.001:
            tb5 = 0.02*r4
        else:
            tb5 = 0.001

        Ct5 = 0 # => it is predetermined that Alpha5=0
        Alpha5 = 0
        Cm5_0    = 10
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

            C5      = np.sqrt(Cm5i**2+Ct5**2)
            U5      = U4*Rr5r4
            W5      = np.sqrt(Cm5i**2+(U5-Ct5)**2)
            Beta5   = np.arccos(Cm5i/W5)

            Bk5     = (tb5*0.05*NR)/(2*np.pi*r5*np.cos(Beta5))
            Cm5ii       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ssi)*(1-Bk4)/(1-Bk5) *Cm4
            # Cm4ii       = mflow/(2*np.pi()*b5*)
            h5ss         = h05ss-1/2*(Cm5ii**2+Ct5**2)
            rho5ssii     = Props('D','H',h5ss,'S',s05ss,fluid)
            errorCm5    = np.abs((Rb5b4*Rr5r4*(rho5ssii/rho4s)*(Cm5ii/Cm4)*((1-Bk5)/(1-Bk4)))-1)
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
            
            C5      = np.sqrt(Cm5**2+Ct5**2)
            U5      = U4*Rr5r4
            W5       = np.sqrt(Cm5**2+(U5-Ct5)**2)
            Beta5   = np.arccos(Cm5/W5)

            Bk5     = (tb5*0.05*NR)/(2*np.pi*r5*np.cos(Beta5))

            Cm5       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ss)*(1-Bk4)/(1-Bk5) *Cm4 # -_- -_- -_-
            h5ss       = h05ss-1/2*(Cm5**2+Ct5**2)
            rho5ss     = Props('D','H',h5ss,'S',s05ss,fluid)
            if np.abs(1-Cm5/Props('A','H',h5ss,'S',s05ss,fluid)) < 5*1e-3:
                choked5 = True
                break
            errorCm5  = np.abs((Rb5b4*Rr5r4*(rho5ss/rho4s)*(Cm5/Cm4)*((1-Bk5)/(1-Bk4)))-1)
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
        # r4      = U4/np.radians(rpm*6)
        # r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s)) # mflow sebagai input
        # angvel  = U4/r4
        # rpm     = angvel*(60/(2*np.pi))
        # r5      = Rr5r4*r4
        # b4      = Rb4r4*r4
        # b5      = Rb5b4*b4
        # rs5     = (2*r5+b5)/2
        # rh5     = rs5-b5
        # if rh5 < 0.0015:
        #     print(f"For flow coeff ={flow_coeff} and work coeff={work_coeff} rh5 too small(<1.5mm), adjust gparams")
        
        Zr      = RZrr4*r4

        # mflow   = 2*np.pi*b5*r5*rho5ss*Cm5

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
        tesdict     = dict()

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
        for i in ('Beta4','Beta5','b4','r4','Zr','rs5','rh5'):
            proceeddict[i]  = globals()[i]
        tesdict['rpm'] = rpm

        outputdict  = {
            'geometry'  : geomdict,
            'thermo'    : thermodict,
            'velocity'  : veltridict,
            'efficiency': effdict,
            'losses'    : lossdict,
            'proceed'   : proceeddict,
            'tes'       : tesdict
            }

        return Effts*-1
    except ValueError:
        return 0


def rh5constrvarm(x):

    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu

    mflow   = x

    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    # cycledict=whichcycle (cyclenum)       # The cycle to be computed
    # globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
    # P_1 = P_1*10**6 sudah diubah jadi pa di fungsi whichcycle
    # P_5 = P_5*10**6
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


    h4s    = h04s-1/2*C4**2
    rho4s   = Props('D','H',h4s,'S',s04s,fluid)
    rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)

    initr4guess=0.045 # => initial guess r4 = 4.5cm
    r40     = initr4guess
    r4      = 0.04  # => set asal untuk memulai loop
    while np.abs(r40-r4)/r40 > 0.01: # residual/error harus lebih kecil dari 1%
        r40     = r4
        if 0.04*r4>0.001:
            tb4 = 0.04*r4
        else:
            tb4 = 0.001
        Bk4     = (tb4*0.05*NR)/(2*np.pi*r40*np.cos(Beta4))
        r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s*(1-Bk4))) # mflow sebagai input
    angvel  = U4/r4
    rpm     = angvel*(60/(2*np.pi))
    r5      = Rr5r4*r4
    b4      = Rb4r4*r4
    b5      = Rb5b4*b4
    rs5     = (2*r5+b5)/2
    rh5     = rs5-b5
    if 0.02*r4>0.001:
        tb5 = 0.02*r4
    else:
        tb5 = 0.001


    return rh5-0.0075

def rs5constrvarm(x):
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu

    mflow   = x

    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    # cycledict=whichcycle (cyclenum)       # The cycle to be computed
    # globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
    # P_1 = P_1*10**6 sudah diubah jadi pa di fungsi whichcycle
    # P_5 = P_5*10**6
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
    h4s    = h04s-1/2*C4**2
    rho4s   = Props('D','H',h4s,'S',s04s,fluid)
    rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)
    
    initr4guess=0.045 # => initial guess r4 = 4.5cm
    r40     = initr4guess
    r4      = 0.04  # => set asal untuk memulai loop
    while np.abs(r40-r4)/r40 > 0.01: # residual/error harus lebih kecil dari 1%
        r40     = r4
        if 0.04*r4>0.001:
            tb4 = 0.04*r4
        else:
            tb4 = 0.001
        Bk4     = (tb4*0.05*NR)/(2*np.pi*r40*np.cos(Beta4))
        r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s*(1-Bk4))) # mflow sebagai input
    angvel  = U4/r4
    rpm     = angvel*(60/(2*np.pi))
    r5      = Rr5r4*r4
    b4      = Rb4r4*r4
    b5      = Rb5b4*b4
    rs5     = (2*r5+b5)/2
    rh5     = rs5-b5
    if 0.02*r4>0.001:
        tb5 = 0.02*r4
    else:
        tb5 = 0.001
    
    return (r4-rs5)/r4-0.2

def sonicconstrvarm(x):
    
    global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
    global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
    global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
    global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5
    global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
    global Effts,Efftt,Efftspred,Reaction,vNu

    mflow   = x

    flow_coeff=tenflow_coeff/10
    work_coeff=tenwork_coeff/10
    # cycledict=whichcycle (cyclenum)       # The cycle to be computed
    # globals().update(cycledict)
    # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
    # globals().update(gparamdict)
    # rpm =whatrpm(m)          # rpm at m will be used
    # P_1 = P_1*10**6 sudah diubah jadi pa di fungsi whichcycle
    # P_5 = P_5*10**6
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
    h4s    = h04s-1/2*C4**2
    rho4s   = Props('D','H',h4s,'S',s04s,fluid)
    a4      = Props('A','H',h4s,'S',s04s,fluid)

    return (a4-C4)/a4-0.025




def GetMFlow(gparams,coeffs):
    global  Rr5r4,Rb5b4,Rb4r4,RZrr4,NR,\
            tenflow_coeff,tenwork_coeff

    Rr5r4   = gparams[0]
    Rb5b4   = gparams[1]
    Rb4r4   = gparams[2]
    RZrr4   = gparams[3]
    NR      = gparams[4]

    tenflow_coeff = coeffs [0]
    tenwork_coeff = coeffs [1]
    


    mflowb  = (0.5,2)
    bnds    = (mflowb)
    
    constrrh5   = {'type': 'ineq', 'fun': rh5constrvarm}
    constrrs5   = {'type': 'ineq', 'fun': rs5constrvarm}
    constrsonic = {'type': 'ineq', 'fun': sonicconstrvarm}
    constr      = [constrrh5,constrrs5,constrsonic]
    initval     = [1]

    opteffts_mflow = optimize.minimize_scalar(ComputeR4varm,initval,method='bounded',bounds=mflowb)

    outputdict  = {'success':opteffts_mflow.success,'Effts':opteffts_mflow.fun,'mflow':opteffts_mflow.x}
    return outputdict


def ComputeR4all(cyclenum,mflows,coeffs:list,gparams:list) -> 'Effts':

    try:
        
        global T_1,T_5,P_1,P_5,p04s,p04,p4s,p4,p05ss,p5ss,p05,p5,T05ss,T05,T5ss,T5,rho4s,rho5ss,mflow,h4s,h04s
        global r4,r5,rs5,rh5,b4,b5,Zr,NR,tb4,tb5
        global C4,Ct4,Cm4,W4,U4,Alpha4,Beta4,C5,Ct5,Cm5,W5,U5,Alpha5,Beta5,Beta4opt,Beta4opt2
        global Cm5didconverge1,Cm5didconverge2,k1Cm5,k2Cm5,errorC5,Ma4s,rpm
        global TotalLoss,LossInc,LossInc0,LossPass,LossTip,LossWind,LossTE,LossExit,S5,O5,Effreductbladeloading
        global Effts,Efftt,Efftspred,Reaction,vNu
        global Beta4,Beta5,b4,r4,Zr,rs5,rh5

        cycledict=whichcycle (cyclenum)       # The cycle to be computed
        globals().update(cycledict)

        mflow   = mflows

        tenflow_coeff   = coeffs[0]
        tenwork_coeff   = coeffs[1]

        Rr5r4   = gparams[0]
        Rb5b4   = gparams[1]
        Rb4r4   = gparams[2]
        RZrr4   = gparams[3]
        NR      = gparams[4]

        flow_coeff=tenflow_coeff/10
        work_coeff=tenwork_coeff/10

        # cycledict=whichcycle (cyclenum)       # The cycle to be computed
        # globals().update(cycledict)
        
        # gparamdict=whichgparamset (l)  # geometry parameter set to be used is the l
        # globals().update(gparamdict)
        # rpm =whatrpm(m)          # rpm at m will be used
        # P_1 = P_1*10**6   sudah diubah jadi pa di fungsi whichcycle
        # P_5 = P_5*10**6

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

        h4s    = h04s-1/2*C4**2
        rho4s   = Props('D','H',h4s,'S',s04s,fluid)
        rho05ss = Props('D','H',h05ss,'S',s05ss,fluid)
        a4      = Props('A','H',h4s,'S',s04s,fluid)
        Ma4s    = C4/a4

        initr4guess=0.045 # => initial guess r4 = 4.5cm
        r40     = initr4guess
        r4      = 0.04  # => set asal untuk memulai loop
        while np.abs(r40-r4)/r40 > 0.01: # residual/error harus lebih kecil dari 1%
            r40     = r4
            if 0.04*r4>0.001:
                tb4 = 0.04*r4
            else:
                tb4 = 0.001
            Bk4     = (tb4*0.05*NR)/(2*np.pi*r40*np.cos(Beta4))
            r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s*(1-Bk4))) # mflow sebagai input

        angvel  = U4/r4
        rpm     = angvel*(60/(2*np.pi))
        r5      = Rr5r4*r4
        b4      = Rb4r4*r4
        b5      = Rb5b4*b4
        rs5     = (2*r5+b5)/2
        rh5     = rs5-b5

        if 0.02*r4>0.001:
            tb5 = 0.02*r4
        else:
            tb5 = 0.001

        Ct5 = 0 # => it is predetermined that Alpha5=0
        Alpha5 = 0
        Cm5_0    = 10
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
            rho5ssi     = rho5ssii

            C5      = np.sqrt(Cm5i**2+Ct5**2)
            U5      = U4*Rr5r4
            W5      = np.sqrt(Cm5i**2+(U5-Ct5)**2)
            Beta5   = np.arccos(Cm5i/W5)

            Bk5     = (tb5*0.05*NR)/(2*np.pi*r5*np.cos(Beta5))
            Cm5ii       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ssi)*(1-Bk4)/(1-Bk5) *Cm4
            # Cm4ii       = mflow/(2*np.pi()*b5*)
            h5ss         = h05ss-1/2*(Cm5ii**2+Ct5**2)
            rho5ssii     = Props('D','H',h5ss,'S',s05ss,fluid)
            errorCm5    = np.abs((Rb5b4*Rr5r4*(rho5ssii/rho4s)*(Cm5ii/Cm4)*((1-Bk5)/(1-Bk4)))-1)
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
            
            C5      = np.sqrt(Cm5**2+Ct5**2)
            U5      = U4*Rr5r4
            W5       = np.sqrt(Cm5**2+(U5-Ct5)**2)
            Beta5   = np.arccos(Cm5/W5)

            Bk5     = (tb5*0.05*NR)/(2*np.pi*r5*np.cos(Beta5))

            Cm5       = (1/(Rb5b4*Rr5r4))*(rho4s/rho5ss)*(1-Bk4)/(1-Bk5) *Cm4 # -_- -_- -_-
            h5ss       = h05ss-1/2*(Cm5**2+Ct5**2)
            rho5ss     = Props('D','H',h5ss,'S',s05ss,fluid)
            if np.abs(1-Cm5/Props('A','H',h5ss,'S',s05ss,fluid)) < 5*1e-3:
                choked5 = True
                break
            errorCm5  = np.abs((Rb5b4*Rr5r4*(rho5ss/rho4s)*(Cm5/Cm4)*((1-Bk5)/(1-Bk4)))-1)
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
        # r4      = U4/np.radians(rpm*6)
        # r4      = np.sqrt(mflow/(2*np.pi*Rb4r4*Cm4*rho4s)) # mflow sebagai input
        # angvel  = U4/r4
        # rpm     = angvel*(60/(2*np.pi))
        # r5      = Rr5r4*r4
        # b4      = Rb4r4*r4
        # b5      = Rb5b4*b4
        # rs5     = (2*r5+b5)/2
        # rh5     = rs5-b5
        # if rh5 < 0.0015:
        #     print(f"For flow coeff ={flow_coeff} and work coeff={work_coeff} rh5 too small(<1.5mm), adjust gparams")
        
        Zr      = RZrr4*r4

        # mflow   = 2*np.pi*b5*r5*rho5ss*Cm5

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
        # a01     = Props('A','P',p01,'T',T01,fluid)
        # a4s     = Props('A','P',p4s,'T',T4s,fluid)
        # a4      = Props('A','P',p4,'T',T4,fluid)
        # Ma4s    = C4/a4s
        # Ma4     = C4/a4
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
        inputdict     = dict()

        for i in ('r4','r5','rs5','rh5','b4','b5','Zr','NR','tb4','tb5'):
            geomdict[i]     = globals()[i]
        for i in ('T_1','T_5','P_1','P_5','p04s','p04','p4s','p4','p5ss','p5','p05ss','p05','T05ss','T05','T5ss','T5','rho4s','rho5ss','mflow','h4s','h04s'):
            thermodict[i]   = globals()[i]
        for i in ('C4','Ct4','Cm4','W4','U4','Alpha4','Beta4','C5','Ct5','Cm5','W5','U5','Alpha5','Beta5','Beta4opt','Beta4opt2','Cm5didconverge1','Cm5didconverge2','k1Cm5','k2Cm5','Ma4s'):
            veltridict[i]   = globals()[i]
        for i in ('Reaction','Effts','Efftt','rpm'):
            effdict[i]      = globals()[i]
        for i in ('LossInc0','LossInc','LossPass','LossTip','LossWind','LossTE','Effreductbladeloading','LossExit'):
            lossdict[i]     = globals()[i]
        for i in ('Beta4','Beta5','b4','r4','Zr','rs5','rh5','r5','b5','tb4','tb5',\
            'Cm4','U4','Ct4'):
            proceeddict[i]  = globals()[i]
        tonozzlelist        = NR,r4,Alpha4,b4,Ct4,rho4,mflow
        # for i in ('mflows','coeffs','gparams'):
            # inputdict       = locals()[i]
        inputdict['mflows'] = mflows
        inputdict['coeffs'] = coeffs
        inputdict['gparams']= gparams
        outputdict  = {
            'geometry'  : geomdict,
            'thermo'    : thermodict,
            'velocity'  : veltridict,
            'efficiency': effdict,
            'losses'    : lossdict,
            'proceed'   : proceeddict,
            'tonozzle'    : tonozzlelist,
            'input'     : inputdict
            }

        return outputdict
    except ValueError:
        return 0


def optimizeR(cyclenum):
    global  tenflow_coeff,tenwork_coeff,\
            Rr5r4,Rb5b4,Rb4r4,RZrr4,NR, \
            mflow,T_1,P_1,T_5,P_5    
    ''' initial '''
    # cyclenum        = 10
    cycledict       = whichcycle (cyclenum)    
    globals().update(cycledict)
    tenflow_coeff   = 1.2       # initial 1.2 range(0.8,)
    tenwork_coeff   = 15        # 15
    coeffs      = [tenflow_coeff,tenwork_coeff]

    ''' ------- '''

    gparamssol  = GetGParams(coeffs=coeffs,mflows=mflow)
    coeffssol   = GetCoeffs(gparams=gparamssol['gparams'],mflows=mflow)
    print(f"Iteration 0 : Effts = {gparamssol['Effts']*-1} %")
        # coeffsol    = GetCoeffs(gparams=gparamssol['gparams'],mflows=mflow) 

    for i in range(1,3+1):
        
        mflowsol    = GetMFlow  (gparams=gparamssol['gparams'], coeffs=coeffssol['coeffs'])
        coeffssol   = GetCoeffs (gparams=gparamssol['gparams'], mflows=mflowsol['mflow'])
        gparamssol  = GetGParams(coeffs=coeffssol['coeffs'],    mflows=mflowsol['mflow']) #tenflow_coeff,tenwork_coeff

        print(f"Iteration {i} : Effts = {gparamssol['Effts']*-1} % ")
        print(f"mflow       : {mflowsol['success']} {mflowsol['mflow']} ")
        print(f"coeffssol   : {coeffssol['success']} {coeffssol['coeffs']}")
        print(f"gparams     : {gparamssol['success']} {gparamssol['gparams']} ")

        # Rr5r4   = gparamssol['gparams'][0]
        # Rb5b4   = gparamssol['gparams'][1]
        # Rb4r4   = gparamssol['gparams'][2]
        # RZrr4   = gparamssol['gparams'][3]
        # NR      = gparamssol['gparams'][4]
        # coeffssol   = GetCoeffs() # gparamssol['gparams']

    optresult = ComputeR4all(cyclenum=cyclenum,mflows=mflowsol['mflow'],coeffs=coeffssol['coeffs'],gparams=gparamssol['gparams'])
    askprint = str(input("Do you want to print optimized result? [Y/N]: "))
    if askprint.upper() == 'Y' or askprint.upper() == 'YES':
        print(f"Optimized result for cycle {cyclenum} :")
        for i in optresult:
            print(i)
            print(optresult[i])
    
    for i in range(7):
        print(f"\x1b[2KCounting down...({7-i}s)",end="\r")
        sleep(1)

    return optresult['proceed']
    
    

    # print(coeffssol)
    # print(gparamssol['gparams'][1])
    # print((gparamssol['gparams']))
# optimizeR(10)