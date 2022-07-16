from Compute_Rotor import*
from Regress_Graph import*
from Compute_Nozzle import*
from Compute_Volute import*

def QuasiNorm(b4,r4,Zr,rs5,rh5,n,Z5):
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
        Z.append(Z+dZ,dtype=object)
    Epsi = Z/(dZ-b4)
    r=rs5+(r4-rs5)*Epsi**n
    fz=np.sqrt(1+(C2*n*(Z-Z5)**(n-1))**2)
    m.append((dZ/3)*(fz+4)/fz) 
    for i in (fz):
        m.append(m+2*(i)/3)            #perbaiki value fz #check OK

    #Hub Sections
    Zh=[]
    Zh.append(0)
    mh=[]
    mh.append(0)
    for i in range(0,QuasiSec):
        Zh.append(Zh+dZ)
    Ephi = Zh/(dZ-b4)
    rh=rh5+Zr-np.sqrt((Zr**2)-(Zh**2))
    fzh=np.sqrt(1+(C2*n*(Z-Z5)**(n-1))**2)
    mh.append((dZ/3)*(fzh+4)/fzh) 
    for i in (fzh):
        mh.append(mh+2*(i)/3)           #perbaiki value mh #check OK


    Zm=(Z-Zh)/2
    rm=(r-rh)/2
    Betha5  = np.arccos(Cm5/W5)
    Betha5s = np.arctan(np.tan(Betha5)/rs5)
    Betha5h = np.arctan(np.tan(Betha5)/rh5)
    ms      = m[50]
    Ash     = 1/np.tan(Betha5s)/rs5
    Bsh     = (1/np.tan(Betha4)-1/np.tan(Betha5s))/(ms**2)
    Csh     = -Bsh/(2*ms)
    Tetha4  = ms/2*(1/np.tan(Betha4)/(r4+1/np.tan(Betha5))/rs5)
    mhh     = mh[50]
    Dsh     = 1/np.tan(Betha5h)/rh5
    Esh     = 3*Tetha4/(mhh**2)-(1/mhh)*(2*1/np.tan(Betha5h))
    Fsh     = 1/mhh**2*(1/np.tan(Betha5h)/rh5+1/np.tan(Betha4)/r4)-2*Tetha4/mhh**3
    return(Z,r,Zh,rh,m,Ash,Bsh,Csh,Dsh,Esh,Fsh)
    
def Zrregs(Z,r,Zh,rh,m,Ash,Bsh,Csh,Dsh,Esh,Fsh):
    #Shroud
    for i in range(0,len(Z)):
        dms=np.sqrt((Z-Z)**2+(r-r)**2)
    msi=[]
    msi.append(0)
    for i in range(0,len(Z)):
        msi=msi+dms
    polynom = 6             #Coeficient harus 6 (Belum otomatis)
    mymodel = np.poly1d(np.polyfit(Z,r,polynom))
    #need to regress the mymodel
    #Define the Function of phi
    coeffphis=[]
    for coeff in mymodel:
        coeffphis.append(coeff)
    phis=6*coeffphis[5]+5*coeffphis[4]+4*coeffphis[3]+3*coeffphis[2]+2*coeffphis[1]+coeffphis[0]
    tethas=msi*(Ash+Bsh**3+Csh**4)
    Bethacs=1/np.arctan(r*(Ash+3*Bsh**2+4*Csh**3))

    #Hub
    for i in range(0,len(Zh)):
        dmh=np.sqrt((Zh-Zh)**2+(rh-rh)**2)    
    mhi=[]
    mhi.append(0)
    for i in (dmh):
        mhi=mhi+i
    polynomh = 6
    mymodelh = np.poly1d(np.polyfit(Zh,rh,polynomh))
    coeffphih=[]
    for coeffh in mymodelh:
        coeffphih.append(coeffh)
    phih=6*coeffphih[5]+5*coeffphih[4]+4*coeffphih[3]+3*coeffphih[2]+2*coeffphih[1]+coeffphih[0]
    tethah=1/np.tan(rh*(mhi*(Dsh+Esh**2+Fsh**3)))
    Bethach=1/np.arctan(r*(Dsh+3*Esh**2+4*Fsh**3))
    return(phis,tethas,Bethacs,phih,tethah,Bethach)