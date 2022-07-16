from Compute_Rotor import*
from Regress_Graph import*
from Compute_Nozzle import*
from Compute_Volute import*
import csv

fluid='R245fa'

#Compute Rotor (work,flow,fluid)
O=Compute(.2,1.1,fluid)
NR,r4,Alpha4,b4,Ct4,rho4,mflow,H_1,T_1,P_1,Zr,rs5,rh5=O[2],O[3],O[4],O[5],O[6],O[7],O[8],O[9],O[10],O[11],O[12],O[13],O[14]
# print('=================================================================================')
# print(NR,r4,Alpha4,b4,Ct4,rho4,mflow)

#Compute Nozzle
VarNozz=VariableNozz(NR,r4,Alpha4,b4,Ct4,rho4,mflow)
Betha2,Nn,r3,Ct3,Cm3,r2,s3,o3,Alpha3,t2Oc,t3Oc,tmaxOc,dOc=VarNozz[0],VarNozz[1],VarNozz[2],VarNozz[3],VarNozz[4],VarNozz[5],VarNozz[6],VarNozz[7],VarNozz[8],VarNozz[9],VarNozz[10],VarNozz[11],VarNozz[12]
ResultsNozz=np.array(VariantSearchNozz(),dtype=object)
SListB=np.array(ResultsNozz[1],dtype=object)
SListC=np.array(ResultsNozz[2],dtype=object)
SListA=np.array(ResultsNozz[0],dtype=object)
ErrorList=np.array(ResultsNozz[3],dtype=object)
Errmin=(np.min(ErrorList))
ErrminLoc=np.where(ErrorList==Errmin)
anew=SListA[ErrminLoc[0]]
bnew=SListB[ErrminLoc[0]]
cnew=SListC[ErrminLoc[0]]
anew=anew.astype(float)
bnew=bnew.astype(float)
cnew=cnew.astype(float)

#Evaluate Nozzle
NozzRess=CheckValNozz(anew,bnew,cnew)
# print('=================================================================================')
# print(NozzRess)

#Compute Volute
VoluteRess=ComputeVol(fluid,r2,r3,Cm3,Betha2,H_1,T_1,P_1,mflow)
# print('=================================================================================')
# print(VoluteRess)
# print('=================================================================================')

#Evaluasi Perhitungan dan Efisiensi

#Create 2D Quasi
n=4
Z5=0
Quas=QuasiNorm(b4,r4,Zr,rs5,rh5,n,Z5)

#Regress Zr
Z,r,Zh,rh,m,Ash,Bsh,Csh,Dsh,Esh,Fsh=Quas[0],Quas[1],Quas[2],Quas[3],Quas[4],Quas[5],Quas[6],Quas[7],Quas[8],Quas[9],Quas[10]
Zreg=Zrregs(Z,r,Zh,rh,m,Ash,Bsh,Csh,Dsh,Esh,Fsh)

#Define Meridional Plane
phis,tethas,Bethacs,phih,tethah,Bethach=Zreg[0],Zreg[1],Zreg[2],Zreg[3],Zreg[4],Zreg[5]
Mer=meridional(phis,tethas,Bethacs,phih,tethah,Bethach)

#Define Vector Component
X,Y,Z,Xh,Yh,Zh=Mer[0],Mer[1],Mer[2],Mer[3],Mer[4],Mer[5]
Vec=VectorComp(X,Y,Z,Xh,Yh,Zh)

#Define Thickness of blade and Coordinate
Txs,Tys,Tzs,Txh,Tyh,Tzh,tb4,tb5=Vec[0],Vec[1],Vec[2],Vec[3],Vec[4],Vec[5],Vec[6],Vec[7]
Cor3D=Coord3D(Txs,Tys,Tzs,Txh,Tyh,Tzh,tb4,tb5)

#Report Hasil Koordinat Rotor
# print()
# print('=================================================================================')
# print('Koordinat X Shroud')
# print(Cor3D[0])
# print('=================================================================================')
# print('Koordinat Y Shroud')
# print(Cor3D[1])
# print('=================================================================================')
# print('Koordinat Z Shroud')
# print(Cor3D[2])
# print('=================================================================================')
# print('Koordinat X Hub')
# print(Cor3D[3])
# print('=================================================================================')
# print('Koordinat Y Hub')
# print(Cor3D[4])
# print('=================================================================================')
# print('Koordinat Z Hub')
# print(Cor3D[5])
# print('=================================================================================')

#Report Hasil Koordinat Nozzle
Noz3D=nozzleconts(anew,bnew,cnew,t2Oc,t3Oc,dOc)
# print()
# print('=================================================================================')
# print('Koordinat X Atas')
# print(Noz3D[0])
# print('=================================================================================')
# print('Koordinat Y Atas')
# print(Noz3D[1])
# print('=================================================================================')
# print('Koordinat X Bawah')
# print(Noz3D[2])
# print('=================================================================================')
# print('Koordinat Y Bawah')
# print(Noz3D[3])

#Writing to TXT
with open('KoordRotor.txt', 'w') as f:
    f.write('=================================================================================\n')
    f.write('Koordinat X Shroud\n')
    f.write("%s\n" % Cor3D[0])
    f.write('=================================================================================\n')
    f.write('Koordinat Y Shroud\n')
    f.write("%s\n" % Cor3D[1])
    f.write('=================================================================================\n')
    f.write('Koordinat Z Shroud\n')
    f.write("%s\n" % Cor3D[2])
    f.write('=================================================================================\n')
    f.write('Koordinat X Hub\n')
    f.write("%s\n" % Cor3D[3])
    f.write('=================================================================================\n')
    f.write('Koordinat Y Hub\n')
    f.write("%s\n" % Cor3D[4])
    f.write('=================================================================================\n')
    f.write('Koordinat Z Hub\n')
    f.write("%s\n" % Cor3D[5])
    f.write('=================================================================================\n')

with open('KoordNozzle.txt', 'w') as g:
    g.write('=================================================================================\n')
    g.write('Koordinat X Atas\n')
    g.write("%s\n" % Noz3D[0])
    g.write('=================================================================================\n')
    g.write('Koordinat Y Atas\n')
    g.write("%s\n" % Noz3D[1])
    g.write('=================================================================================\n')
    g.write('Koordinat X Bawah\n')
    g.write("%s\n" % Noz3D[2])
    g.write('=================================================================================\n')
    g.write('Koordinat Y Bawah\n')
    g.write("%s\n" % Noz3D[3])

with open('TurbineSpecs.txt','w') as h:
    h.write('Main Rotor Specs\n')
    h.write('Number of Rotor Blade      : %2d pcs\n' % (NR))
    h.write('Alpha4                     : %2d Degree\n' % (Alpha4))
    h.write('Radius Inlet Rotor (r4)    : %2f m\n' % (r4))
    h.write('Radius Outlet Rotor (r5)   : %2f m\n' % ((rh5+rs5)/2))
    h.write('Inlet Width (b4)           : %2f m\n' % (b4))
    h.write('Outlet Width (b5)          : %2f m\n' % (rs5-rh5))
    h.write('Axial Length (Zr)          : %2f m\n' % (Zr))
    h.write('Tangential Speed Inlet     : %2f m/s\n' % (Ct4))
    h.write('Massflow                   : %2f kg/s\n' % (mflow))
    h.write('Inlet Density              : %2f kg/m3\n' % (rho4))
    h.write('=================================================================================\n')
    h.write('Main Nozzle Specs\n')
    h.write('Number of Nozzle Blade     : %2d pcs\n' % (Nn))
    h.write('Alpha3                     : %2f Degree\n' % (Alpha3))
    h.write('Radius Inlet Nozzle(r2)    : %2f m\n' % (r2))
    h.write('Radius Outlet Nozzle(r3)   : %2f m\n' % (r3))
    h.write('Inlet Width (b3)           : %2f m\n' % (b4))
    h.write('Tangential Speed Inlet     : %2f m/s\n' % (Ct3))
    h.write('a                          : %2f m\n' % (anew))
    h.write('b                          : %2f m\n' % (bnew))
    h.write('c                          : %2f m\n' % (cnew))
    h.write('=================================================================================\n')
    h.write('Main Volute Specs\n')
    h.write('Va                         : %2f m\n' % (VoluteRess[0]))
    h.write('Vb                         : %2f m\n' % (VoluteRess[1]))
    h.write('Radius Volute (r1)         : %2f m\n' % (VoluteRess[2]))
    h.write('=================================================================================\n')
    h.write('Volute Radius per Area\n')
    DegreeVolute=[0,45,90,135,180,225,270,315,360]
    for i in (DegreeVolute):
        h.write('%2d Degree\n' %(i))
        h.write('Va                : %2f m\n' % (i/360*VoluteRess[0]))
        h.write('Vb                : %2f m\n' % (i/360*VoluteRess[1]))
        h.write('=================================================================================\n')

f=open('RotorCoordinate.csv','w')
writer=csv.writer(f)
header=['X','Y','Z']
writer.writerow(header)
Koordinatx=list(map(lambda x:[x],Cor3D[0]))
Koordinaty=list(map(lambda y:[y],Cor3D[1]))
Koordinatz=list(map(lambda z:[z],Cor3D[2]))
data=(Koordinatx,Koordinaty,Koordinatz)
writer.writerow(data)
