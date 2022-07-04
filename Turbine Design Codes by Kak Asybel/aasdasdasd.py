from Compute_Rotor import*
from Regress_Graph import*
from Compute_Nozzle import*
from Compute_Volute import*
NR,r4,Alpha4,b4,Ct4,rho4,mflow=22,0.0291,12.58,3.25064/1000,70.7,32.0081,0.3

#VariableNozzle
VariableNozz(NR,r4,Alpha4,b4,Ct4,rho4,mflow)
VarNozz=VariableNozz(NR,r4,Alpha4,b4,Ct4,rho4,mflow)
ResultsNozz=np.array(VariantSearchNozz())

SListB=np.array(ResultsNozz[0],dtype=object)
SListC=np.array(ResultsNozz[1],dtype=object)
SListA=np.array(ResultsNozz[2],dtype=object)
ErrorList=np.array(ResultsNozz[3],dtype=object)

# print('=================================================================================')
# print('A:')
# print(SListA)
# print('=================================================================================')
# print('B:')
# print(SListB)
# print('=================================================================================')
# print('C:')
# print(SListC)
# print('=================================================================================')
print('Error:')
print(ErrorList)
print('=================================================================================')


print(min(ErrorList))