from ifrturbinepackage.definitions import *
from ifrturbinepackage.inputs import *


from ifrturbinepackage.rotor import *


# from ifrturbinepackage.nozzle import *
# from ifrturbinepackage.volute import *

# print(T_1)
# cycledict=whichcycle(10)
# print(cycledict['fluid'])
# global T_1,T_5,P_1,P_5,fluid,mflow
# k=10
# dfcycle=pd.read_csv(os.path.join(ROOT_DIR,"Inputs\cyclelist.csv"),skiprows=1,header=0,index_col=0)
# T_1 = dfcycle.iloc[k-1]['T_1']
# P_1 = dfcycle.iloc[k-1]['P_1']
# T_5 = dfcycle.iloc[k-1]['T_5']
# T_5 = dfcycle.iloc[k-1]['P_5']
# fluid = dfcycle.iloc[k-1]['fluid']
# print(T_1)
n = 0

p = whichfitfun(n)

def mfitefftsforn(x):
    tenflow_coeff   = x[0]
    tenwork_coeff   = x[1]
    effts=0
    for i in range(0,6):
        for j in range (0,6):
            effts=effts+p[i][j]*tenflow_coeff**i*tenwork_coeff**j
    return effts*-1

tenflowb=(1,12)  # Constraint: 0<=10*flowcoeff<=9
tenworkb=(0,25)  # Constraint: 0<=10*workcoeff<=20

    
    
constr1 = {'type': 'ineq', 'fun':constraint1 }
constr2 = {'type': 'ineq', 'fun':constraint2 }
constr  = [constr1,constr2]
bnds=(tenflowb,tenworkb)
initval=[0.5,3]
opteffts_test=optimize.minimize(mfitefftsforn,initval,method='SLSQP',bounds=bnds,constraints=constr)
print(opteffts_test)
print(mfitefftsforn([2,12]))
# p=whichfitfun(0)
# print(p)