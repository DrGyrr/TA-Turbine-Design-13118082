from ifrturbinepackage.definitions import *
from ifrturbinepackage.inputs import *


from ifrturbinepackage.rotor import *


# from ifrturbinepackage.nozzle import *
# from ifrturbinepackage.volute import *

# print(T_1)
cycledict=whichcycle(10)
print(cycledict['fluid'])
# global T_1,T_5,P_1,P_5,fluid,mflow
# k=10
# dfcycle=pd.read_csv(os.path.join(ROOT_DIR,"Inputs\cyclelist.csv"),skiprows=1,header=0,index_col=0)
# T_1 = dfcycle.iloc[k-1]['T_1']
# P_1 = dfcycle.iloc[k-1]['P_1']
# T_5 = dfcycle.iloc[k-1]['T_5']
# T_5 = dfcycle.iloc[k-1]['P_5']
# fluid = dfcycle.iloc[k-1]['fluid']
# print(T_1)

