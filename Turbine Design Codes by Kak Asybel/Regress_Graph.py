import numpy as np
import CoolProp
from CoolProp.CoolProp import PropsSI as Pr
from scipy.optimize import curve_fit
import random as random
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import random
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import cross_val_score
from sklearn import linear_model
from matplotlib import cm
from Compute_Rotor import Compute

def RegressG(MinVar1,MaxVar1,MinVar2,MaxVar2):
    global Efficiency_list1
    Efficiency_List1 = []
    Efficiency_List2 = []
    Efficiency_List3 = []
    Efficiency_List4 = []
    Efficiency_List5 = []
    Efficiency_List6 = []
    NSAMPLE=200
    Var1_List = []
    Var2_List = []
    R_List=[]
    for i in range(0,NSAMPLE):
        var1  = random.uniform(MinVar1,MaxVar1)
        var2  = random.uniform(MinVar2,MaxVar2)
        var3 = 'R245fa'
        var4 = 'R134a'
        var5 = 'R1234yf'
        Var1_List.append(var1)
        Var2_List.append(var2)
        Efficiency1 = Compute(var1,var2,var3)[0]
        Efficiency2 = Compute(var1,var2,var4)[0]
        Efficiency3 = Compute(var1,var2,var5)[0]
        Efficiency_List1.append(Efficiency1*100)
        Efficiency_List2.append(Efficiency2*100)
        reaction= Compute(var1,var2,var3)[1]
        R_List.append(reaction)

    #Regresi 2 Variabel
    datapoints = np.array(np.vstack((Var1_List,Var2_List,Efficiency_List1)).T)

    X = datapoints[:,0:2]
    Y = datapoints[:,-1]
    #  degree polynomial features
    deg_of_poly = 3
    poly = PolynomialFeatures(degree=deg_of_poly)
    X_ = poly.fit_transform(X)
    # Fit linear model
    clf = linear_model.LinearRegression()
    clf.fit(X_, Y)

    # The test set, or plotting set
    N = NSAMPLE
    LengthVar1 = MaxVar1
    LengthVar2 = MaxVar2
    predict_x0, predict_x1 = np.meshgrid(np.linspace(MinVar1, MaxVar1, N), 
                                        np.linspace(MinVar2, MaxVar2, N))
    predict_x = np.concatenate((predict_x0.reshape(-1, 1), 
                                predict_x1.reshape(-1, 1)), 
                            axis=1)
    predict_x_ = poly.fit_transform(predict_x)
    predict_y = clf.predict(predict_x_)

    # Plot
    fig = plt.figure(figsize=(40, 10))
    ax1 = fig.add_subplot(121, projection='3d')
    surf = ax1.plot_surface(predict_x0, predict_x1, predict_y.reshape(predict_x0.shape), 
                            rstride=1, cstride=1, cmap=cm.jet, alpha=0.5)
    ax1.scatter(datapoints[:, 0], datapoints[:, 1], datapoints[:, 2], c='b', marker='o')

    ax1.set_xlim((MinVar1, MaxVar1))
    ax1.set_ylim((MinVar2, MaxVar2))
    fig.colorbar(surf, ax=ax1)
    ax2 = fig.add_subplot(122)
    cs = ax2.contourf(predict_x0, predict_x1, predict_y.reshape(predict_x0.shape))
    ax2.contour(cs, colors='k')
    fig.colorbar(cs, ax=ax2)
    plt.show()
    KOEFISIEN=clf.coef_
    VARIABEL=poly.get_feature_names_out()
    print(KOEFISIEN)
    print(VARIABEL)