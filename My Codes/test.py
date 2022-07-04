# IMPORTS
import numpy as np
import CoolProp
from CoolProp.CoolProp import PropsSI as Props
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.Plots import PropertyPlot

plot = PropertyPlot('HEOS::R245fa', 'TS', unit_system='EUR', tp_limits='ORC')
plot.calc_isolines(CoolProp.iQ, num=11)
plot.calc_isolines(CoolProp.iP, iso_range=[1,50], num=10, rounding=True)
plot.draw()
plot.isolines.clear()
plot.props[CoolProp.iP]['color'] = 'green'
plot.props[CoolProp.iP]['lw'] = '0.5'
plot.calc_isolines(CoolProp.iP, iso_range=[1,50], num=10, rounding=False)
plot.show()