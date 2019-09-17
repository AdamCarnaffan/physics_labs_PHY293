#%% Imports
from pandas import read_csv
import matplotlib.pyplot as plot
import numpy as np

def get_ls_line(x, y):
    """This returns the Slope & Y-intercept for the least-squares method of fitting"""
    pts_len = len(x)
    sum_x = sum([xpt for xpt in x])
    sum_y = sum([ypt for ypt in y])
    sum_xy = sum([xpt*ypt for xpt, ypt in zip(x, y)])
    sum_xsq = sum([xpt**2 for xpt in x])
    m = ((pts_len*sum_xy - sum_x*sum_y)/(pts_len*sum_xsq - sum_x**2))
    b = ((sum_y-m*sum_x)/pts_len)
    return m, b

def get_fit_quality(y, fit_y, y_acc):
    """This returns the chi-squared value for the data & fit line"""
    y_diff = [ypt - fitpt for ypt, fitpt in zip(y, fit_y)]
    print(np.array(y_diff)**2)
    return sum([(ypt - fitpt)**2/acc**2 for ypt, fitpt, acc in zip(y, fit_y, y_acc)])

#%% Get Data for Option 1
data_1 = read_csv('lab_0~Exercise/src1.txt')
plot_data = np.array(data_1)

#%% Data Fitting for Option 1
m, b = get_ls_line(plot_data[:,0], plot_data[:,2])
ls_fit = [b + m * x for x in plot_data[:,0]]

#%% Plotting Option 1
plot.errorbar(plot_data[:,0], plot_data[:,2], plot_data[:,1], plot_data[:,3], 'ro')
plot.title("Option 1 - V vs. I")
plot.plot(plot_data[:,0], ls_fit)
plot.xlabel("Voltage Measured [V]")
plot.ylabel("Current Measured [mA]")
plot.figure()

#%% Fit Quality for Option 1
q = get_fit_quality(plot_data[:,2], ls_fit, plot_data[:,3])
N = 2 # Always 2 for linear fit, really DOF
print("chi-squared: {}; DOF: {};".format(q, N))

#%% Get Data for Option 2
data_2 = read_csv('lab_0~Exercise/src2.txt')
plot_data = np.array(data_2)

#%% Data Fitting for Option 2
m, b = get_ls_line(plot_data[:,0], plot_data[:,2])
ls_fit = [b + m * x for x in plot_data[:,0]]

#%% Plotting Option 2
plot.errorbar(plot_data[:,0], plot_data[:,2], plot_data[:,1], plot_data[:,3], 'ro')
plot.title("Option 2 - V vs. I")
plot.xlabel("Voltage Measured [V]")
plot.ylabel("Current Measured [mA]")
plot.plot(plot_data[:,0], ls_fit)
plot.figure()


#%% Fit Quality for Option 2
q = get_fit_quality(plot_data[:,2], ls_fit, plot_data[:,3])
N = 2 # Always 2 for linear fit, really DOF
print("chi-squared: {}; DOF: {};".format(q, N))