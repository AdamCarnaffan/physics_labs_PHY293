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

def get_fit_quality_chi_sq(y, fit_y, y_acc):
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
q = get_fit_quality_chi_sq(plot_data[:,2], ls_fit, plot_data[:,3])
N = 2 # Always 2 for linear fit, really DOF
print("reduced chi-squared: {}; chi-squared: {}; DOF: {};".format(q/N, q, N))
# Get uncertainty
pts_len = len(plot_data[:,0])
delta = pts_len*sum([x**2 for x in plot_data[:,0]]) - sum([x for x in plot_data[:,0]])**2
s_yxsq = (1/(pts_len - 2))*sum([(ypt - yest)**2 for ypt, yest in zip(plot_data[:,3], ls_fit)])
s_m = np.sqrt(pts_len*(s_yxsq/delta))
s_b = np.sqrt((s_yxsq*sum([x**2 for x in plot_data[:,0]]))/delta)
print("slope error: {}; intercept error: {};".format(s_m, s_b))

#%% Get Data for Option 2
data_2 = read_csv('lab_0~Exercise/src2.txt')
plot_data = np.array(data_2)

#%% Run Rv Calcs for Option 2
r_vals = [1.2e6, 34e3, 845, 327]
r_err = [1e5, 1e2, 1, 1]
r_v = [v/(i - v/r) for v, i, r in zip(plot_data[:,0], plot_data[:,2], r_vals)]
r_v_err = []
for r_v_pt, v, verr, i, ierr, r, rerr in zip(r_v, plot_data[:,0], plot_data[:,1], plot_data[:,2], plot_data[:,3], r_vals, r_err):
    v_over_r_err = (v/r)*np.sqrt((verr/v)**2 + (rerr/r)**2)
    i_a_diff_err = np.sqrt(ierr**2 + v_over_r_err**2)
    r_v_err += [r_v_pt*np.sqrt((verr/v)**2 + ((i - v/r)/i_a_diff_err)**2)]
print(r_v)
print(r_v_err)

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
q = get_fit_quality_chi_sq(plot_data[:,2], ls_fit, plot_data[:,3])
N = 2 # Always 2 for linear fit, really DOF
print("reduced chi-squared: {}; chi-squared: {}; DOF: {};".format(q/N, q, N))
# Get uncertainty
pts_len = len(plot_data[:,0])
delta = pts_len*sum([x**2 for x in plot_data[:,0]]) - sum([x for x in plot_data[:,0]])**2
s_yxsq = (1/(pts_len - 2))*sum([(ypt - yest)**2 for ypt, yest in zip(plot_data[:,3], ls_fit)])
s_m = np.sqrt(pts_len*(s_yxsq/delta))
s_b = np.sqrt((s_yxsq*sum([x**2 for x in plot_data[:,0]]))/delta)
print("slope error: {}; intercept error: {};".format(s_m, s_b))