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
    return sum([(ypt - fitpt)**2/acc**2 for ypt, fitpt, acc in zip(y, fit_y, y_acc)])

#%% Get Data
data_1 = read_csv('lab_2~Electron/src_v150.txt')
plot_data = np.array(data_1)

#%% Data Fitting for Average
m, b = get_ls_line(plot_data[:,2], plot_data[:,0])
ls_fit = [b + m * x for x in plot_data[:,2]]

#%% Plotting Average
plot.style.use('ggplot')
plot.errorbar(plot_data[:,2], plot_data[:,0], yerr=[0.001]*len(plot_data[:,2]), 
                                            xerr=[0.1]*len(plot_data[:,0]), fmt='ro')
plot.title("Wavelength of Ultrasonic Waves at Different Frequencies", color='k')
plot.plot(plot_data[:,2], ls_fit)
plot.xlabel("1 / Frequency [1s x $10^{-6}$]")
plot.ylabel("Wavelength [m]")
fig = plot.gcf()
plot.figure()

#%% Save Figure
fig.savefig('sexyplot.png', facecolor='w')

#%% Fit Quality for Average
q = get_fit_quality_chi_sq(plot_data[:,0], ls_fit, [0.001]*len(plot_data[:,0]))
N = 2 # Always 2 for linear fit, really DOF
print("reduced chi-squared: {}; chi-squared: {}; DOF: {};".format(q/N, q, N))
# Get uncertainty
pts_len = len(plot_data[:,0])
# delta = pts_len*sum([(1/x)**2 for x in plot_data[:,2]]) - sum([1/x for x in plot_data[:,2]])**2
# s_yxsq = (1/(pts_len - 2))*sum([(ypt - yest)**2 for ypt, yest in zip(avg_lambda_s, ls_fit)])
# s_m = np.sqrt(pts_len*(s_yxsq/delta))
# s_b = np.sqrt((s_yxsq*sum([(1/x)**2 for x in plot_data[:,0]]))/delta)
print("slope: {}; intercept: {};".format(m, b))
# print("slope error: {}; intercept error: {};".format(s_m, s_b))

#%% Plot Residuals
plot.style.use('ggplot')
plot.errorbar(plot_data[:,2], np.array(plot_data[:,0]) - np.array(ls_fit), yerr=[0.001]*len(plot_data[:,0]), 
                                            xerr=[0.1]*len(plot_data[:,2]), fmt='ro')
plot.title("Residuals of Ultrasonic Waves at Different Frequencies", color='k')
plot.plot(plot_data[:,2], [0]*len(plot_data[:,2]))
plot.xlabel("1 / Frequency [1s x $10^{-6}$]")
plot.ylabel("Standardized Residuals for Wavelength [m]")
fig = plot.gcf()
plot.figure()

#%%
