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

def compute_lambda_s(m, theta):
    sodium_lamp_wavelength = 0.000000614361
    return m*sodium_lamp_wavelength/np.sin(theta)

#%% Get Data
data_1 = read_csv('lab_1~Ultrasonic/src1.txt')
plot_data = np.array(data_1)
avg_theta = np.array([sum([plot_data[v][trial] for trial in range(2, len(plot_data[0]))])/(len(plot_data[0]) - 2) for v in range(0, len(plot_data))])
avg_lambda_s = [compute_lambda_s(plot_data[row][1], avg_theta[row]) for row in range(0, len(avg_theta))]

#%% Data Fitting for Average
m, b = get_ls_line(plot_data[:,0], avg_lambda_s)
ls_fit = [b + m * x for x in plot_data[:,0]]

#%% Plotting Average
plot.errorbar(1/plot_data[:,0], avg_lambda_s, yerr=[0.00003926990817]*len(avg_lambda_s), 
                                            xerr=[0.005]*len(plot_data[:,0]), fmt='ro')
plot.title("The Wavelength of Ultrasonic Waves in Water at Different Frequencies")
plot.plot(1/plot_data[:,0], ls_fit)
plot.xlabel("1 / Frequency [s]")
plot.ylabel("Wavelength [mm]")
plot.figure()

#%% Fit Quality for Average
q = get_fit_quality_chi_sq(avg_lambda_s, ls_fit, [0.00003926990817]*len(avg_lambda_s))
N = 2 # Always 2 for linear fit, really DOF
print("reduced chi-squared: {}; chi-squared: {}; DOF: {};".format(q/N, q, N))
# Get uncertainty
pts_len = len(plot_data[:,0])
delta = pts_len*sum([(1/x)**2 for x in plot_data[:,0]]) - sum([1/x for x in plot_data[:,0]])**2
s_yxsq = (1/(pts_len - 2))*sum([(ypt - yest)**2 for ypt, yest in zip(avg_lambda_s, ls_fit)])
s_m = np.sqrt(pts_len*(s_yxsq/delta))
s_b = np.sqrt((s_yxsq*sum([(1/x)**2 for x in plot_data[:,0]]))/delta)
print("slope: {}; intercept: {};".format(m, b))
print("slope error: {}; intercept error: {};".format(s_m, s_b))

#%%
