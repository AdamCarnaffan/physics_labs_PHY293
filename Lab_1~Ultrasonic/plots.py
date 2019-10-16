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
    return m*lambda_L/np.sin(theta)

#%% Get Data
data_1 = read_csv('lab_1~Ultrasonic/src1.txt')
plot_data = np.array(data_1)

# Get uncertainties
lambda_L = 0.000000614361
theta_L = 0.120855
d_sodium = 0.0254/2500

m_err = 0
theta_err = 0.000039269
d_err = 0.0000005
lambda_L_err = lambda_L * np.sqrt((d_err/d_sodium)**2 + (theta_err/theta_L)**2)

# Get averages
avg_theta = np.array([sum([plot_data[v][trial] for trial in range(2, len(plot_data[0]))])/(len(plot_data[0]) - 2) for v in range(0, len(plot_data))])
avg_lambda_s = [compute_lambda_s(plot_data[row][1], avg_theta[row]) for row in range(0, len(avg_theta))]
avg_lambda_err = [lambda_s*np.sqrt((lambda_L_err/lambda_L)**2 + (theta_err/theta)**2) for lambda_s, theta in zip(avg_lambda_s, avg_theta)]

#%% Data Fitting for Average
m, b = get_ls_line(1/plot_data[:,0], avg_lambda_s)
ls_fit = [b + m * x for x in 1/plot_data[:,0]]

#%% Plotting Average
plot.style.use('ggplot')
plot.errorbar(1/plot_data[:,0], avg_lambda_s, yerr=avg_lambda_err, 
                                            xerr=[0.005]*len(plot_data[:,0]), fmt='ro')
plot.title("Wavelength of Ultrasonic Waves at Different Frequencies", color='k')
plot.plot(1/plot_data[:,0], ls_fit)
plot.xlabel("1 / Frequency [1s x $10^{-6}$]")
plot.ylabel("Wavelength [m]")
fig = plot.gcf()
plot.figure()

#%% Save Figure
fig.savefig('sexyplot.png', facecolor='w')

#%% Fit Quality for Average
q = get_fit_quality_chi_sq(avg_lambda_s, ls_fit, avg_lambda_err)
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

#%% Plot Residuals
plot.style.use('ggplot')
plot.errorbar(1/plot_data[:,0], np.array(avg_lambda_s) - np.array(ls_fit), yerr=avg_lambda_err, 
                                            xerr=[0.005]*len(plot_data[:,0]), fmt='ro')
plot.title("Residuals of Ultrasonic Waves at Different Frequencies", color='k')
plot.plot(1/plot_data[:,0], [0]*11)
plot.xlabel("1 / Frequency [1s x $10^{-6}$]")
plot.ylabel("Standardized Residuals for Wavelength [m]")
fig = plot.gcf()
plot.figure()

#%%
