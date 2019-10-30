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
    # y_diff = [ypt - fitpt for ypt, fitpt in zip(y, fit_y)]
    return sum([(ypt - fitpt)**2/acc**2 for ypt, fitpt, acc in zip(y, fit_y, y_acc)])

#%% Get Data
data_150 = np.array(read_csv('lab_2~Electron/src_v150.txt'))
data_200 = np.array(read_csv('lab_2~Electron/src_v200.txt'))
data_250 = np.array(read_csv('lab_2~Electron/src_v250.txt'))
data_300 = np.array(read_csv('lab_2~Electron/src_v300.txt'))
# plot_data = np.array(data_150)

# DEFINE CONSTANTS
B_e = 3.116 * 10**(-6)
R = 30.7/200
n = 130
mu_0 = 4 * np.pi * 10**(-7)
k = (1/np.sqrt(2)) * (4/5)**(3/2) * mu_0 * n / R
I_o = B_e / k
V_savg = np.sqrt((150.3 + 200.2 + 249.5 + 301.2)/4)

V_err = 0.1
I_err = 0.01
k_err = 0.0001

#%% Make Calculations
dts = [data_150, data_200, data_250, data_300]
plot_data = np.zeros((2, 8))
I_vals = []
avgs = []

for d in range(12, 4, -1):
    for st in dts:
        for row in st:
            if row[2] == d:
                if len(avgs)-1 < 12-d: avgs = avgs + [[]]
                avgs[12-d] += [np.sqrt(row[1])/((row[0]/100 - I_o)*1000)]
                if len(I_vals)-1 < 12-d: I_vals = I_vals + [[]]
                I_vals[12-d] += [row[0]]
                # print(row[0]/100 - I_o)

for d in range(12, 4, -1):
    plot_data[0,12-d] = d
    plot_data[1,12-d] = sum(avgs[12-d])
    I_vals[12-d] = sum(I_vals[12-d])/len(I_vals[12-d])

plot_data = np.transpose(plot_data)
plot_data[:,0] = (plot_data[:,0]/200)**2
plot_data[:,1] = 10*(plot_data[:,1]/k)**2

# print(plot_data)

# Calculate error
x_err = [x * np.sqrt((0.001/(np.sqrt(x)))**2 + (0.001/(np.sqrt(x)))**2) for x in plot_data[:,0]]
y_err = [y * np.sqrt((V_err/V_savg)**2 + (I_err/I_avg)**2 + (k_err/k)**2) for y, I_avg in zip(plot_data[:,1], I_vals)]

print(V_savg)
print(k)
print(I_vals)

#%% Data Fitting for Average
m, b = get_ls_line(plot_data[:,0], plot_data[:,1])
ls_fit = [b + m * x for x in plot_data[:,0]]

#%% Plotting Average
plot.style.use('ggplot')
plot.errorbar(plot_data[:,0], plot_data[:,1], yerr=y_err, 
                                            xerr=x_err, fmt='ro')
plot.title("The Radius of Orbit of an Electron Stimulated by a Magnetic Field", color='k')
plot.plot(plot_data[:,0], ls_fit)
plot.xlabel("Radius Squared [$m^2$]")
plot.ylabel("Power x Coil Characteristic (k) [$C m^2 kg^{-1}$]")
fig = plot.gcf()
fig.set_size_inches(11,8)
plot.figure(figsize=(18,16), dpi=80)

#%% Save Figure
fig.savefig('sexyplot.png', facecolor='w')

#%% Fit Quality for Average
q = get_fit_quality_chi_sq(plot_data[:,1], ls_fit, y_err)
N = 2 # Always 2 for linear fit, really DOF
print("reduced chi-squared: {}; chi-squared: {}; DOF: {};".format(q/N, q, N))
# Get uncertainty
pts_len = len(plot_data[:,0])
delta = pts_len*sum([(1/x)**2 for x in plot_data[:,0]]) - sum([1/x for x in plot_data[:,0]])**2
s_yxsq = (1/(pts_len - 2))*sum([(ypt - yest)**2 for ypt, yest in zip(plot_data[:,1], ls_fit)])
s_m = np.sqrt(pts_len*(s_yxsq/delta))
s_b = np.sqrt((s_yxsq*sum([(1/x)**2 for x in plot_data[:,0]]))/delta)
print("slope: {}; intercept: {};".format(m, b))
print("slope error: {}; intercept error: {};".format(s_m, s_b))

#%% Plot Residuals
plot.style.use('ggplot')
plot.errorbar(plot_data[:,0], np.array(plot_data[:,1]) - np.array(ls_fit), yerr=y_err, 
                                            xerr=x_err, fmt='ro')
plot.title("Risiduals of Radius of Orbit of an Electron Stimulated by a Magnetic Field", color='k')
plot.plot(plot_data[:,0], [0]*len(plot_data[:,0]))
plot.xlabel("Radius Squared [$m^2$]")
plot.ylabel("Risidual of Power x Coil Characteristic (k) [$C m^2 kg^{-1}$]")
fig = plot.gcf()
fig.set_size_inches(11,8)
plot.figure(figsize=(18,16), dpi=80)

#%%
fig.savefig('sexyresiduals.png', facecolor='w')

# %%
