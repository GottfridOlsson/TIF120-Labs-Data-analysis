import numpy as np
import CSV_handler as CSV
import os
import matplotlib.pyplot as plt


def scale_X(x, k, m):
    return k*x+m

# Read CSV #
CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
data = CSV.read(CURRENT_PATH + '\\data from lab\\GROUP1_QMS.csv')
header = CSV.get_header(data)
time    = data[header[0]]
voltage = data[header[1]]

# Take one period #
time_cut    = np.array(time[0:1064])
voltage_cut = np.array(voltage[0:1064])
voltage_cut = voltage_cut + np.min(-voltage_cut) # set min at zero intensity
indexes = np.array(range(0,1064))

# Peaks #
peak_ranges   = [(150,180), (470,500), (650,720), (900,930), (970,1030)]
peak_voltages = np.array([np.max(-voltage_cut[start:end]) for (start, end) in peak_ranges])
peak_times    = np.array([np.where(-voltage_cut == peak_voltage)[0][0] for peak_voltage in peak_voltages])

# Normalize # 
y_scaling = 1/np.max(peak_voltages)
voltage_cut *= y_scaling
peak_voltages *= y_scaling

# Rough approximation: Create x-axis from two assumed elements at certain peaks #
index1 = 167 
dalton1 = 2
index2 = 680 #999
dalton2 = 28 #44
# y=k*x+m with points (index, Dalton)
k = (dalton2-dalton1)/(index2-index1)
m = dalton1-k*index1

Daltons = scale_X(indexes, k, m)
Da_peak_times = scale_X(peak_times, k, m)


# Finer fit: create x-axis from many assumed elements with clearly defined peaks #
Daltons_expected = [1, 2, 15, 16, 17, 18, 28, 40, 44]
Daltons_expected_labels = ['H+', 'H2+', 'NH+', 'O2.2+ or CH4+', 'OH^+', 'H2O^+', 'N2^+', 'Ar^+', 'CO2^+']

Daltons_clear_peaks    = [2, 28, 40, 44]
indexes_at_clear_peaks = [167, 680, 919, 999]

k_fit, m_fit = np.polyfit(indexes_at_clear_peaks, Daltons_clear_peaks, deg=1)
print(f"Fit to 2 peaks: y_fit = {k}x + {m}")
print(f"Fit to {len(Daltons_clear_peaks)} peaks: y_fit = {k_fit}x + {m_fit}")
x_range = np.arange(0,1000)
y = [x*k+m for x in x_range]
y_fit = [x*k_fit+m_fit for x in x_range]
Daltons_good_fit = scale_X(indexes, k_fit, m_fit)

# Plots #
index_votlage_plot = False
fit_clear_peaks = True
good_fit_Da_plot = True

if index_votlage_plot:
    plt.plot(-voltage_cut)
    plt.xlabel('Index')
    plt.ylabel('Normalized intensity (arbitrary units)')
    plt.show()

if fit_clear_peaks:
    plt.scatter(indexes_at_clear_peaks, Daltons_clear_peaks, label='Clear peaks')
    plt.plot(x_range, y, label='2 peak line fit')
    plt.plot(x_range, y_fit, label='4 peak line fit')
    plt.xlabel('Index')
    plt.ylabel('Dalton (m/q)')
    plt.legend()
    plt.show()    



if good_fit_Da_plot:
    plt.plot(Daltons_good_fit, -voltage_cut)
    plt.vlines(Daltons_expected, 0, 1, 'k')
    plt.xlabel('Dalton (m/q)')
    plt.ylabel('Normalized intensity (arbitrary units)')
    plt.show()

CSV.print_arrays_to_CSV(CURRENT_PATH + '\\formatted data\\QMS_formatted.csv', 'Dalton (m/q)', Daltons_good_fit, 'Normalized intensity (arbitrary units)', -voltage_cut)
