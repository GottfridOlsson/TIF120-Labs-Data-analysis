import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import CSV_handler as CSV
import os

matplotlib.rcParams.update({
    "text.usetex": True,
    "font.family": "serif", 
    "font.serif" : ["Computer Modern Roman"]
    })

axis=11
tick=9
legend=9
matplotlib.rc('font',   size=axis)      
matplotlib.rc('axes',   titlesize=axis) 
matplotlib.rc('axes',   labelsize=axis) 
matplotlib.rc('xtick',  labelsize=tick)
matplotlib.rc('ytick',  labelsize=tick)
matplotlib.rc('legend', fontsize=legend)


fig = plt.figure(figsize=(16/2.54, 12/2.54))
ax = fig.add_subplot(3,1,(2,3))

data = np.genfromtxt("data from lab/GROUP1_6.AES")
energy = data[:,0]
intensity = data[:,1]
intensity = intensity[energy < 600]
energy    = energy[energy < 600]
intensity -= np.min(intensity)
intensity /= np.max(intensity)

# Find positions of peaks
element_data = {
    "Mo": {"peak_range": (165, 185), "sensitivity": (0.33, 0.03), "transition": "MNN"},
    "C":  {"peak_range": (230, 270), "sensitivity": (0.20, 0.04), "transition": "KLL"}, 
    "N":  {"peak_range": (350, 375), "sensitivity": (0.32, 0.07), "transition": "KLL"},
    "O":  {"peak_range": (470, 500), "sensitivity": (0.49, 0.04), "transition": "KLL"}
}

# Get data for peak
for element in element_data:
    
    # Data from measurement 
    start, end = element_data[element]["peak_range"]
    indices = (energy > start) & (energy < end)
    element_data[element]["min_energy"]     = energy[indices][np.argmin(intensity[indices])]
    element_data[element]["max_energy"]     = energy[indices][np.argmax(intensity[indices])]
    element_data[element]["min_intensity"]  = np.min(intensity[indices])
    element_data[element]["max_intensity"]  = np.max(intensity[indices])
    element_data[element]["peak_center"]    = 0.5 * (element_data[element]["min_energy"] + element_data[element]["max_energy"])
    element_data[element]["peak_amplitude"] = element_data[element]["max_intensity"] - element_data[element]["min_intensity"]
    element_data[element]["corrected_amplitude"] = element_data[element]["peak_amplitude"] / element_data[element]["sensitivity"][0]

    # Data from reference
    csv = CSV.read("AES_reference/csv/" + element + ".csv")
    energy_peak = np.array(csv["Energy (eV)"])
    intensity_peak = np.array(csv["Intensity (a.u.)"])
    element_data[element]["reference_energy"] = energy_peak
    element_data[element]["reference_intensity"] = intensity_peak
    element_data[element]["reference_center"] = 0.5 * (energy_peak[np.argmin(intensity_peak)] + energy_peak[np.argmax(intensity_peak)])

# Calculate total amplitude
total_amplitude = 0
for element in element_data:
    total_amplitude +=  element_data[element]["corrected_amplitude"]

# Calculate relative composition
for element in element_data:
    element_data[element]["fraction"] = element_data[element]["corrected_amplitude"] / total_amplitude

# Calculate error margin in fraction
for element in element_data: 
    rel_sensitivity_error = element_data[element]["sensitivity"][1] / element_data[element]["sensitivity"][0]
    amplitude = element_data[element]["corrected_amplitude"]
    relative_variance = rel_sensitivity_error**2 * (total_amplitude - amplitude)**2 / total_amplitude**2
    
    for other_element in element_data:
        if (other_element == element): continue
        rel_sensitivity_error = element_data[other_element]["sensitivity"][1] / element_data[other_element]["sensitivity"][0]
        amplitude = element_data[element]["corrected_amplitude"]
        relative_variance += rel_sensitivity_error**2 * amplitude**2 / total_amplitude**2

    relative_error = np.sqrt(relative_variance)
    element_data[element]["fraction_error"] = element_data[element]["fraction"] * relative_error

for element in element_data: 
    print(f"{element}: ({100*element_data[element]['fraction']:.2f} ± {100*element_data[element]['fraction_error']:.2f})%")

# Calibrate energy
center_measured = np.array([element_data[element]["peak_center"] for element in element_data])
center_expected = np.array([element_data[element]["reference_center"] for element in element_data])
param = np.polyfit(center_measured, center_expected, 1)
energy_calibrated = np.polyval(param, energy)

# Plot calibration line
energy_measured = np.linspace(0, 650)
energy_corrected = np.polyval(param, energy_measured)
#ax.plot(center_measured, center_expected, "o")
#ax.plot(energy_measured, energy_corrected, "k--")
#ax.show()

# Plot reference curves
for i, element in enumerate(element_data):

    energy_peak = element_data[element]["reference_energy"]
    intensity_peak = element_data[element]["reference_intensity"]

    shift = [0.1, 0.6, 0.1, 0.1][i]
    amp = element_data[element]["peak_amplitude"]
    low = element_data[element]["min_intensity"]
    
    label = None
    if not i: label = "Reference"
    ax.plot(           energy_peak, shift + low + amp + amp * intensity_peak, "k", linewidth=1, label=label)
    ax.text(np.max(energy_peak)+10, shift + low + amp + amp * intensity_peak[-1] - 0.03, element, fontsize=legend)

# Plot measured spectrum
ax.plot(energy_calibrated, intensity, "k", linewidth=2, label="Measured")

# Also plot reference curves for some elements with no visible peaks
for i, element in enumerate(["Ar"]):
    csv = CSV.read("AES_reference/csv/" + element + ".csv")
    energy_peak = np.array(csv["Energy (eV)"])
    intensity_peak = np.array(csv["Intensity (a.u.)"])

    label = None
    if not i: label = "Reference"
    ax.plot(energy_peak, 0.2 * intensity_peak-0.1, "k--", linewidth=1, label=label)
    ax.text(np.max(energy_peak)+10, 0.2 * intensity_peak[-1]-0.13, element, fontsize=legend)

# Plot amplitude indicators
for element in element_data:
    max_energy    = np.polyval(param, element_data[element]["max_energy"])
    max_intensity = element_data[element]["max_intensity"]
    min_energy    = np.polyval(param, element_data[element]["min_energy"])
    min_intensity = element_data[element]["min_intensity"]
    lim = max(min_energy, max_energy) + 40
    ax.hlines(max_intensity, max_energy, lim, "r", linestyles="--", linewidth=1)
    ax.hlines(min_intensity, min_energy, lim, "r", linestyles="--", linewidth=1)
    ax.vlines(lim-5, min_intensity, max_intensity, "r", linewidth=1)
    ax.text(lim-30, min_intensity-0.11, f'{element_data[element]["peak_amplitude"]:.2f}', color="r", fontsize=legend)

ax.set_xlabel("Energy / eV")
ax.set_ylabel("d$N$/d$E$")
ax.set_ylim(-0.3, 2.2)
ax.legend(framealpha=1, loc='upper left')
ax.tick_params(labelleft=False)
ax.tick_params(left=False)


# Plot full spectrum for reference
ax = fig.add_subplot(3,1,1)

data = np.genfromtxt("data from lab/GROUP1_D.AES")
energy = data[:,0]
energy_corrected = np.polyval(param, energy)
intensity = data[:,1]
intensity -= np.min(intensity)
intensity /= np.max(intensity)

ax.plot(energy_corrected, intensity, "k-")
ax.set_ylim(0.765, 0.825)
ax.tick_params(labelleft=False)
ax.tick_params(left=False)
ax.set_ylabel("d$N$/d$E$")

fig.tight_layout()

plt.savefig(os.path.abspath(os.path.dirname(__file__)) + '/Figures/TIF120_UHV-lab_AES-spectrum.pdf')
plt.show()