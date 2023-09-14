import numpy as np
import matplotlib.pyplot as plt
import CSV_handler as CSV

reference_elements = ["Ag", "Ar", "Be", "C", "Ca", "Cr", "I", "Mo", "N", "O", "Pd", "Ti"]

data = np.genfromtxt("data from lab/GROUP1_6.AES")
energy = data[:,0]
intensity = data[:,1]
intensity -= np.min(intensity)
intensity /= np.max(intensity)
plt.plot(data[:,0], data[:,1], "k", linewidth=1)

for i, element in enumerate(reference_elements):
    csv = CSV.read("AES_reference/csv/" + element + ".csv")
    energy = np.array(csv["Energy (eV)"])
    intensity = np.array(csv["Intensity (a.u.)"])
    plt.plot(energy, 1.1 + 0.2 * (i + intensity), "k", linewidth=1)

#plt.ylim((-0.1, 3))
plt.xlim((0, 700))
plt.show()
