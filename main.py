import numpy as np
import CSV_handler as CSV
import matplotlib.pyplot as plt


# Constants #
h = 4.135e-15   # Plancks constant [eV*s]
c = 3e8         # Speed of light in vacuum [m/s]


# Read CSV #
data = CSV.read('Mo_Werner_2009_formatted_from_Canvas.csv')
header = CSV.get_header(data)
wavelength = data[header[0]] # micrometer
n = data[header[1]]
k = data[header[2]]


# Convert #
epsilon_real      = n**2 - k**2
epsilon_imaginary = 2*n*k
eV = h*c/(wavelength*1e-6) # [eV]

plt.plot(eV, epsilon_real, label="Re($\\epsilon$)")
plt.plot(eV, epsilon_imaginary, label="Im($\\epsilon$)")
plt.ylabel("Relative permittivity")
plt.xlabel("Energy (eV)")
plt.ylim(-10, 10)
plt.legend()
plt.show()


# Polarizability method (Barnes and Le Ru) #

# Bulk polarizability of thin slab according to Barnes
alpha_V_real = 1 - epsilon_real / (epsilon_real**2 + epsilon_imaginary**2)
alpha_V_imaginary = epsilon_imaginary / (epsilon_real**2 + epsilon_imaginary**2)

plt.plot(eV, alpha_V_real, label="Re($\\alpha_V$)")
plt.plot(eV, alpha_V_imaginary, label="Im($\\alpha_V$)")
plt.legend()
plt.xlabel("Energy (eV)")
plt.ylabel("Polarizability of thin slab (units of $\\epsilon_0$)")
plt.show()

# Polarizability of small sphere
alpha_sphere_real = ((epsilon_real - 1) * (epsilon_real + 2) + epsilon_imaginary**2) \
                    / (2 * epsilon_imaginary**2)
alpha_sphere_imaginary = ((epsilon_real + 2) * epsilon_imaginary - (epsilon_real - 1) * epsilon_imaginary) \
                        / ((epsilon_real + 2)**2 + epsilon_imaginary**2)
alpha_sphere_magnitude = np.sqrt(alpha_sphere_real**2 + alpha_sphere_imaginary**2)

plt.plot(eV, alpha_sphere_real, label="Re($\\alpha_{Sphere}$)")
plt.plot(eV, alpha_sphere_imaginary, label="Im($\\alpha_{Sphere}$)")
plt.plot(eV, alpha_sphere_magnitude, label="Mag($\\alpha_{Sphere}$)")
plt.legend()
plt.xlabel("Energy (eV)")
plt.ylabel("Polarizability of sphere (a.u.)")
plt.show()

# Energy loss function (Weaver) #
els_bulk = epsilon_imaginary / (epsilon_real**2 + epsilon_imaginary**2)
els_surface = epsilon_imaginary / ((1+epsilon_real)**2 + epsilon_imaginary**2)

plt.plot(eV, els_bulk, label="Energy loss function, bulk")
plt.plot(eV, els_surface, label="Energy loss function, surface")
plt.legend()
plt.xlabel("Energy (eV)")
plt.ylabel("Energy loss functinon (arbitrary unit)")
plt.show()

CSV.print_arrays_to_CSV("output/ELS_theoretical_bulk.csv", "Energy (eV)", eV, "Loss (a.u.)", els_bulk)
CSV.print_arrays_to_CSV("output/ELS_theoretical_surface.csv", "Energy (eV)", eV, "Loss (a.u.)", els_surface)
CSV.print_arrays_to_CSV("output/polariability_of_sphere.csv", 
                        "Energy (eV)", eV,
                        "Re(alpha) (a.u.)", alpha_sphere_real,
                        "Im(alpha) (a.u.)", alpha_sphere_imaginary,
                        "abs(alpha) (a.u.)", alpha_sphere_magnitude)

'''
Compare and discuss the bulk and surface plasmon energies for Mo measured by EELS
with the corresponding plasmon energies calculated on the basis of the complex
dielectric function of Mo (data is found on Canvas in the file Mo_Werner_2009). The
theoretical estimation, using the data of the dielectric function of Mo, should be done
as described in the papers by Barnes and Le Ru (polarizability) and more rigorously as
described by Weaver (energy loss function). Use both ways of estimating the plasmon
energies theoretically and compare these results with the experimental ones (both
your own results and the other experimental results presented in the preparatory reading)
'''