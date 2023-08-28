import numpy as np
import CSV_handler as CSV


# Constants #
h = 4.135e-15   # Plancks constant [eV*s]
c = 3e8         # Speed of light in vacuum [m/s]


# Read CSV #
data = CSV.read('Mo_Werner_2009_formatted.csv')
header = CSV.get_header(data)
wavelength = data[header[0]] # micrometer
n = data[header[1]]
k = data[header[2]]


# Convert #
epsilon_real = n**2 - k**2
epsilon_imaginary = 2*n*k
eV = h*c/(wavelength*1e-6) # [eV]


# Polarizability method (Barnes and Le Ru) #


# Energy loss function (Weaver) #


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