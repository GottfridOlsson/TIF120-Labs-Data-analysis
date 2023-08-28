import numpy as np
import CSV_handler as CSV

# Read CSV #
data = CSV.read('Mo_Werner_2009_formatted.csv')
header = CSV.get_header(data)
wavelength = data[header[0]] # micrometer
n = data[header[1]]
k = data[header[2]]

# Calculate dielectric function #
