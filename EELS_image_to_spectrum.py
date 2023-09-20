import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import CSV_handler as CSV

# Script for extracting spectra from Figure 6 in the paper 
# "Optical properties of V, Ta, and Mo from 0.1 to 35 eV" by Minni et.al 1974

images = ["Minni_a",
          "Minni_b",
          "Minni_c"]

for image in images:

    img = Image.open("EELS_reference/image/" + image + '.png').getdata()
    width, height = img.size
    img_data = np.array(img)
    img_data = 1 - np.sum(img_data[:,0:3], axis=1)/(3*255)
    img_data.resize((height, width))

    # Start and endpoint of axis in Minnis paper
    start = 0.0
    end   = 80.0

    energy = np.linspace(start, end, width)
    value  = np.argmax(img_data, axis=0).astype(np.double)

    # Cut into solid line sections
    defined_indices = value != 0
    chunk_label = np.zeros_like(defined_indices, dtype=np.int)

    chunk_count = 0
    for i in range(1, len(defined_indices)):
        if (defined_indices[i] and not defined_indices[i-1]):
            chunk_count += 1
        if (defined_indices[i]): 
            chunk_label[i] = chunk_count

    for chunk in range(1, np.max(chunk_label)+1):
        chunk_indices = chunk_label == chunk
       
        value -= np.min(value)
        value /= np.max(value)
        value = 1-value

        CSV.print_arrays_to_CSV(
        "EELS_reference/csv/" + image + "_" + str(chunk) + ".csv", 
        "Energy (eV)", energy[chunk_indices], 
        "Intensity (a.u.)", value[chunk_indices])
    
