import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import CSV_handler as CSV

images = ["I_400_600",
          "O_450_530",
          "Ag_200_400",
          "Ar_150_230",
          "Be_50_150",
          "C_200_300",
          "Ca_220_350",
          "Cr_250_600",
          "I_400_600",
          "Mo_80_230",
          "N_300_440",
          "Pd_200_400",
          "Ti_300_500"]

for image in images:


    img = Image.open("AES_reference/image/" + image + '.png').getdata()
    width, height = img.size
    img_data = np.array(img)
    img_data = 1 - np.sum(img_data[:,0:3], axis=1)/(3*255)
    img_data.resize((height, width))
    
    element, start, end = image.split("_")
    start = float(start)
    end   = float(end)

    energy = np.linspace(start, end, width)
    value  = np.argmax(img_data, axis=0).astype(np.double)

    value -= np.min(value)
    value /= np.max(value)

    #plt.imshow(img_data)
    #plt.plot(value, "k--")
    #plt.title(element)
    #plt.show()

    CSV.print_arrays_to_CSV(
        "AES_reference/csv/" + element + ".csv", 
        "Energy (eV)", energy, 
        "Intensity (a.u.)", value)