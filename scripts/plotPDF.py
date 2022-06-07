import numpy as np
import matplotlib.pyplot as plt

# Read in file
info = [('pdf_oPs_output.txt', 'o-Ps'),
        ('pdf_e+_output.txt', 'e+'),
        ('pdf_e-_output.txt', 'e-')]

for filename, particle in info:
    file = open(filename, 'r')
    lines = file.readlines()

    times = np.array(lines[0].replace('[', '', 1).replace(']', '', 1).split(', '))
    times = times.astype(float)
    probs = np.array(lines[1].replace('[', '', 1).replace(']', '', 1).split(', '))
    probs = probs.astype(float)

    file.close()
    plt.plot(times, probs, label=particle)

#plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$t_{res}$ (ns)')
plt.ylabel('probability density')
plt.title('1 MeV e+/e- Pulse Shapes (~10 000 events each)')
#plt.xlim([-10, 100])
plt.legend(loc='best')
plt.show()