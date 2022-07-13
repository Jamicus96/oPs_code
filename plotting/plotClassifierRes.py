import numpy as np
import matplotlib.pyplot as plt

# Read in file

oPs_file = open('ClassifierRes_oPs_output.txt', 'r')
oPs_lines = oPs_file.readlines()
oPs_delays = np.array(oPs_lines[0].replace('[', '', 1).replace(']', '', 1).split(', '))
oPs_delays = oPs_delays.astype(float)
oPs_logLdiff = np.array(oPs_lines[1].replace('[', '', 1).replace(']', '', 1).split(', '))
oPs_logLdiff = oPs_logLdiff.astype(float)
oPs_nhits = np.array(oPs_lines[2].replace('[', '', 1).replace(']', '', 1).split(', '))
oPs_nhits = oPs_nhits.astype(float)
oPs_file.close()

elec_file = open('ClassifierRes_e-_output.txt', 'r')
elec_lines = elec_file.readlines()
elec_delays = np.array(elec_lines[0].replace('[', '', 1).replace(']', '', 1).split(', '))
elec_delays = elec_delays.astype(float)
elec_logLdiff = np.array(elec_lines[1].replace('[', '', 1).replace(']', '', 1).split(', '))
elec_logLdiff = elec_logLdiff.astype(float)
elec_nhits = np.array(elec_lines[2].replace('[', '', 1).replace(']', '', 1).split(', '))
elec_nhits = elec_nhits.astype(float)
elec_file.close()

oPs_scaled_logLdiff = oPs_logLdiff / oPs_nhits
elec_scaled_logLdiff = elec_logLdiff / elec_nhits

############ LL vs delay #############

plt.scatter(oPs_delays, oPs_scaled_logLdiff, label='o-Ps')
plt.scatter(elec_delays, elec_scaled_logLdiff, label='e-')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('decay time (ns)')
plt.ylabel('log-Likelyhood difference divided by nhit')
plt.title('1 MeV Likelyhood o-Ps vs e+ For Different Decay Times')
plt.legend(loc='best')
#plt.xlim([-10, 100])
plt.ylim([-0.14, 0.065])

############ LL #############

# bins = np.arange(-0.1, 0.06, 0.001)
# plt.hist(x = oPs_scaled_logLdiff, bins = bins, color = 'b', alpha=0.5, rwidth=0.85, label = 'o-Ps')
# plt.hist(x = elec_scaled_logLdiff, bins = bins, color = 'r', alpha=0.5, rwidth=0.85, label = 'e-')
# plt.legend(loc='best')
# #plt.xlim([-50, 50])
# plt.xlabel('log-Likelyhood difference divided by nhit')

############ Receiver Operating Characteristic (ROC) curve #############

# N = 100
# logL_perNHit = np.linspace(-0.09, 0.05, N)
# signal = np.zeros(N)
# noise = np.zeros(N)

# for n in range(N):
#     cut = logL_perNHit[n]

#     # Ratio of signal (o-Ps) that passes the cut
#     for i in range(len(oPs_scaled_logLdiff)):
#         if oPs_scaled_logLdiff[i] > cut:
#             signal[n] += 1
#     signal[n] *= 100 / len(oPs_scaled_logLdiff)

#     # Ratio of noise (e-) that passes the cut
#     for j in range(len(elec_scaled_logLdiff)):
#         if elec_scaled_logLdiff[j] > cut:
#             noise[n] += 1
#     noise[n] *= 100 / len(elec_scaled_logLdiff)

# plt.plot(noise, signal)
# plt.xlabel('Percentage of noise (e-) that passes cut')
# plt.ylabel('Percentage of signal (o-Ps) that passes cut')
# plt.title('Positronium Classifier ROC curve')

######################################

plt.show()