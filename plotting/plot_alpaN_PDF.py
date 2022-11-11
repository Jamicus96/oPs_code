from ssl import ALERT_DESCRIPTION_CLOSE_NOTIFY
import numpy as np
import matplotlib.pyplot as plt

# Read in file info
file_original = '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/alpha-n_without_oPs/PDFs_alpaN_original.txt'
file_without_oPs = '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/oldPDF_without_oPs/PDFs_Positronium_alphaN13C_0.5-9.0.txt'
file_with_oPs = '/Users/jp643/Documents/Studies/PhD/Antinu/Positronium/Results/oldPDF_with_oPs/spectra/PDFs_Positronium_alphaN13C_0.5-9.0.txt'

f_or = open(file_original, 'r')
f_wo = open(file_without_oPs, 'r')
f_wi = open(file_with_oPs, 'r')

data_original = f_or.readlines()
data_without = f_wo.readlines()
data_with = f_wi.readlines()
f_or.close()
f_wo.close()
f_wi.close()

# Convert to arrays
times_or = np.array(data_original[0].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)
alpha1_or = np.array(data_original[4].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)
alpha2_or = np.array(data_original[5].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)
alpha3_or = np.array(data_original[6].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)

times_wo = np.array(data_without[0].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)
alpha1_wo = np.array(data_without[4].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)
alpha2_wo = np.array(data_without[5].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)
alpha3_wo = np.array(data_without[6].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)

times_wi = np.array(data_with[0].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)
alpha1_wi = np.array(data_with[4].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)
alpha2_wi = np.array(data_with[5].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)
alpha3_wi = np.array(data_with[6].replace('[', '').replace(']', '').replace('\n', '').split(', '), dtype=float)

# Plot
plt.plot(times_or, alpha1_or, '-', label='Original, E < 3.5MeV', color='r')
plt.plot(times_or, alpha2_or, '--', label='Original, 3.5MeV < E < 5.4MeV', color='r', )
plt.plot(times_or, alpha3_or, ':', label='Original, E > 5.4MeV', color='r')

plt.plot(times_wo, alpha1_wo, '-', label='Without o-Ps, E < 3.5MeV', color='b')
plt.plot(times_wo, alpha2_wo, '--', label='Without o-Ps, 3.5MeV < E < 5.4MeV', color='b')
plt.plot(times_wo, alpha3_wo, ':', label='Without o-Ps, E > 5.4MeV', color='b')

plt.plot(times_wi, alpha1_wi, '-', label='With o-Ps, E < 3.5MeV', color='g')
plt.plot(times_wi, alpha2_wi, '--', label='With o-Ps, 3.5MeV < E < 5.4MeV', color='g')
plt.plot(times_wi, alpha3_wi, ':', label='With o-Ps, 5.4MeV < E', color='g')

plt.legend(loc='best')
plt.xlabel(r'$t_{res}$ (ns)')
plt.ylabel('Probability Density')
plt.title('AlphaN Prompt Event PDFs')
plt.xlim(-30, 230)

plt.show()