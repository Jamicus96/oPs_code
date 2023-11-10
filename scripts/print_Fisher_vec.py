import numpy as np
import matplotlib.pyplot as plt


def fisher_discriminant(x_1, x_2, N_1, N_2):

    # Compute means
    mu_1 = np.sum(x_1, axis=0) / len(x_1)
    mu_2 = np.sum(x_2, axis=0) / len(x_2)

    # Compute covariance matrices
    S_1 = np.cov(x_1.transpose())
    S_2 = np.cov(x_2.transpose())

    # Compute inverse of weigthed sum of covariance matrices
    S_inv = np.linalg.inv((float(N_1) / float(N_1 + N_2)) * S_1 + (float(N_2) / float(N_1 + N_2)) * S_2)

    # Compute projection vector
    a = np.matmul(S_inv, mu_1 - mu_2)

    # Work out if a.mu_1 returns higher or not
    higher = False
    if (np.sum(a * mu_1) > np.sum(a * mu_2)):
        higher = True
    
    # Return projection vector a, threshold t (half-way between projected means), and higher bool
    t = 0.5 * (np.sum(a * mu_1) + np.sum(a * mu_2))

    return a, t, higher

def read_file(file_address):

    f = open(file_address)
    lines = f.readlines()
    f.close()

    x = np.zeros((len(lines), len(lines[0].replace('\n', '').split(', '))))
    for i in range(len(lines)):
        line_splt = lines[i].replace('\n', '').split(', ')
        x[i] = np.array(line_splt, dtype=float)
    
    return x

def main():

    IBD_file = '/mnt/lustre/scratch/epp/jp643/antinu/Positronium/labppo_2p2_scintillator/flat/hists/hist_IBD_flat.txt'
    alphaN_file = '/mnt/lustre/scratch/epp/jp643/antinu/Positronium/labppo_2p2_scintillator/flat/hists/hist_alphaN_13C_flat.txt'
    print('Reading in:', IBD_file)
    x_IBD = read_file(IBD_file)
    print('Reading in:', alphaN_file)
    x_alphaN = read_file(alphaN_file)

    a, t, higher = fisher_discriminant(x_IBD, x_alphaN, 1, 1)

    a_str = '['
    for i in range(len(a) - 1):
        a_str += '%.4f, ' % a[i]
    a_str += '%.4f]' % a[-1]

    print('a =', a_str)
    print('t =', t)
    print('higher =', higher)
    

if __name__ == '__main__':
    main()