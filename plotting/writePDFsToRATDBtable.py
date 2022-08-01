import json

json_file_address = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/Positronium/results/PDFs/PDF_stats.json'
ratdb_file_address = '/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/rat/data/POSITRONIUM.ratdb'


def rewriteTABLE():
    f = open(ratdb_file_address, "r")
    ratdb_lines = f.readlines()
    f.close()

    json_file = open(json_file_address)
    data = json.load(json_file)

    # Copy lines before probs
    pre_prob_lines = []
    for line in ratdb_lines:
        if 'oPs_recon_energies' in line:
            break
        else:
            pre_prob_lines.append(line)

    energies = [('0.5', '1.5'), ('1.0', '2.0'), ('2.5', '3.5'), ('3.0', '4.0'),\
                ('3.5', '4.5'), ('4.0', '5.0'), ('4.5', '5.5'), ('5.0', '6.0'),\
                ('5.5', '6.5'), ('6.0', '7.0'), ('6.5', '7.5'), ('7.0', '8.0'),\
                ('7.5', '8.5'), ('8.0', '9.0'), ('8.5', '9.5'), ('9.0', '10.0')]

    # Create lines with prob info to print to file
    oPs_recon_energies = '    oPs_recon_energies: ['
    elec_recon_energies = '    elec_recon_energies: ['
    oPs_lines = ['    oPs_probabilities: {\n']
    elec_lines = ['    electron_probabilities: {\n']
    for oPs_energy, elec_energy in energies:
        oPs_PDF = data['o-Ps'][oPs_energy]['PDF_Fitted']
        elec_PDF = data['e-'][elec_energy]['PDF_Fitted']

        oPs_str = '        "' + str(float(oPs_energy) + 1.0) + '": ['
        for prob in oPs_PDF:
            oPs_str += str(prob) + ', '
        oPs_str = oPs_str[:-2] + ']'

        elec_str = '        "' + elec_energy + '": ['
        for prob in elec_PDF:
            elec_str += str(prob) + ', '
        elec_str = elec_str[:-2] + ']'

        oPs_recon_energies += '"' + str(float(oPs_energy) + 1.0) + '"'
        elec_recon_energies += '"' + elec_energy + '"'
        if oPs_energy != energies[len(energies) - 1][0]:
            oPs_str += ','
            elec_str += ','
            oPs_recon_energies += ', '
            elec_recon_energies += ', '
        
        elec_lines.append(elec_str + '\n')
        oPs_lines.append(oPs_str + '\n')

    oPs_recon_energies += '],\n'
    elec_recon_energies += '],\n'
    oPs_lines.append('    },\n')
    elec_lines.append('    }\n')

    # Put it all together
    new_lines = pre_prob_lines + [oPs_recon_energies] + ['\n'] + [elec_recon_energies] + ['\n'] + oPs_lines + ['\n'] + elec_lines + ['  }\n']

    # Overwrite file
    new_ratdb_table = "".join(new_lines)
    f = open(ratdb_file_address, "w")
    f.write(new_ratdb_table)
    f.close()


if __name__ == '__main__':
    rewriteTABLE()