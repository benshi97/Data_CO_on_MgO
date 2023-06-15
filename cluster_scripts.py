import pandas as pd

import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 9})
import extrapolate
from scipy.stats import linregress
from ase import units

Hartree = 27.211386245988
Bohr = 0.5291772105638411
cm1_to_eV = 1/8065.54429
hundredcm1 = 100 * cm1_to_eV * 1000

color_dict = {'red':'#e6194b',
'green': '#3cb44b',
'yellow': '#ffe119',
'blue': '#4363d8',
'orange': '#f58231',
'purple': '#911eb4',
'cyan':  '#42d4f4',
'magenta': '#f032e6',
'lime':  '#bfef45',
'pink': '#fabed4',
'teal': '#469990',
'lavendar': '#dcbeff',
'brown': '#9A6324',
'beige':'#fffac8',
'maroon':'#800000',
'mint': '#aaffc3',
'olive': '#808000',
'apricot':'#ffd8b1', 
'navy':'#000075',
'grey': '#a9a9a9',
'white': '#ffffff', 
'black':'#000000'}

plt.rcParams["axes.prop_cycle"] = plt.cycler(color=['#4363d8', '#e6194B', '#3cb44b', '#f58231', '#ffe119', '#911eb4', '#42d4f4', \
'#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', \
'#000075', '#a9a9a9', '#ffffff', '#000000'])

def read_vib_freq(filename, lines=None):
    """Read vibrational frequencies.

    Returns list of real and list of imaginary frequencies."""
    freq = []
    i_freq = []

    with open(filename, 'r',encoding = "ISO-8859-1") as f:
        lines = f.readlines()

    for line in lines:
        data = line.split()
        if 'THz' in data:
            if 'f/i=' not in data:
                freq.append(float(data[-2]))
            else:
                i_freq.append(float(data[-2]))
    return freq, i_freq

def get_vibrational_energy_contribution(vib_energies, temperature):
    """Calculates the change in internal energy due to vibrations from
    0K to the specified temperature for a set of vibrations given in
    eV and a temperature given in Kelvin. Returns the energy change
    in eV."""
    kT = units.kB * temperature
    dU = 0.
    for energy in vib_energies:
        dummy_ene = energy*0.001
        dU += dummy_ene / (np.exp(dummy_ene / kT) - 1.)
        # print(dummy_ene / (np.exp(dummy_ene / kT) - 1.))
    return dU

def get_ZPE_correction(vib_energies):
    """Returns the zero-point vibrational energy correction in eV."""
    zpe = 0.
    for energy in vib_energies:
        dummy_ene = energy*0.001
        zpe += 0.5 * dummy_ene
    return zpe


def get_quasi_rrho(r_freq,i_freq,T):
    k = units.kB
    combined_freq = r_freq + [0.0001]*len(i_freq)
    kT = (k*T*1000)

    dU = 0
    zpe = 0.
    eth = 0.
    for i in combined_freq:
        omega = 1/(1+((hundredcm1/i)**4))
        dURRho = i / (np.exp(i / kT) - 1.) + 0.5*i
        zpe += omega*0.5*i
        eth += omega*i / (np.exp(i / kT) - 1.) + (1 - omega)*0.5*kT
        dU += omega * dURRho + (1 - omega)*0.5*kT
        # print(omega,omega * dURRho + (1 - omega)*0.5*kT)

    return dU,eth,zpe,kT

def get_ad_rrho(r_freq,i_freq,dof,T):
    k = units.kB
    combined_freq = r_freq + [0.0001]*len(i_freq)
    kT = (k*T*1000)
    counter = 0

    dU = 0
    for i in combined_freq:
        if counter < dof:
            dU += i / (np.exp(i / kT) - 1.) + 0.5*i
        else:
            dU += 0.5*kT
        
        counter += 1
    
    return dU



def find_co_gamma(eads_list, tot_atom_list = [6, 22, 34, 42, 58, 82, 84, 100, 108]):
    power_list = []

    for i in range(1000):
        power_list += [np.round(0.1 + 0.01*i,3)]
    
    msd_list = []
    msd_best = []
    for k in range(1000):
        length = len(eads_list)
        slope, intercept, r, p, se = linregress([1/(x**(power_list[k])) for x in tot_atom_list[:length]],\
            [float(x) for x in eads_list[:length]])
        predicted_values = [intercept + slope/(x**(power_list[k])) for x in tot_atom_list[:length]]
        msd_list += [np.sqrt(np.mean(np.square(np.array(predicted_values) - np.array(eads_list[:length]))))] 
    msd_best = [power_list[msd_list.index(min(msd_list))], min(msd_list), msd_list.index(min(msd_list))]

    return msd_best[0]

def get_eads(filename,code_format='mrcc',typ='ccsdt',structs=['AD_SLAB','SLAB_CP','AD_CP']):
    if code_format == 'mrcc':
        a = find_energy(filename+ '/{0}/mrcc.out'.format(structs[0]),code_format=code_format,typ = typ)
        b = find_energy(filename+ '/{0}/mrcc.out'.format(structs[1]),code_format=code_format,typ = typ)
        c = find_energy(filename+ '/{0}/mrcc.out'.format(structs[2]),code_format=code_format,typ = typ)
        eads = (a-b-c)
    elif 'orca' in code_format:
        a = find_energy(filename+ '/{0}/orca.out'.format(structs[0]),code_format=code_format,typ = typ)
        b = find_energy(filename+ '/{0}/orca.out'.format(structs[1]),code_format=code_format,typ = typ)
        c = find_energy(filename+ '/{0}/orca.out'.format(structs[2]),code_format=code_format,typ = typ)
        eads = (a-b-c)
    elif 'vasp' in code_format:
        a = find_energy(filename+ '/{0}/OUTCAR'.format(structs[0]),code_format=code_format,typ = typ)
        b = find_energy(filename+ '/{0}/OUTCAR'.format(structs[1]),code_format=code_format,typ = typ)
        c = find_energy(filename+ '/{0}/OUTCAR'.format(structs[2]),code_format=code_format,typ = typ)
        eads = (a-b-c)        
    return eads  


def find_energy(filename,typ='ccsdt',code_format='mrcc'):

    if code_format=='mrcc':
        if typ == 'lccsdt':
            search_word = 'CCSD(T) correlation energy + MP2 corrections [au]:'
        elif typ == 'lccsdt_tot':
            search_word = 'Total LNO-CCSD(T) energy with MP2 corrections [au]'
        elif typ == 'lccsdt_lmp2_tot':
            search_word = 'Total LMP2 energy [au]'
        elif typ == 'ccsdt_mp2_tot':
            search_word = 'Total MP2 energy [au]'
        elif typ == 'ccsdt':
            search_word = 'CCSD(T) correlation energy [au]:'
        elif typ == 'ccsdt_tot':
            search_word = 'Total CCSD(T) energy'
        elif typ == 'hf':
            search_word = 'Reference energy [au]:    '
        elif typ == 'lmp2':
            search_word = 'LMP2 correlation energy [au]:         '
        elif typ == 'lmp2_tot':
            search_word = 'DF-MP2 energy [au]:       '
        elif typ == 'lmp2_corr':
            search_word = 'DF-MP2 correlation energy [au]:   '
        elif typ == 'mp2':
            search_word = 'MP2 correlation energy [au]:   '
        elif typ == 'lccsd':
            search_word = 'CCSD correlation energy + 0.5 MP2 corrections [au]:'
        elif typ == 'lccsdt_lccsd_tot':
            search_word = 'Total LNO-CCSD energy with MP2 corrections [au]:'
        elif typ == 'fnoccsdt_tot':
            search_word = 'Total CCSD(T+) energy + MP2 + PPL corr. [au]'
        elif typ == 'fnoccsd_tot':
            search_word = 'Total CCSD energy + MP2 + PPL corr. [au]:'
        elif typ == 'fnoccsdt_mp2_tot':
            search_word = 'DF-MP2 energy [au]:'
        elif typ == 'fnoccsdt':
            search_word = 'CCSD(T+) correlation en. + MP2 + PPL corr. [au]:'
        elif typ == 'fnoccsd':
            search_word = 'CCSD correlation energy + MP2 + PPL corr. [au]:'
        elif typ == 'fnoccsdt_mp2':
            search_word = 'DF-MP2 correlation energy'

        elif typ == 'ccsd':
            search_word = 'CCSD correlation energy [au]: '
        elif typ == 'ccsd_tot':
            search_word = 'Total CCSD energy [au]: '
        elif typ == 'dft':
            search_word = '***FINAL KOHN-SHAM ENERGY:'
        elif typ == 'B2PLYP':
            search_word = 'MP2 contribution [au]:'
        elif typ == 'DSDPBEP86':
            search_word = 'SCS-MP2 contribution [au]:'

        
        
        with open(filename, "r") as fp:
            a = [line for line in fp if search_word in line]
        if len(a) == 0:
            return 0.0
        else:
            if typ == 'dft':
                return float(a[-1].split()[-2])
            elif 'fnoccsdt_mp2' in typ:
                return float(a[0].split()[-1])

            else:
                return float(a[-1].split()[-1])
    
    elif code_format=='vasp':
        search_word= 'energy  without entropy='
        with open(filename, "r",encoding = "ISO-8859-1") as fp:
            a = [line for line in fp if search_word in line]
        if len(a) == 0:
            return 0.0
        else:
            return float(a[-1].split()[-1])
    
    elif code_format == 'quantum_espresso':
        search_word = '!    total energy              ='
        with open(filename, "r",encoding = "ISO-8859-1") as fp:
            a = [line for line in fp if search_word in line]
        if len(a) == 0:
            return 0.0
        else:
            return float(a[-1].split()[-2])
        
    elif code_format == 'cc4s':
        if typ == 'CCSD corr':
            search_word = 'Ccsd correlation energy:'
        elif typ == 'CCSD FS':
            search_word = 'Finite-size energy correction:'
        elif typ == 'CCSD BSIE':
            search_word = 'Ccsd-Bsie energy correction:'
        elif typ == 'HF':
            search_word = 'energy  without entropy='
        elif typ == '(T) corr':
            search_word = '(T) correlation energy:'
        elif typ == 'MP2 FS':
            search_word = 'Finite-size energy correction:'
        elif typ == 'MP2 corr':
            search_word = 'converged values  '

        with open(filename, "r",encoding = "ISO-8859-1") as fp:
            a = [line for line in fp if search_word in line]
        if len(a) == 0:
            return 0.0
        else:
            return float(a[-1].split()[-1])

    elif code_format=='vasp_wodisp':
        search_word= 'energy without entropy ='
        with open(filename, "r",encoding = "ISO-8859-1") as fp:
            a = [line for line in fp if search_word in line]
        if len(a) == 0:
            return 0.0
        else:
            return float(a[-1].split()[-1])

    elif code_format=='dftd3':
        search_word= ' Edisp /kcal,au'
        with open(filename, "r",encoding = "ISO-8859-1") as fp:
            a = [line for line in fp if search_word in line]
        if len(a) == 0:
            return 0.0
        else:
            return float(a[-1].split()[-1])        

    elif code_format=='orca':
        if typ == 'lccsdt':
            search_word = 'Final correlation energy'
        elif typ == 'mp2':
            search_word = 'E(MP2)'
        elif typ == 'ccsd':
            search_word = 'E(CORR)'
        elif typ == 'ccsdt':
            search_word = 'Final correlation energy'
        elif typ == 'hf':
            search_word = 'E(0)'
        elif typ == 'hf_lmp2':
            search_word = 'Total Energy'
        elif typ == 'lmp2':
            search_word = 'E(SL-MP2) including corrections'
        elif typ == 'ri-mp2':
            search_word = 'RI-MP2 CORRELATION ENERGY'
        # elif typ == 'mp2':
        #     search_word = 'MP2 correlation energy [au]:   '
        elif typ == 'lccsd':
            search_word = 'E(CORR)(corrected)'
        elif typ == 'dft':
            search_word = 'FINAL SINGLE POINT ENERGY'
        # elif typ == 'ccsd':
        #     search_word = 'CCSD correlation energy [au]: '        
        with open(filename, "r") as fp:
            a = [line for line in fp if search_word in line]
        if len(a) == 0:
            return 0.0
        elif typ == 'hf_lmp2':
            return float(a[-1].split()[-4])
        elif typ == 'ri-mp2':
            return float(a[-1].split()[-2])            
        else:
            return float(a[-1].split()[-1])

    elif code_format=='orca_mp2':
        if typ == 'lccsdt':
            search_word = 'Final correlation energy'
        # elif typ == 'ccsdt':
        #     search_word = 'CCSD(T) correlation energy [au]:'
        elif typ == 'hf':
            search_word = 'Total energy after final integration'
        elif typ == 'lmp2':
            search_word = 'DLPNO-MP2 CORRELATION ENERGY'
        # elif typ == 'mp2':
        #     search_word = 'MP2 correlation energy [au]:   '
        elif typ == 'lccsd':
            search_word = 'E(CORR)(corrected)'
        # elif typ == 'ccsd':
        #     search_word = 'CCSD correlation energy [au]: '

        
        with open(filename, "r") as fp:
            a = [line for line in fp if search_word in line]
        if len(a) == 0:
            return 0.0
        else:
            return float(a[-1].split()[-2])

def get_energy(filepath,system='MgO',basis_list =['SVP','TZVPP', 'QZVPP', 'CVDZ','CVTZ','CVQZ','CV5Z','VDZ','VTZ','VQZ','V5Z'],code_format='mrcc'):
    if code_format=='mrcc':
        if system=='TiO2':
            basis_list = ['SVP','TZVPP', 'QZVPP', 'CVTZ','CVQZ','CV5Z','VDZ','VTZ','VQZ','V5Z']
            #['SVP','TZVPP', 'QZVPP', 'CVTZ','CVQZ','CV5Z','VDZ','VTZ','VQZ','V5Z']

        energies_dict = dict.fromkeys(basis_list)

        for i in energies_dict:
            energies_dict[i] = {'perfect': {'hf': 0.0, 'lmp2': 0.0,'lccsd': 0.0,'lccsdt': 0.0, 'total_lmp2': 0.0, 'total_lccsd': 0.0, 'total_lccsdt': 0.0},
            'defect': {'hf': 0.0, 'lmp2': 0.0,'lccsd': 0.0, 'lccsdt': 0.0, 'total_lmp2': 0.0, 'total_lccsd': 0.0, 'total_lccsdt': 0.0},
            'O': {'hf': 0.0, 'lmp2': 0.0,'lccsd': 0.0, 'lccsdt': 0.0, 'total_lmp2': 0.0, 'total_lccsd': 0.0, 'total_lccsdt': 0.0},
            'vac_energy': {'hf': 0.0, 'lmp2': 0.0,'lccsd': 0.0, 'lccsdt': 0.0, 'total_lmp2': 0.0, 'total_lccsd': 0.0, 'total_lccsdt': 0.0}}

        
        for index, i in enumerate(basis_list):
            for j in ['perfect','defect', 'O']:
                for k in ['hf','lmp2','lccsd','lccsdt']:
                    if j == 'O' and k == 'lccsdt':
                        energies_dict[i][j][k] = find_energy('{0}/{1}/{2}/mrcc.out'.format(filepath,i,j),typ='ccsdt')
                    elif j == 'O' and k == 'lmp2':
                        energies_dict[i][j][k] = find_energy('{0}/{1}/{2}/mrcc.out'.format(filepath,i,j),typ='mp2')
                    elif j == 'O' and k == 'lccsd':
                        energies_dict[i][j][k] = find_energy('{0}/{1}/{2}/mrcc.out'.format(filepath,i,j),typ='ccsd')
                    else:
                        energies_dict[i][j][k] = find_energy('{0}/{1}/{2}/mrcc.out'.format(filepath,i,j),typ=k)

                energies_dict[i][j]['total_lmp2'] = energies_dict[i][j]['hf'] + energies_dict[i][j]['lmp2']
                energies_dict[i][j]['total_lccsd'] = energies_dict[i][j]['hf'] + energies_dict[i][j]['lccsd']
                energies_dict[i][j]['total_lccsdt'] = energies_dict[i][j]['hf'] + energies_dict[i][j]['lccsdt']

            energies_dict[i]['vac_energy']['hf'] = (energies_dict[i]['defect']['hf'] + energies_dict[i]['O']['hf'] - energies_dict[i]['perfect']['hf'])*Hartree
            for l in ['lmp2','lccsd','lccsdt']:
                energies_dict[i]['vac_energy'][l] = (energies_dict[i]['defect'][l] \
                    + energies_dict[i]['O'][l] - energies_dict[i]['perfect'][l])*Hartree
                energies_dict[i]['vac_energy']['total_{0}'.format(l)] =  energies_dict[i]['vac_energy']['hf'] + energies_dict[i]['vac_energy'][l]

        return energies_dict

    elif code_format=='orca':
        if system=='TiO2':
            basis_list = ['SVP','TZVPP', 'QZVPP', 'CVTZ','CVQZ','VTZ','VQZ','SVPD','TZVPPD','QZVPPD']
            #['SVP','TZVPP', 'QZVPP', 'CVTZ','CVQZ','CV5Z','VDZ','VTZ','VQZ','V5Z']

        energies_dict = dict.fromkeys(basis_list)

        for i in energies_dict:
            energies_dict[i] = {'perfect': {'hf': 0.0, 'lmp2': 0.0,'lccsd': 0.0,'lccsdt': 0.0, 'total_lmp2': 0.0, 'total_lccsd': 0.0, 'total_lccsdt': 0.0},
            'defect': {'hf': 0.0, 'lmp2': 0.0,'lccsd': 0.0, 'lccsdt': 0.0, 'total_lmp2': 0.0, 'total_lccsd': 0.0, 'total_lccsdt': 0.0},
            'O': {'hf': 0.0, 'lmp2': 0.0,'lccsd': 0.0, 'lccsdt': 0.0, 'total_lmp2': 0.0, 'total_lccsd': 0.0, 'total_lccsdt': 0.0},
            'vac_energy': {'hf': 0.0, 'lmp2': 0.0,'lccsd': 0.0, 'lccsdt': 0.0, 'total_lmp2': 0.0, 'total_lccsd': 0.0, 'total_lccsdt': 0.0}}
        
        for index, i in enumerate(basis_list):
            for j in ['perfect','defect', 'O']:
                for k in ['hf','lmp2','lccsd','lccsdt']:
                    energies_dict[i][j][k] = find_energy('{0}/{1}/{2}/orca.out'.format(filepath,i,j),typ=k,code_format='orca')
                energies_dict[i][j]['total_lmp2'] = energies_dict[i][j]['hf'] + energies_dict[i][j]['lmp2']

                energies_dict[i][j]['total_lccsd'] = energies_dict[i][j]['hf'] + energies_dict[i][j]['lccsd']
                energies_dict[i][j]['total_lccsdt'] = energies_dict[i][j]['hf'] + energies_dict[i][j]['lccsdt']

            energies_dict[i]['vac_energy']['hf'] = (energies_dict[i]['defect']['hf'] + energies_dict[i]['O']['hf'] - energies_dict[i]['perfect']['hf'])*Hartree
            for l in ['lmp2','lccsd','lccsdt']:
                energies_dict[i]['vac_energy'][l] = (energies_dict[i]['defect'][l] \
                    + energies_dict[i]['O'][l] - energies_dict[i]['perfect'][l])*Hartree
                energies_dict[i]['vac_energy']['total_{0}'.format(l)] =  energies_dict[i]['vac_energy']['hf'] + energies_dict[i]['vac_energy'][l]
        return energies_dict

    if code_format=='mrcc_can':
        if system=='TiO2':
            basis_list = ['SVP','TZVPP', 'QZVPP', 'CVTZ','VTZ','VQZ']
            #['SVP','TZVPP', 'QZVPP', 'CVTZ','CVQZ','CV5Z','VDZ','VTZ','VQZ','V5Z']

        energies_dict = dict.fromkeys(basis_list)

        for i in energies_dict:
            energies_dict[i] = {'perfect': {'hf': 0.0, 'ccsdt': 0.0, 'total': 0.0},
            'defect': {'hf': 0.0, 'ccsdt': 0.0, 'total': 0.0},
            'O': {'hf': 0.0, 'ccsdt': 0.0, 'total': 0.0},
            'vac_energy': {'hf': 0.0, 'ccsdt': 0.0, 'total': 0.0}}

        
        for index, i in enumerate(basis_list):
            for j in ['perfect','defect', 'O']:
                for k in ['hf','ccsdt']:
                    energies_dict[i][j][k] = find_energy('{0}/{1}/{2}/mrcc.out'.format(filepath,i,j),typ=k)
                energies_dict[i][j]['total'] = energies_dict[i][j]['hf'] + energies_dict[i][j]['ccsdt']

            energies_dict[i]['vac_energy']['hf'] = (energies_dict[i]['defect']['hf'] + energies_dict[i]['O']['hf'] - energies_dict[i]['perfect']['hf'])*Hartree
            energies_dict[i]['vac_energy']['ccsdt'] = (energies_dict[i]['defect']['ccsdt'] + energies_dict[i]['O']['ccsdt'] - energies_dict[i]['perfect']['ccsdt'])*Hartree
            energies_dict[i]['vac_energy']['total'] =  energies_dict[i]['vac_energy']['hf'] + energies_dict[i]['vac_energy']['ccsdt']
        return energies_dict


# Scripts for calculating the cWFT corrections at the different levels of theory and systems

def get_corr_mrcc(filename,method):
    data_smaller_def2 = get_energy('{0}'.format(filename),basis_list=['TZVPP','QZVPP'])
    data_smaller_def2_cbs = extrapolate.get_cbs(data_smaller_def2['TZVPP']['vac_energy']['hf'],data_smaller_def2['TZVPP']['vac_energy'][method],\
        data_smaller_def2['QZVPP']['vac_energy']['hf'],data_smaller_def2['QZVPP']['vac_energy'][method],X=3,Y=4,family='def2',convert_Hartree=False ,shift=-5.21/2,output=False)

    data_smaller_cc = get_energy('{0}'.format(filename),basis_list=['CVTZ','CVQZ'])
    data_smaller_cc_cbs = extrapolate.get_cbs(data_smaller_cc['CVTZ']['vac_energy']['hf'],data_smaller_cc['CVTZ']['vac_energy'][method],\
            data_smaller_cc['CVQZ']['vac_energy']['hf'],data_smaller_cc['CVQZ']['vac_energy'][method],X=3,Y=4,family='mixcc',convert_Hartree=False ,shift=-5.21/2,output=False)
    # print(data_smaller_cc_cbs[2], data_smaller_def2_cbs[2],data_smaller_cc_cbs[2] - data_smaller_def2_cbs[2])
    return data_smaller_cc_cbs[2], data_smaller_def2_cbs[2],data_smaller_cc_cbs[2] - data_smaller_def2_cbs[2]

def get_corr_mrcc_b2plyp(folder):
    ene_vac_hf = []
    ene_vac_mp2 = []
    i='B2PLYP'
    for j in ['TZ','QZ']:
        ene_perfect = find_energy('{0}/CV{1}/perfect/mrcc.out'.format(folder,j),typ='hf')
        ene_defect = find_energy('{0}/CV{1}/defect/mrcc.out'.format(folder,j),typ='hf')
        ene_O = find_energy('{0}/CV{1}/O/mrcc.out'.format(folder,j),typ='hf')
        ene_vac_hf += [(ene_defect + ene_O - ene_perfect)*Hartree]
        ene_perfect = find_energy('{0}/CV{1}/perfect/mrcc.out'.format(folder,j),typ='{0}'.format(i))
        ene_defect = find_energy('{0}/CV{1}/defect/mrcc.out'.format(folder,j),typ='{0}'.format(i))
        ene_O = find_energy('{0}/CV{1}/O/mrcc.out'.format(folder,j),typ='{0}'.format(i))
        ene_vac_mp2 += [(ene_defect + ene_O - ene_perfect)*Hartree]
    # print(ene_vac_hf,ene_vac_mp2)
    data_smaller_cc_cbs = extrapolate.get_cbs(ene_vac_hf[0],ene_vac_mp2[0],ene_vac_hf[1],ene_vac_mp2[1],X=3,Y=4,family='mixcc',output=False)
    ene_vac_hf = []
    ene_vac_mp2 = []
    i='B2PLYP'
    for j in ['TZ','QZ']:
        ene_perfect = find_energy('{0}/{1}VPP/perfect/mrcc.out'.format(folder,j),typ='hf')
        ene_defect = find_energy('{0}/{1}VPP/defect/mrcc.out'.format(folder,j),typ='hf')
        ene_O = find_energy('{0}/{1}VPP/O/mrcc.out'.format(folder,j),typ='hf')
        ene_vac_hf += [(ene_defect + ene_O - ene_perfect)*Hartree]
        ene_perfect = find_energy('{0}/{1}VPP/perfect/mrcc.out'.format(folder,j),typ='{0}'.format(i))
        ene_defect = find_energy('{0}/{1}VPP/defect/mrcc.out'.format(folder,j),typ='{0}'.format(i))
        ene_O = find_energy('{0}/{1}VPP/O/mrcc.out'.format(folder,j),typ='{0}'.format(i))
        ene_vac_mp2 += [(ene_defect + ene_O - ene_perfect)*Hartree]
    # print(ene_vac_hf,ene_vac_mp2)

    data_smaller_def2_cbs = extrapolate.get_cbs(ene_vac_hf[0],ene_vac_mp2[0],ene_vac_hf[1],ene_vac_mp2[1],X=3,Y=4,family='def2',output=False)
    return (data_smaller_cc_cbs[2], data_smaller_def2_cbs[2],data_smaller_cc_cbs[2] - data_smaller_def2_cbs[2])    

def get_energy_extrapolation(energies_dict,system='MgO',output=True,code_format='mrcc'):

    X = 'SVP'
    Y = 'TZVPP'

    a,b,cbs1 = extrapolate.get_cbs(energies_dict['SVP']['vac_energy']['hf'],energies_dict['SVP']['vac_energy']['lccsdt'],\
    energies_dict['TZVPP']['vac_energy']['hf'],energies_dict['TZVPP']['vac_energy']['lccsdt'],X=2,Y=3,family='def2',convert_Hartree=False,shift=0.0,output=False)

    X = 'TZVPP'
    Y = 'QZVPP'

    a,b,cbs2 = extrapolate.get_cbs(energies_dict['TZVPP']['vac_energy']['hf'],energies_dict['TZVPP']['vac_energy']['lccsdt'],\
    energies_dict['QZVPP']['vac_energy']['hf'],energies_dict['QZVPP']['vac_energy']['lccsdt'],X=3,Y=4,family='def2' ,convert_Hartree=False,shift=0.0,output=False)

    X = 'VDZ'
    Y = 'VTZ'

    a,b,cbs3 = extrapolate.get_cbs(energies_dict['VDZ']['vac_energy']['hf'],energies_dict['VDZ']['vac_energy']['lccsdt'],\
    energies_dict['VTZ']['vac_energy']['hf'],energies_dict['VTZ']['vac_energy']['lccsdt'],X=2,Y=3,family='mixcc',convert_Hartree=False,shift=0.0,output=False)

    X = 'VTZ'
    Y = 'VQZ'

    a,b,cbs4 = extrapolate.get_cbs(energies_dict[X]['vac_energy']['hf'],energies_dict[X]['vac_energy']['lccsdt'],\
    energies_dict[Y]['vac_energy']['hf'],energies_dict[Y]['vac_energy']['lccsdt'],X=3,Y=4,family='mixcc',convert_Hartree=False,shift=0.0,output=False)
    
    X = 'VQZ'
    Y = 'V5Z'

    a,b,cbs5 = extrapolate.get_cbs(energies_dict[X]['vac_energy']['hf'],energies_dict[X]['vac_energy']['lccsdt'],\
    energies_dict[Y]['vac_energy']['hf'],energies_dict[Y]['vac_energy']['lccsdt'],X=4,Y=5,family='mixcc',convert_Hartree=False,shift=0.0,output=False)

    if system=='MgO':
        X = 'CVDZ'
        Y = 'CVTZ'

        a,b,cbs6 = extrapolate.get_cbs(energies_dict[X]['vac_energy']['hf'],energies_dict[X]['vac_energy']['lccsdt'],\
        energies_dict[Y]['vac_energy']['hf'],energies_dict[Y]['vac_energy']['lccsdt'],X=2,Y=3,family='mixcc' ,convert_Hartree=False,shift=0.0,output=False)


    X = 'CVTZ'
    Y = 'CVQZ'

    a,b,cbs7 = extrapolate.get_cbs(energies_dict[X]['vac_energy']['hf'],energies_dict[X]['vac_energy']['lccsdt'],\
    energies_dict[Y]['vac_energy']['hf'],energies_dict[Y]['vac_energy']['lccsdt'],X=3,Y=4,family='mixcc',convert_Hartree=False ,shift=0.0,output=False)

    X = 'CVQZ'
    Y = 'CV5Z'

    a,b,cbs8 = extrapolate.get_cbs(energies_dict[X]['vac_energy']['hf'],energies_dict[X]['vac_energy']['lccsdt'],\
    energies_dict[Y]['vac_energy']['hf'],energies_dict[Y]['vac_energy']['lccsdt'],X=4,Y=5,family='mixcc',convert_Hartree=False ,shift=0.0,output=False)

    if output==True:
        print('{0:12s}: {1:12.7f}'.format('SVP/TZVPP',cbs1))
        print('{0:12s}: {1:12.7f}'.format('TZVPP/QZVPP',cbs2))
        print('{0:12s}: {1:12.7f}'.format('VDZ/VTZ',cbs3))
        print('{0:12s}: {1:12.7f}'.format('VTZ/VQZ',cbs4))
        print('{0:12s}: {1:12.7f}'.format('VQZ/V5Z',cbs5))

        if system=='MgO':
            print('{0:12s}: {1:12.7f}'.format('CVDZ/CVTZ',cbs6))
        
        print('{0:12s}: {1:12.7f}'.format('CVTZ/CVQZ',cbs7))
        print('{0:12s}: {1:12.7f}'.format('CVQZ/CV5Z',cbs8))

    if system=='MgO':
        return [cbs1,cbs2,cbs3,cbs4,cbs5,cbs6,cbs7,cbs8],['SVP/TZVPP','TZVPP/QZVPP','VDZ/VTZ','VTZ/VQZ','VQZ/V5Z','CVDZ/CVTZ','CVTZ/CVQZ','CVQZ/CV5Z']
    elif system=='TiO2' and code_format=='mrcc':
        return [cbs1,cbs2,cbs3,cbs4,cbs5,cbs7,cbs8],['SVP/TZVPP','TZVPP/QZVPP','VDZ/VTZ','VTZ/VQZ','VQZ/V5Z','CVTZ/CVQZ','CVQZ/CV5Z']


# Scripts for parsing energy from MRCC and ORCA DFT calculations
def get_dft_vac_energy(folder):
    perfect_ene = get_mrcc_ene('{0}/perfect/mrcc.out'.format(folder))
    defect_ene = get_mrcc_ene('{0}/defect/mrcc.out'.format(folder))
    O_ene = get_mrcc_ene('{0}/O/mrcc.out'.format(folder))
    vac_ene = (defect_ene + O_ene - perfect_ene)*Hartree
    return vac_ene

def get_dft_vac_energy_orca(folder):
    perfect_ene = get_orca_ene('{0}/perfect/orca.out'.format(folder))
    defect_ene = get_orca_ene('{0}/defect/orca.out'.format(folder))
    O_ene = get_mrcc_ene('{0}/O/mrcc.out'.format(folder))
    vac_ene = (defect_ene + O_ene - perfect_ene)*Hartree
    return vac_ene

def get_dft_vac_energy_orca_mgo_bulk(folder):
    perfect_ene = get_orca_ene('{0}/perfect/orca.out'.format(folder))
    defect_ene = get_orca_ene('{0}/defect/orca.out'.format(folder))
    O_ene = get_orca_ene('{0}/O/orca.out'.format(folder))
    vac_ene = (defect_ene + O_ene - perfect_ene)*Hartree
    return vac_ene


def get_mrcc_ene(filename):
    with open(filename) as f:
        a = []
        for line in f.readlines():
                if 'FINAL KOHN' in line:
                    a += [float(line.split()[3])]
    return a[-1]

def get_orca_ene(filename):
    with open(filename) as f:
        a = []
        for line in f.readlines():
                if 'FINAL' in line:
                    a += [float(line.split()[-1])]
    return a [-1]  
