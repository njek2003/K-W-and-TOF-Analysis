
"""
@author: Noah Kalicki
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
# Set global font properties
mpl.rcParams['font.family'] = 'serif'  # Use serif fonts
mpl.rcParams['font.serif'] = 'Cambria'  # Specific serif font
mpl.rcParams['font.size'] = 12  # Set default font size
mpl.rcParams['axes.titlesize'] = 20  # Font size for axis titles
mpl.rcParams['axes.labelsize'] = 16  # Font size for axis labels
mpl.rcParams['xtick.labelsize'] = 15  # Font size for x-tick labels
mpl.rcParams['ytick.labelsize'] = 15  # Font size for y-tick labels
mpl.rcParams['legend.fontsize'] = 15  # Font size for legend
import numpy as np
import pandas as pd
import os
from scipy.optimize import curve_fit
from scipy.stats import maxwell
from scipy.constants import k
import sys

#Check other Drafts of this code for comments
#-------------------#
## Insert your directories here ##
directory1 = input("Please paste your excel file directory: ") #Excel File Directory
directory2 = input("Please paste your TOF file directory: ") #TOF File(s) Directory
directory3 = input("Please paste your K&W file directory: ") #KW File(s) Directory
excel_file_name = input("Please paste your excel file *name*: ") #'0724_TOF_KW'
sheet1 = input("TOF file sheet name? ")
sheet2 = input("K&W file sheet name? ") 
print('\n')

excel_file = f"{excel_file_name}.xlsx"
filepath = rf"{directory1}\{excel_file}"
df1 = pd.read_excel(filepath, sheet_name = f'{sheet1}')
df2 = pd.read_excel(filepath, sheet_name = f'{sheet2}')
columns_df1 = df1.columns
columns_df2 = df2.columns
f1 = np.array(df1[columns_df1[0]]) #for TOF file name
f2 = np.array(df2[columns_df2[0]]) #for KW file name 
x = np.array(df1[columns_df1[3]])  #for molecule name
x2 = np.array(df1[columns_df1[2]]) #for nozzle temperature during TOF
x3 = np.array(df2[columns_df2[1]]) #for nozzle temperature during KW
x4 = np.array(df1[columns_df1[4]]) if len(columns_df1)-1 > 3 else None #for backing P during TOF
x5 = np.array(df2[columns_df2[2]]) if len(columns_df1)-1 > 3 else None #fir backing P during KW
x6 = np.array(df1[columns_df1[1]]) #for molar mass m/z

sticking = []
tfinal = []
vfinal = []
efinal = []
indcheck = []
ind = []
Temp = []
unc = []
key = []
mm=[]

for s in range(len(f2)):
    ## Quick note: This next for loop is about to say "OK... for each KW file (or multiple if needed),
    # we're going to correspond to it exactly ONE molecule TOF file, which has an index 
    # in the excel file data you provided. You can check which KW file index it was corresponded to
    # by just appending s as well as i into the list... where s is the index of the KW file 
    # and i is the index of the TOF file."
    len1 = len(key)
    for i in range(len(x)):
        if x[i] == 'CO2' and np.abs(x2[i] - x3[s]) <=15:
            if x4 is not None and x5 is not None:
                if np.abs(x4[i] - x5[s]) <= 15:
                    key.append((s,i))
            else:
                key.append((s,i))
    j=0
    while j < (len(key)-1):
        if key[j][0] == key[j+1][0]:
            key.pop(j) #Can be either j or j+1
        else:
            j += 1
    len2 = len(key)
    if len1 == len2:
        continue #This means that key had nothing new appended to it, so it would be useless to continue with this iteration
    for p in key:
        if p[0] == s:
            file_index = p[1]
            file_path = rf"{directory2}\{f1[file_index]}.txt"
            break
    filename1 = os.path.basename(file_path)
    TOF = pd.read_csv(file_path, header=None, delim_whitespace=True, names=['Time', 'Counts'])
    time = (TOF['Time'] / 1000).to_numpy()
    counts = TOF['Counts'].to_numpy()
    if s == 0:
        print("\n")
    else:
        None
    print(f"THIS IS THE OUTPUT FOR TOF ANALYSIS OF {filename1} WHICH CORRESPONDS TO {x[file_index]}")
    print("\n")
    #Fit Counts vs. time to a Maxwell-Boltzmann Distribution Model
    def maxwell_boltzmann(t, a, loc, scale, offset, slope):
        return a * maxwell.pdf(t, loc=loc, scale=scale) + offset + slope * t
        #A baseline correction was necessary to obtain a 'better looking' fit here.
    def main():
        global time_fine, counts_fine_fit
        try:
            a_guess = max(counts) #- min(counts)
            loc_guess = time[np.argmax(counts)]
            offset_guess = min(counts)
            scale_guess = (max(time)-min(time)) / 100 #for a guess of an arbitrary fraction of the time range
            slope_guess = 10                            #Feel free to change dividing factor based on fit quality...
            guesses = [a_guess,loc_guess,scale_guess,offset_guess,slope_guess]
            popt,pcov = curve_fit(maxwell_boltzmann, time, counts, p0=guesses)
            #perr= np.sqrt(np.diag(pcov))
            time_fine = np.linspace(min(time), max(time), 1000)
            counts_fine_fit = maxwell_boltzmann(time_fine, *popt)
            show_graphs = 'y'
            if show_graphs == 'y':
                plt.figure(figsize=(8, 6))
                plt.scatter(time, counts, label='Raw Data', color='k', s=15)
                plt.plot(time_fine, counts_fine_fit, '-', label='Maxwell-Boltzmann Fit', color = 'm')
                plt.xlabel('Time (sec)')
                plt.ylabel('Counts')
                plt.title(f'Time-of-Flight Data Fitted to Maxwell-Boltzmann Distribution         {filename1}')
                plt.legend()
                plt.grid(True, color = 'k', alpha = 0.2)
                plt.minorticks_on()
                plt.show()
            elif show_graphs == 'n':
                print("wow. Francisco worked so hard on these")
                print("\n")
            else:
                raise ValueError("Invalid input. Please enter 'y' or 'n'.")
        except ValueError as ve:
            print(f"ValueError: {ve}")
            sys.exit(0)
    if __name__ == "__main__":
          main()
    """
    '''Velocity Transform'''
    """
    distance = 1.162 #[m]
    velocity_fine = np.zeros(len(time_fine))
    for i in range(len(time_fine)):
        if time_fine[i] == 0:
            velocity_fine[0] = 0
        else:
            velocity_fine[i] = (distance/time_fine[i])
    #-------
    v_counts_TOF = np.zeros(len(velocity_fine))
    for i in range(len(velocity_fine)):
        if velocity_fine[i] == 0:
            v_counts_TOF[0] = 0
        else:
            v_counts_TOF[i] = (counts_fine_fit[i]/((velocity_fine[i])**3)) #From Jacobian
          
    ## Define gradient of v_counts 
    v_counts_derivative = np.gradient(v_counts_TOF)
    cross =[] 
    for i in range(len(v_counts_TOF)):
        if v_counts_derivative[i] >= 0 and v_counts_derivative[i-1] <= 0: 
            cross.append(i)
    if cross:
        print(f"Crossover indices in v_prob_der: {cross}")
        for index in cross: 
            print(f"Index: {index}")
            print("Value at index (v_prob): ",'%0.2g' %v_counts_derivative[index])
            print("Value at index - 1: ",'%0.2g' %v_counts_derivative[index-1])
    else:
        raise ValueError("No Crossover indices found")
    
    #Make the cutoff and filter it -- remember that the arrays are "backwards" 
    velocity_redefine = velocity_fine[cross[0]-150:cross[0]] #Because we want this cutoff change to permanent.
    v_counts = v_counts_TOF[cross[0]-150:cross[0]] #Because we want this cutoff change to permanent.
    #I made the redefinition dependent on the cross[0] so that a modification of this would (hopefully) not be necessary.
    
    
    ## Fit velocity to theoretical velocity distribution of a supersonic beam
    ### Define velocity function
    def ss_velocity_redefine(v, v_0, alpha, c):
        return c * (v**3) * np.exp(-(v-v_0)**2 / (alpha**2)) #From book.
    
        return c * (v**3) * np.exp(-(v-v_0)**2 / (alpha**2)) #From book.
    for p in key:
        if p[0] == s:
            if x[p[1]] == 'CO2':
                molar_mass = x6[p[1]]
                if molar_mass == 44:
                    molecule = f'$^{{12}}$CO$_{x[p[1]][2]}$' 
                if molar_mass == 45:
                    molecule = f'$^{{13}}$CO$_{x[p[1]][2]}$'

                mm.append(molar_mass)
                Temp.append(x2[p[1]])
                indcheck.append((p[0])) #KW file to use
                ind.append((p[1])) #TOF file to use
                T = x2[p[1]] #Temperature of Nozzle
            elif x[p[1]] == 'He':
                molar_mass = 4
                print("He found during CO2 Analysis ")
                sys.exit(1)
 
    print("\n")
    print(f"MM: {molar_mass} g/mol")
    print(f"T: {T} K")
    
    ### Make educated guesses for what v_0 (horizontal shift), alpha, and c (scale) should be
    m = (molar_mass / 1000) / (6.02 * 10**23)
    print("\n")
    v_max_index = np.argmax(v_counts)
    v_0_guess = velocity_redefine[v_max_index]
    v_alpha_guess = np.sqrt(2 * k * T / m)
    v_initial_guess = [v_0_guess, v_alpha_guess, max(v_counts)]
    
    
    ### Fit the data and save counts as new adjusted value
    v_popt, v_pcov = curve_fit(ss_velocity_redefine, velocity_redefine, v_counts, p0=v_initial_guess, maxfev=500000)
    v_0, alpha, v_c = v_popt
    print("The velocity fit parameters are v_0 = %10.2E, alpha = %10.2E, c = %10.2E" % (v_0, alpha, v_c))
    v_fit_counts = np.array(ss_velocity_redefine(velocity_redefine, *v_popt))
        
    ## Normalize the data
    v_area = np.trapz(v_fit_counts, velocity_redefine)
    v_prob = v_counts / np.abs(v_area)
    v_fit_prob = v_fit_counts / np.abs(v_area)
    
    #Calculate R^2 (velocity)
    r_squaredv = 1 - ((np.sum((v_prob - v_fit_prob)**2))/np.sum((v_prob - np.mean(v_prob))**2))
    
    # Energy transform
    ## Convert velocity to energy using E = 1/2 mv^2 (J)
    energy = 0.5 * m * (velocity_redefine**2)
    
    ## Convert counts in velocity domain to counts in energy domain
    e_counts = np.zeros(len(velocity_redefine))
    for index in range(len(velocity_redefine)):
        if velocity_redefine[index] == 0:
            e_counts[index] = (0)
        else:
            e_counts[index] = (v_fit_counts[index] / (m * velocity_redefine[index]))
    
    
    ## Fit energy to theoretical distribution
    ### Define energy distribution function (see Francisco for how this was derived)
    def ss_energy(E, d, E_0, E_alpha):
        return d * E * np.exp(-(np.sqrt(2 * E / m) - E_0)**2 / E_alpha**2)
    
    ### Make educated guesses on parameters
    d_guess = v_c * 2 / m**2
    E_0_guess = v_0
    E_alpha_guess = alpha
    E_initial_guess = [d_guess, E_0_guess, E_alpha_guess]
    
    e_popt, e_pcov = curve_fit(ss_energy, energy, e_counts, p0=E_initial_guess, maxfev=500000)
    E_0, E_alpha, E_c = e_popt
    print("The energy fit parameters are E_0 = %10.2E, E_alpha = %10.2E, E_c = %10.2E \n" % (E_0, E_alpha, E_c))
    e_fit_counts = np.array(ss_energy(energy, *e_popt))
    
    ## Normalize the data
    E_area = np.trapz(e_fit_counts, energy)
    e_prob = e_counts / np.abs(E_area)
    e_fit_prob = e_fit_counts / np.abs(E_area)
    
    ## Convert to meV
    energy_meV = energy * 6.24 * 10**18 * 1000
    
    #Calculate R^2 (energy)
    r_squarede = 1 - ((np.sum((e_prob - e_fit_prob)**2))/np.sum((e_prob - np.mean(e_prob))**2))
    
    # Final calculations #
    ## Write a function to calculate FWHM and help visualize their location on a graph
    def fwhm(x, y):
        half_max = 0.5 * max(y)
        indices_greater = np.where(y >= half_max)[0]
        if len(indices_greater) < 2:
            return np.nan, np.nan, np.nan
        else:
            left_index = indices_greater[0]
            right_index = indices_greater[-1]
            fwhm = np.abs(x[right_index] - x[left_index])
        return fwhm, x[left_index], x[right_index]
    ## Store FWHM calculations
    v_FWHM, v_FWHM_left, v_FWHM_right = fwhm(velocity_redefine, v_fit_prob)
    e_FWHM, e_FWHM_left, e_FWHM_right = fwhm(energy_meV, e_fit_prob)
    
    ## Store max velocity and energy
    v_max_index = np.where(max(v_fit_prob) == v_fit_prob)
    v_max = velocity_redefine[v_max_index][0]
    e_max_index = np.where(max(e_fit_prob) == e_fit_prob)
    e_max = energy_meV[e_max_index][0]
    
    ## Store FWHM calculations
    v_FWHM, v_FWHM_left, v_FWHM_right = fwhm(velocity_redefine, v_fit_prob)
    e_FWHM, e_FWHM_left, e_FWHM_right = fwhm(energy_meV, e_fit_prob)
    
    ## Calculate delta V / V and delta E / E
    deltav_v = v_FWHM / v_max
    deltae_e = e_FWHM / e_max
    
    # Print final values 
    print('V (m/sec): %0.5g' %v_max)
    print('E (meV): %0.5g' %e_max)
    print('Velocity FWHM: %0.5g' %v_FWHM)
    print('Energy FWHM: %0.5g' %e_FWHM)
    print('delta V / V: %0.5g' %deltav_v)
    print('delta E / E: %0.5g' %deltae_e)
    
    showgraphs2 = 'y'#input('Want to see the graphs? (y/n) ')
    if showgraphs2 == 'y':
        for p in key:
            if p[0] == s:
                plt.figure(figsize=(8, 6))
                plt.scatter(velocity_redefine, v_prob, label='Data', color='black', s=15) #CHECK BACK AT THESE
                plt.plot(velocity_redefine, v_fit_prob, '--m', label = 'Fit')
                plt.axvline(v_FWHM_left, color='r', linestyle='--', label=f'FWHM Left: {v_FWHM_left:.2E}')
                plt.axvline(v_FWHM_right, color='r', linestyle='--', label=f'FWHM Right: {v_FWHM_right:.2E}')
                plt.xlabel('Velocity (m/s)')
                plt.ylabel('Probability')
                plt.text(velocity_redefine[15-len(velocity_redefine)],np.mean(v_fit_prob), f"$R^2$: {r_squaredv:.4f}", color = 'm')
                plt.text(velocity_redefine[15-len(velocity_redefine)],np.mean(v_fit_prob)*0.85, f"{x[p[1]]} and {T}", color = 'm')
                plt.title(f'Filtered Velocity Data versus Probability         {filename1}')
                plt.legend()
                plt.grid(True, color = 'k', alpha = 0.2)
                plt.minorticks_on()
                plt.show()
                
                plt.figure(figsize=(8, 6))
                plt.scatter(energy_meV, e_prob, label='Data', color='black', s=15)
                plt.plot(energy_meV, e_fit_prob, '--m', label = 'Fit')
                plt.axvline(e_FWHM_left, color='r', linestyle='--', label=f'FWHM Left: {e_FWHM_left:.2E}')
                plt.axvline(e_FWHM_right, color='r', linestyle='--', label=f'FWHM Right: {e_FWHM_left:.2E}')
                plt.xlabel('Energy (meV)')
                plt.ylabel('Probability')
                plt.text(energy_meV[30-len(energy_meV)],np.mean(e_fit_prob), f"$R^2$: {r_squarede: .4f}", color = 'm')
                plt.text(energy_meV[30-len(energy_meV)],np.mean(e_fit_prob)*0.85, f"{x[p[1]]} and {T}", color = 'm')
                plt.title(f'Filtered Energy Data versus Probability         {filename1}')
                plt.legend()
                plt.grid(True, color = 'k', alpha = 0.2)
                plt.minorticks_on()
                plt.show()
    else:
        plt.figure(figsize=(10, 6))
        plt.scatter(energy_meV, e_prob, label='Data', color='white', s=15)
        plt.axis('off')
        plt.text(np.mean(energy_meV),np.mean(e_fit_prob), "wow. Francisco worked so hard on these", color = 'gray', fontsize = '15')
    for p in key:
        if p[0] == s:
            file_path2 = rf'{directory3}\{f2[p[0]]}.txt'
    filename2 = os.path.basename(file_path2)
    # Assuming the data starts after 23 lines of text
    skip_rows = 23
    df = pd.read_csv(file_path2, skiprows=skip_rows, names=['Time(s)', 'Data'], delimiter=',', skipinitialspace=True)
    d1 = (df['Time(s)'] ).values
    d2 = (df['Data']).values
    

    # Conditions
    der = np.gradient(d2)
    p2 = []
    c = []
    p222 = []
    
    #Min-Max Normalize the data to help only look for changes in signal
    D2 = (d2 - np.min(d2))/(np.max(d2)-np.min(d2))
    
    for i in range(len(d1)): #On 07/25/24, found a versatility flaw that will hopefully be taken care of with these new conditions
            if d1[i] > 300 and d1[i] < 750 and D2[i] - D2[i-50] >= 0.5: #looking for sig. change in amplitude of signal
                p2.append(i)
            if d1[i] > 300 and d1[i] < 750 and D2[i] - D2[i-50] <= -0.5: #looking for sig. change in amplitude of signal
                p2.append(i)
            if d1[i] > 920 and d1[i] < 975: #Impose a hard cutoff.
                c.append(i)
    #Truncate the values of these lists so that we can manually put cutoffs of KW in the desired places
    p2 = [np.min(p2),np.max(p2)]
    c = [np.min(c)]
    d2 = d2[:c[0]] #To just cut off when closing the gate valve
    d1 = d1[:c[0]]

    
    
    ## Plotting
    p21 = d1[p2[0]-500] #Condition for end of p1
    for i in range(len(d1)): #Condition for beginning of p3
        if d1[i]>700 and d1[i] < 800 and (der[i] - der[i-50]) > 2e-10:
            p22 = d1[p2[1]+1500]
            ye1 = p2[1]+1500
            ye2  = None
            break
        else:
            p22 = d1[p2[1]+500]
            ye1 = None
            ye2 = p2[1]+500
    for i in range(len(d1)):
        if d1[i]>850 and d1[i] < 1000 and np.abs(der[i] - der[i-50]) < 2e-10:
            d1c = d1[c[0]-100]
            ye3 = c[0]-100
            ye4 = None
            break
        else:
            d1c = d1[c[0]-500]
            ye3 = None
            ye4 = c[0]-500
    # d1c = d1[c[0]-500]  #Condition for end of p3
    p24 = d1[p2[1]-100] #Condition for end of p2
    p23 = d1[p2[0]+750] #Condition for beginning of p2
    plt.figure(figsize=(8,6))
    plt.plot(d1,d2, color = 'black',  linestyle = '--')
    plt.axvline(p21, color ='r', linewidth = 1)
    plt.axvline(p22, color ='r', linewidth = 1)
    plt.axvline(d1c, color = 'r', linewidth = 1)
    plt.axvline(p24, color ='r', linewidth = 1)
    plt.axvline(p23, color ='r', linewidth = 1) # Shifted to account for the behavior in the data
    plt.title('King & Wells ')
    plt.text(np.max(d1)*0.65, np.max(d2)*1.10, f'{filename2} (KW) <-> {filename1} (TOF)', fontsize = 10)
    plt.xlabel('Time(s)')
    plt.ylabel('Pressure (Pa)')
    plt.minorticks_on()
    
    
    try1 = np.mean(d2[0:p2[0]-500])
    std1 = np.std(d2[0:p2[0]-500])
    try2 = np.mean(d2[p2[0]+750:p2[1]-100])
    std2 = np.std(d2[p2[0]+750:p2[1]-100])
    try3 = np.mean(d2[ye1 if ye1 is not None else ye2:ye3 if ye3 is not None else ye4])
    std3 = np.std(d2[ye1 if ye1 is not None else ye2:ye3 if ye3 is not None else ye4])
    
    plt.axhline(try1,color ='r', linewidth = 1)
    plt.axhline(try2,color ='r', linewidth = 1)
    plt.axhline(try3,color ='r', linewidth = 1)
    plt.show()
    
    #Calculate sticking probability
    stick = (try2 - try3)/(try2 - try1)
    unc_stick = np.sqrt((((((try2-try3)/((try2-try1)**2))**2)*(std1**2))) + ((((
        (try3-try1)/((try2-try1)))**2)*(std2**2))) + (((std3**2)/((try2-try1)**2))))
    print("\n")
    print(f"{filename2} sticking prob.: ","%0.3g" %stick)
    print("\n")
    #Here we want to append the information we care about for making a visualization of the sticking probability
    sticking.append(round(stick,5))
    unc.append(round(unc_stick,5))
    vfinal.append(v_max)
    efinal.append(e_max)
    if s == len(x3) - 1:
        print("END of LOOP")
        print('')
v = np.array(vfinal)
e = np.array(efinal)
st = np.array(sticking)
unce = np.array(unc)

plt.figure(figsize=(8,6))
plt.errorbar(Temp, st, color = 'blue', yerr=unce, xerr=unce, fmt='o', ecolor='goldenrod', capsize=5, label='Temperature with error')
plt.title('Sticking Probability vs. Nozzle Temperature')
plt.xlabel('Nozzle Temperature (K)')
plt.ylabel('Sticking Probability')
if len(key) == 1:
    plt.ylim([0.85,.985])
    plt.xlim([0.00,np.max(Temp)*2.5])
    plt.text(np.min(Temp)*0.10, 0.875, f"File Date: {filename1[0:2]}/{filename1[2:4]}/{filename1[4:6]}", fontsize = 15)
    plt.text(np.min(Temp)*0.10,0.885, f"{molecule} Sticking", fontsize = 15)
    plt.text(np.min(Temp)*0.10,0.868, f"Number of Sticking Points: {len(st)}", fontsize = 15)
else:
    plt.ylim([0.85,.985])
    plt.text(np.min(Temp),0.875, f"File Date: {filename1[0:2]}/{filename1[2:4]}/{filename1[4:6]}", fontsize = 15)
    plt.text(np.min(Temp),0.885, f"{molecule} Sticking", fontsize = 15)
    plt.text(np.min(Temp),0.868, f"Number of Sticking Points: {len(st)}", fontsize = 15)
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(1.25)
plt.legend(loc = 'lower left')
plt.minorticks_on()
ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1.15, top=True,right=True)
plt.show()

plt.figure(figsize=(8,6))
plt.errorbar(v, st,color = 'blue', yerr=unce, xerr=unce, fmt='o', ecolor='goldenrod', capsize=5, label='Velocity with error')
plt.title('Sticking Probability vs. Velocity')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Sticking Probability')
if len(key) == 1:
    plt.ylim([0.85,.985])
    plt.xlim([0.00,np.max(v)*2.5])
    plt.text(np.min(v)*0.10, 0.875, f"File Date: {filename1[0:2]}/{filename1[2:4]}/{filename1[4:6]}", fontsize = 15)
    plt.text(np.min(v)*0.10,0.885, f"{molecule} Sticking", fontsize = 15)
    plt.text(np.min(v)*0.10,0.868, f"Number of Sticking Points: {len(st)}", fontsize = 15)
else:
    plt.ylim([0.85,.985])
    plt.text(np.min(v),0.875, f"File Date: {filename1[0:2]}/{filename1[2:4]}/{filename1[4:6]}", fontsize = 15)
    plt.text(np.min(v),0.885, f"{molecule} Sticking", fontsize = 15)
    plt.text(np.min(v),0.868, f"Number of Sticking Points: {len(st)}", fontsize = 15)
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(1.25)
plt.legend(loc = 'lower left')
plt.minorticks_on()
ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1.15, top=True,right=True)
plt.show()

plt.figure(figsize=(8,6))
plt.errorbar(e, st, color = 'blue',yerr=unce, xerr=unce, fmt='o', ecolor='goldenrod', capsize=5, label='Energy with error')
plt.title('Sticking Probability vs. Energy')
plt.xlabel('Energy (meV)')
plt.ylabel('Sticking Probability')
if len(key) == 1:
    plt.ylim([0.85,.985])
    plt.xlim([0.00,np.max(e)*2.5])
    plt.text(np.min(e)*0.10, 0.875, f"File Date: {filename1[0:2]}/{filename1[2:4]}/{filename1[4:6]}", fontsize = 15)
    plt.text(np.min(e)*0.10,0.885, f"{molecule} Sticking", fontsize = 15)
    plt.text(np.min(e)*0.10,0.868, f"Number of Sticking Points: {len(st)}", fontsize = 15)
else:
    plt.ylim([0.85,.985])
    plt.text(np.min(e),0.875, f"File Date: {filename1[0:2]}/{filename1[2:4]}/{filename1[4:6]}", fontsize = 15)
    plt.text(np.min(e),0.885, f"{molecule} Sticking", fontsize = 15)
    plt.text(np.min(e),0.868, f"Number of Sticking Points: {len(st)}", fontsize = 15)
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(1.25)
plt.legend(loc = 'lower left')
plt.minorticks_on()
ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1.15, top=True,right=True)
plt.show()

if len(ind)>1 and any(ind[i] == ind[i-1] for i in range(len(ind))):
    print("Here, should expect multiple sticking points in one region.") 
    print('')
else:
    None
print("These error bars were calculated by assuming the standard deviation of the points within the range of P1, P2, and P3 was the respective uncertainty for P1, P2, and P3.")
    

    
    
























