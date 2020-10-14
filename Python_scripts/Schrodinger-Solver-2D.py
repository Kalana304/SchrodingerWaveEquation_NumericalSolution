# Author - K.G. Abeywardena
# Date - 29/01/2020

import sys
import numpy as np 
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import matplotlib
from mpl_toolkits.mplot3d import Axes3D


plt.style.use('classic')


#Function for creating the Hamiltonian matrix
def Hamiltonian(N, dx, potential, h_bar, mass):          
    """
        This fucntion creates the Hamiltoninan matrix that is required to calculate the numerical solution
        for the Schrodinger's 2-D Wave Equation for a given 2-D potential curve. Assumes the Seperability
        of Variables. 

        Inputs:
            N           :   Number of data points to consider
            dx          :   Differnce between two consequetive x values
            potential   :   Potential function 
            h_bar       :   -(m*)/(2h)  where m* = equivalent mass and h = plank constant 
            mass        :   Mass of the particle 

        Outputs:
            norm vectors    :   Normalized eigen vectors of the Hamiltonian Matrix
            eigen energy    :   Eigen Values of the Hamiltonian Matrix
            probability     :   Probability Density Distribution 
    """                            

    Laplacian = (-2*np.diag(np.ones((N)),0) + np.diag(np.ones((N-1)),1) + np.diag(np.ones((N-1)),-1))/(dx**2)
    H = -(h_bar**2/(2*mass))*Laplacian + np.diag(potential)
    eigen_energy, eigen_vectors = eigh(H); print(eigen_energy.shape)

    probability = np.abs(eigen_vectors)**2
    K = np.sum(probability, axis=0)*dx
    norm_vectors = eigen_vectors/np.sqrt(K)             #Normalizing the eigen vectors
    probability /= K                                    #Normalizing the PDFs
    return norm_vectors, eigen_energy, probability

#Plotting Functions
def plotWavefunc(x, potential, num_states, wave_vectsX, wave_vectsY, energy, func_name):
    """
        This function is for the visualization of the Potential Energy Functions and 
        Corresponding Wave Functions which are the solutions of to the wave equation.

        Inputs:
            x           :   Data points
            potential   :   Potential Energy Function
            num_states  :   Number of energy primary states
            wave_vects  :   Wave function - as a set of vectors 
            energy      :   Energy values of each state 
            func_name   :   Function name of potential energy curve 

        Outputs:
            None
    """
    fig = plt.figure('WaveFunc', facecolor='w', edgecolor='k')
    
    for i in range(num_states-1, -1, -1):
        
        ax = fig.add_subplot(num_states//2, num_states//2, i+1,  projection='3d')

        X, Y = np.meshgrid(x, x)
    
        wave_vect = np.reshape(wave_vectsX[:,i], (len(wave_vectsX),1))*np.transpose(np.reshape(wave_vectsY[:,i], (len(wave_vectsX),1)))
        ax.plot_surface(X, Y, wave_vect + energy[i])
        
        plt.xlabel('Distance (x)', fontdict={'weight': 'bold', 'size': 8, 'color': 'black'})
        plt.ylabel('Distance (y)', fontdict={'weight': 'bold', 'size': 8, 'color': 'black'})
        ax.set_title('Energy {} =  {:10.4f}eV'.format(i,energy[i]), fontdict={'fontsize': 10, 'fontweight': 'bold'})
        plt.grid(b=True)
        plt.autoscale(tight=True)
        
    plt.suptitle(f'Wave Functions for {func_name}', fontsize = 22, fontweight = 'bold')       
    plt.savefig(f'wavePlot-{func_name}.png')
    plt.close('all')
    return

def plotProbfunc(x, potential, num_states, wave_vectsX, wave_vectsY, energy, func_name):
    """
        This function is for the visualization of the Potential Energy Functions and 
        Corresponding PDFs which are the solutions of to the wave equation.

        Inputs:
            x           :   Data points
            potential   :   Potential Energy Function
            num_states  :   Number of energy primary states
            prob_vects  :   PDF of each wave - as a set of vectors 
            energy      :   Energy values of each state 
            func_name   :   Function name of potential energy curve 

        Outputs:
            None
    """
    fig = plt.figure('PDF', facecolor='w', edgecolor='k')

    for i in range(num_states-1, -1, -1):

        ax = fig.add_subplot(num_states//2, num_states//2, i+1,  projection='3d')
        X, Y = np.meshgrid(x, x)

        wave_vect = np.reshape(wave_vectsX[:,i], (len(wave_vectsX),1))*np.transpose(np.reshape(wave_vectsY[:,i], (len(wave_vectsX),1)))
        probability = np.abs(wave_vect)**2
        
        ax.plot_surface(X, Y, probability + energy[i])

        plt.xlabel('Distance (x)', fontdict={'weight': 'bold', 'size': 8, 'color': 'black'})
        plt.ylabel('Distance (y)', fontdict={'weight': 'bold', 'size': 8, 'color': 'black'})
        ax.set_title('Energy {} =  {:10.4f}eV'.format(i,energy[i]), fontdict={'fontsize': 10, 'fontweight': 'bold'})
        plt.grid(b=True)
        plt.autoscale(tight=True)
        
    plt.suptitle(f'PDF Functions for {func_name}', fontsize = 22, fontweight = 'bold')       
    plt.savefig(f'PDFPlot-{func_name}.png')
    plt.close('all')
    return

def infiniteWell(N, dx, h_bar, mass, x, n_states):
    n = np.floor(N/3).astype('int32')
    V_x = np.ones((N)) * 2000
    V_x[n:2*n] = 0
    wave_vec_x, energy_x, _ = Hamiltonian(N,dx, V_x, h_bar, mass)   # Solution for x-direction
    wave_vec_y, energy_y, _ = Hamiltonian(N,dx, V_x, h_bar, mass)   # Solution for y-direction

    energy = energy_x + energy_y; print(energy.shape)               # Total Eigen Energy values
    plotWavefunc(x, V_x, n_states, wave_vec_x, wave_vec_y, energy, 'Near Infinte Potential Well')
    plotProbfunc(x, V_x, n_states, wave_vec_x, wave_vec_y, energy, 'Near Infinte Potential Well')
    return

def finiteWell(N, dx, h_bar, mass, x, n_states):
    n = np.floor(N/3).astype('int32')
    V_x = np.ones((N)) * 10
    V_x[n:2*n] = 0
    wave_vec_x, energy_x, _ = Hamiltonian(N,dx, V_x, h_bar, mass)   # Solution for x-direction
    wave_vec_y, energy_y, _ = Hamiltonian(N,dx, V_x, h_bar, mass)   # Solution for y-direction

    energy = energy_x + energy_y; print(energy.shape)               # Total Eigen Energy values
    plotWavefunc(x, V_x, n_states, wave_vec_x, wave_vec_y, energy, 'Finte Potential Well')
    plotProbfunc(x, V_x, n_states, wave_vec_x, wave_vec_y, energy, 'Finte Potential Well')
    return

def HarmonicOcilator(N, dx, h_bar, mass, x, n_states):
    V_x = x**2
    wave_vec_x, energy_x, _ = Hamiltonian(N,dx, V_x, h_bar, mass)
    wave_vec_y, energy_y, _ = Hamiltonian(N,dx, V_x, h_bar, mass)

    energy = energy_x + energy_y
    plotWavefunc(x, V_x, n_states, wave_vec_x, wave_vec_y, energy, 'Harmonic Ocillator Potential Function')
    plotProbfunc(x, V_x, n_states, wave_vec_x, wave_vec_y, energy, 'Harmonic Ocillator Potential Function')
        
    return

def StepBarrier(N, dx, h_bar, mass, x, n_states):
    V_x = np.zeros((N))
    V_x[N//2:] += 10 
    wave_vec_x, energy_x, _ = Hamiltonian(N,dx, V_x, h_bar, mass)
    wave_vec_y, energy_y, _ = Hamiltonian(N,dx, V_x, h_bar, mass)

    energy = energy_x + energy_y
    plotWavefunc(x, V_x, n_states, wave_vec_x, wave_vec_y, energy, 'Step Potential Function')
    plotProbfunc(x, V_x, n_states, wave_vec_x, wave_vec_y, energy, 'Step Potential Function')
    
    
    return

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('| python file | x limit (L) | number of points (N)| number of eigen states |')
    else:
        xLimit = int(sys.argv[1])                       #user input for x-co-ordinate limit
        N = int(sys.argv[2])                            #user input for the number of samples
        eigen_states = int(sys.argv[3])                 #user input for the number of eigen states to look at

        print('Enter- 0:Free particle | 1: Near infinite well | 2: Finite well | 3: Harmonic Oscillator ')
        
        try:
            potential_func = int(input('Enter the potential function number (0-7):'))         #user input for type of potential energy
        except:
            print('Enter- 0:Free particle | 1: Near infinite well | 2: Finite well | 3: Harmonic Oscillator')

        x = np.linspace(-xLimit, xLimit, N)
        dx = x[1] - x[0]
        h_bar = 1
        mass = 1

        plt.rcParams['figure.figsize'] = (10,9)
        plt.rc('legend', fontsize=12)
        plt.rcParams['legend.frameon'] = True

        if potential_func == 1:
            infiniteWell(N, dx, h_bar, mass, x, eigen_states)

        elif potential_func == 2:
            finiteWell(N, dx, h_bar, mass, x, eigen_states)

        elif potential_func == 3:
            HarmonicOcilator(N, dx, h_bar, mass, x, eigen_states)
