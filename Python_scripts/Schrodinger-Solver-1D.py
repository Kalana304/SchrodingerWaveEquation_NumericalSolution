# Author - K.G. Abeywardena
# Date - 29/01/2020

import sys
import numpy as np 
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import matplotlib

plt.style.use('classic')

#Function for creating the Hamiltonian matrix
def Hamiltonian(N, dx, potential, h_bar, mass):    
    """
        This fucntion creates the Hamiltoninan matrix that is required to calculate the numerical solution
        for the Schrodinger's 1-D Wave Equation for a given 1-D potential curve

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
    eigen_energy, eigen_vectors = eigh(H)

    probability = np.abs(eigen_vectors)**2
    K = np.sum(probability, axis=0)*dx
    norm_vectors = eigen_vectors/np.sqrt(K)             #Normalizing the eigen vectors
    probability /= K                                    #Normalizing the PDFs
    return norm_vectors, eigen_energy, probability

#Plotting Functions
def plotWavefunc(x, potential, num_states, wave_vects, energy, func_name):
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
    plt.figure('WaveFunc')
    Labels = []

    for i in range(num_states-1, -1, -1):
        plt.plot(x, wave_vects[:,i] + energy[i])
        Labels.append('E = {:10.4f}eV'.format(energy[i]))

    Labels.append('Potential Energy Function')

    if 'Linear' in func_name:
        plt.plot(x, np.transpose(potential)*max(max(energy[:num_states+1]), max(potential))/max(potential+10**-40), 'k--')
    else:
        plt.plot(x, np.transpose(potential)*min(max(energy[:num_states+1]), max(potential))/max(potential+10**-40), 'k--')
        
    plt.xlabel('Distance (x)', fontdict={'weight': 'bold', 'size': 16, 'color': 'black'})
    plt.ylabel('Energy (eV)', fontdict={'weight': 'bold', 'size': 16, 'color': 'black'})

    plt.legend(Labels, loc='upper right', fontsize = 'medium')
    plt.title(f'Wave Functions for {func_name}', fontdict={'weight': 'bold', 'size': 18, 'color': 'black'})
    plt.grid(b=True)
    plt.autoscale(tight=True, axis='both')
    plt.savefig(f'wavePlot-{func_name}.png')
    plt.autoscale(tight=True)
    plt.close('all')
    return

def plotProbfunc(x, potential, num_states, prob_vects, energy, func_name):
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
    plt.figure('PDF')
    Labels = []

    for i in range(num_states-1, -1, -1):
        plt.plot(x, prob_vects[:,i] + energy[i])
        Labels.append('E = {:10.4f}eV'.format(energy[i]))

    Labels.append('Potential Energy Function')

    if 'Linear' in func_name:
        plt.plot(x, np.transpose(potential)*max(max(energy[:num_states+1]), max(potential))/max(potential+10**-40), 'k--')
    else:
        plt.plot(x, np.transpose(potential)*min(max(energy[:num_states+1]), max(potential))/max(potential+10**-40), 'k--')
            
    plt.xlabel('Distance (x)', fontdict={'weight': 'bold', 'size': 16, 'color': 'black'})
    plt.ylabel('Energy (eV)', fontdict={'weight': 'bold', 'size': 16, 'color': 'black'})
    plt.legend(Labels, loc='upper right', fontsize = 'medium')
    plt.title(f'PDF for {func_name}', fontdict={'weight': 'bold', 'size': 16, 'color': 'black'})
    plt.grid(b=True)
    plt.autoscale(tight=True, axis='both')
    plt.savefig(f'probPlot-{func_name}.png')
    plt.autoscale(tight=True)
    plt.close('all')
    return

#Potential Enegery Functions
def FreeParticle(N, dx, h_bar, mass, x, n_states):
    V_x = np.zeros((N))
    wave_vec, energy, prob = Hamiltonian(N,dx, V_x, h_bar, mass)
    plotWavefunc(x, V_x, n_states, wave_vec, energy, 'Free-Particle')
    plotProbfunc(x, V_x, n_states, prob, energy, 'Free-Particle')
    return

def infiniteWell(N, dx, h_bar, mass, x, n_states):
    n = np.floor(N/3).astype('int32')
    V_x = np.ones((N)) * 2000
    V_x[n:2*n] = 0
    wave_vec, energy, prob = Hamiltonian(N,dx, V_x, h_bar, mass)
    plotWavefunc(x, V_x, n_states, wave_vec, energy, 'Near Infinte Potential Well')
    plotProbfunc(x, V_x, n_states, prob, energy, 'Near Infinte Potential Well')
    return

def finiteWell(N, dx, h_bar, mass, x, n_states):
    n = np.floor(N/3).astype('int32')
    V_x = np.ones((N)) * 10
    V_x[n:2*n] = 0
    wave_vec, energy, prob = Hamiltonian(N,dx, V_x, h_bar, mass)
    plotWavefunc(x, V_x, n_states, wave_vec, energy, 'Finte Potential Well')
    plotProbfunc(x, V_x, n_states, prob, energy, 'Finte Potential Well')
    return

def linearfunc(N, dx, h_bar, mass, x, n_states):
    V_x = x
    wave_vec, energy, prob = Hamiltonian(N,dx, V_x, h_bar, mass)
    plotWavefunc(x, V_x, n_states, wave_vec, energy, 'Linear Potential Function')
    plotProbfunc(x, V_x, n_states, prob, energy, 'Linear Potential Function')
    return

def HarmonicOcilator(N, dx, h_bar, mass, x, n_states):
    V_x = x**2
    wave_vec, energy, prob = Hamiltonian(N,dx, V_x, h_bar, mass)
    plotWavefunc(x, V_x, n_states, wave_vec, energy, 'Harmonic Ocillator Potential Function')
    plotProbfunc(x, V_x, n_states, prob, energy, 'Harmonic Ocillator Potential Function')
    return

def StepBarrier(N, dx, h_bar, mass, x, n_states):
    V_x = np.zeros((N))
    V_x[N//2:] += 10 
    wave_vec, energy, prob = Hamiltonian(N,dx, V_x, h_bar, mass)
    plotWavefunc(x, V_x, n_states, wave_vec, energy, 'Step Potential Function')
    plotProbfunc(x, V_x, n_states, prob, energy, 'Step Potential Function')
    return

def Triangular(N, dx, h_bar, mass, x, n_states):
    n = np.floor(N/3).astype('int32')
    V_x = np.ones((N)) * 100
    V_x[n:2*n] = x[n:2*n]*40
    wave_vec, energy, prob = Hamiltonian(N,dx, V_x, h_bar, mass)
    plotWavefunc(x, V_x, n_states, wave_vec, energy, 'Triangular Potential Function')
    plotProbfunc(x, V_x, n_states, prob, energy, 'Traingular Potential Function')
    return

def CustomBarrier(N, dx, h_bar, mass, x, n_states):
    n = np.floor(N/20).astype('int32')
    V_x = np.zeros((N))
    V_x[n*10:n*11] = 5
    wave_vec, energy, prob = Hamiltonian(N,dx, V_x, h_bar, mass)
    plotWavefunc(x, V_x, n_states, wave_vec, energy, 'Custom Potential Function')
    plotProbfunc(x, V_x, n_states, prob, energy, 'Custom Potential Function')
    return

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('| python file | x limit (L) | number of points (N)| number of eigen states |')
    else:
        xLimit = int(sys.argv[1])                       #user input for x-co-ordinate limit
        N = int(sys.argv[2])                            #user input for the number of samples
        eigen_states = int(sys.argv[3])                 #user input for the number of eigen states to look at

        print('Enter- 0:Free particle | 1: Near infinite well | 2: Finite well | 3: Linear | 4: Harmonic Oscillator | 5: Step Barrier | 6: Triangular | 7: Custom Barrier')
        
        try:
            potential_func = int(input('Enter the potential function number (0-7):'))         #user input for type of potential energy
        except:
            print('Enter- 0:Free particle | 1: Near infinite well | 2: Finite well | 3: Linear | 4: Harmonic Oscillator | 5: Step Barrier | 6: Triangular | 7: Custom Barrier')

        x = np.linspace(-xLimit, xLimit, N)
        dx = x[1] - x[0]
        h_bar = 1
        mass = 1

        plt.rcParams['figure.figsize'] = (10,8)
        plt.rc('legend', fontsize=12)
        plt.rcParams['legend.frameon'] = True

        if potential_func == 0:
            FreeParticle(N, dx, h_bar, mass, x, eigen_states )

        elif potential_func == 1:
            infiniteWell(N, dx, h_bar, mass, x, eigen_states )
        
        elif potential_func == 2:
            finiteWell(N, dx, h_bar, mass, x, eigen_states )

        elif potential_func == 3:
            linearfunc(N, dx, h_bar, mass, x, eigen_states )

        elif potential_func == 4:
            HarmonicOcilator(N, dx, h_bar, mass, x, eigen_states )

        elif potential_func == 5:
            StepBarrier(N, dx, h_bar, mass, x, eigen_states )

        elif potential_func == 6:
            Triangular(N, dx, h_bar, mass, x, eigen_states )

        elif potential_func == 7:
            CustomBarrier(N, dx, h_bar, mass, x, eigen_states)

        else:
            print('Invalida input')

