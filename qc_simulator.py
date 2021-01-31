import numpy as np
from collections import Counter
from itertools import product


class QuantumCircuit():
    
    ### Class Variables
    
    # Computational Basis
    Zero = np.array([[1.0],[0.0]])
    One = np.array([[0.0],[1.0]])
    
    # Projection operators
    P0 = Zero @  Zero.T
    P1 = One @ One.T

    # Single Qubit gates
    X = np.array([[0,1],[1,0]],dtype = 'complex_')
    Y = np.array([[0,0-1.j],[0+1.j,0]],dtype = 'complex_')
    Z = np.array([[1,0],[0,1]],dtype = 'complex_')
    H = (1/np.sqrt(2))*np.array([[1,1],[1,-1]],dtype = 'complex_')
    I2 = np.array([[1,0],[0,1]],dtype = 'complex_')

    def __init__(self,nqubits:int):
        self.nqubits = nqubits
           
    def get_groundstate(self):
        return self.ntensor_product(*[QuantumCircuit.Zero]*self.nqubits)
    
    def get_std_basis(self):
        """Returns a list of standard basis states for a composite system of #nqubits"""

        # obtain a list of all combinations of (0,1) foe given no.of qubits
        basis_states_list = list(product(['0', '1'], repeat=self.nqubits)) 
        # concat the entries in states_list to get basisvectors
        basis_states = [''.join(element) for element in basis_states_list] 
        return list(map(lambda x : '|'+str(x)+'>',basis_states)) # modify into ket vectors

    
    def get_finalstate(self,program,istate):
        """Given a program (quantum-circuit), execute it on a given initial state"""

        fstate = np.copy(istate) # copy the istate vector for fstate
        for operation in program:
            #print(self.nqubits)
            #print(type(operation['gate']),operation['target'])
            matrix = self.get_operator(operation['gate'],operation['target'])
            fstate = matrix @ fstate  # matrix product
        return fstate

    def get_operator(self,gatename,target_qubits):
        """Return the matrix for a specified Q-gate, no of qubits and target_qubit"""

        single_qgates = {'I2':QuantumCircuit.I2,'X':QuantumCircuit.X,
                         'Y':QuantumCircuit.Y,'Z':QuantumCircuit.Z,
                         'H':QuantumCircuit.H}
        
        
        if gatename in single_qgates:

            # create an array of Identity operations of len = no of qubits
            
            multiqubitgate_list = np.array([single_qgates['I2'] for i in range(self.nqubits)])
            # replace Id by Unitary Gate for the target qubit
            multiqubitgate_list[target_qubits[0]] = single_qgates[gatename]
            result = QuantumCircuit.ntensor_product(*multiqubitgate_list)


        elif gatename == 'CNOT': 
            # for CNOT we need two inputs the (control,target) qubits

            control, target = target_qubits

            # create two arrays of Identity operations for CNOT
            multiqubitgate_list1 = np.array([single_qgates['I2'] for i in range(self.nqubits)])
            multiqubitgate_list2 = np.array([single_qgates['I2'] for i in range(self.nqubits)])

            # replace control qubit in list1 and target qubit in list2
            multiqubitgate_list1[control] = QuantumCircuit.P0
            multiqubitgate_list2[control] = QuantumCircuit.P1
            multiqubitgate_list2[target] = single_qgates['X']

            result = QuantumCircuit.ntensor_product(*multiqubitgate_list1) + QuantumCircuit.ntensor_product(*multiqubitgate_list2) 
            
        return result
        
    def measurement(self,statevector):
        """Perform measurement on given statevector"""
        probability_dist = abs((statevector**2).real)
        return np.random.choice(range(len(statevector)),1, p=probability_dist.flatten())[0]
    

    def get_counts(self,statevector,shots):
        """Get Counts for given number of measurement shots"""

        measurement_results = [self.measurement(statevector) for i in range(shots)] #get measurements for all shots

        counter = Counter(measurement_results)  #creates a dictionary with index of basis state and frequency
        values = counter.values() 
        
        #get the ket vectors for the indices in counter
        newkeys_for_counter = [self.get_std_basis()[int(i)] for i in counter.keys()] 

        return dict(zip(newkeys_for_counter,values))
    
    @classmethod
    def ntensor_product(cls,*args):
        """Calculate the Tensor product over a variable number of inputs"""
        result = np.array([[1.0]])
        for element in args:
            result = np.kron(result, element)
        return result

