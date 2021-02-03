import numpy as np
from collections import Counter
from itertools import product
from typing import Optional, Dict


class QuantumCircuit():
    """Implement a simple Quantum Circuit Compiler with 
    Basic Single Qubit Gates and the CNOT gate"""

    # Define Class Variables

    # Computational Basis
    Zero: np.ndarray = np.array([[1.0], [0.0]])
    One: np.ndarray = np.array([[0.0], [1.0]])

    # Projection operators
    P0: np.ndarray = Zero @  Zero.T
    P1: np.ndarray = One @ One.T

    # Single Qubit gates
    X: np.ndarray = np.array([[0, 1], [1, 0]], dtype='complex_')
    Y: np.ndarray = np.array([[0, 0-1.j], [0+1.j, 0]], dtype='complex_')
    Z: np.ndarray = np.array([[1, 0], [0, 1]], dtype='complex_')
    H: np.ndarray = (1/np.sqrt(2))*np.array([[1, 1], [1, -1]], dtype='complex_')
    I2: np.ndarray = np.array([[1, 0], [0, 1]], dtype='complex_')
    S: np.ndarray = np.array([[1, 0], [0, 0+1.j]], dtype='complex_')
    T: np.ndarray = np.array([[1, 0], [0, np.exp(0+1.j*np.pi/4]], dtype='complex_')

    def __init__(self, nqubits: int) -> None:
        if nqubits <= 0:
            raise ValueError(f'Number of qubits must be positive. Given {nqubits}')
        self.nqubits = nqubits


    def get_groundstate(self)-> np.ndarray :
        """Tensor product of #nqubit systems in Zero state"""
        return self.ntensor_product(*[QuantumCircuit.Zero]*self.nqubits)

    def get_std_basis(self)-> np.ndarray :
        """Returns a list of standard basis states for a composite system of #nqubits"""

        # obtain a list of all combinations of (0,1) foe given no.of qubits
        basis_states_list = list(product(['0', '1'], repeat=self.nqubits)) 
        # concat the entries in states_list to get basisvectors
        basis_states = [''.join(element) for element in basis_states_list] 
        return list(map(lambda x : '|'+str(x)+'>',basis_states)) # modify into ket vectors
    
    def get_finalstate(self,
                       program:Dict,
                       initial_state:Optional[np.ndarray] = None)-> np.ndarray:
        """Given a program (quantum-circuit), execute it on a given initial state and return final state"""
        
        if initial_state is None:
            initial_state = self.get_groundstate()
        
        fstate = np.copy(initial_state) # copy the initialstate vector for finalstate
        
        for operation in program:
            matrix = self.get_operator(operation['gate'],operation['target'])
            fstate = matrix @ fstate  # matrix product
        return fstate

    
    def get_operator(self,gatename:str,target_qubits:list) -> np.ndarray:
        """Return the matrix for a specified Q-gate, no of qubits and target_qubit"""

        
        single_qgates = {'I2':QuantumCircuit.I2,'X':QuantumCircuit.X,
                         'Y':QuantumCircuit.Y,'Z':QuantumCircuit.Z,
                         'H':QuantumCircuit.H}

        if target_qubits[0] < 0 or target_qubits[0] > self.nqubits:
            raise IndexError(f'Incorrect entry for target qubit. Please try a value between 0 and {self.nqubits}')
            
        #if gatename not in single_qgates or gatename != 'CNOT':
        #    raise ValueError(f'Incorrect gate name. Expecting only {single_qgates.keys()} gates or CNOT')

        elif gatename in single_qgates:

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

        
    def measurement(self,statevector:np.ndarray) -> int :
        """Perform measurement on given statevector"""
        
        probability_dist = abs((statevector**2).real)
        return np.random.choice(range(len(statevector)),1, p=probability_dist.flatten())[0]
    

    def get_counts(self,statevector:np.ndarray,shots:int) -> Dict[str, int]:
        """Get Counts for given number of measurement shots"""

        measurement_results = [self.measurement(statevector) for i in range(shots)] #get measurements for all shots

        counter = Counter(measurement_results)  #creates a dictionary with index of basis state and frequency
        values = counter.values() 
        
        #get the ket vectors for the indices in counter
        newkeys_for_counter = [self.get_std_basis()[int(i)] for i in counter.keys()] 

        return dict(zip(newkeys_for_counter,values))
    
    @classmethod
    def ntensor_product(cls,*args) -> np.ndarray:
        """Calculate the Tensor product over a variable number of inputs"""
        result = np.array([[1.0]])
        for element in args:
            result = np.kron(result, element)
        return result
    
    @staticmethod
    def normalize_state(vector:np.ndarray) -> np.ndarray:
        """Returns a normalized vector"""
        
        return vectore/np.sum(vector**2)

