import os
import sys
import matplotlib.pyplot as plt

from numpy.lib.function_base import place

sys.path.append(os.path.dirname(__file__))
# Import dataclass library
from dataclasses import dataclass

# Import list library
import typing

# Import numpy library
import numpy as np

# Import functions
from ModelInput import ModelInput
from Material import Material
from Node import Node
from Element import Element
from Layer import Layer

########################
# FINITE ELEMENT CLASS #
########################

# Finite element class initialization
@dataclass
class FEM:
    """
    Main class that include layers
    """
    layers: typing.List[Layer]

    # Method multiplier ( 0 = explicit / 0.5 = Crank Nicholson / 1.0 = implicit )
    beta: float

    S: 'np.ndarray'

    step_size: float = 0.001

    current_step: int = 0

    # Initialization
    def __init__(self, model: ModelInput, layers: typing.List[Layer], steps: int):
        self.layers = layers

        # Method multiplier
        self.beta = model.beta

        # duration
        self.total_duration = model.EFDP

        # model
        self.model = model

        # steps
        self.steps = steps

        self.del_phi = self.model.endLifeFluence / self.steps

        self.phi = 0.0

    # Assembly finite element problem
    def assembly(self):
        # CALCULATE DIMENSIONS of matrix
        num_layers = len(self.layers)
        mat_size = int(num_layers * (self.model.Nelements + 1))
        s_dim = (mat_size, mat_size)

        s = np.zeros(s_dim)

        # assemble layers and add them to S matrix
        for i, layer in enumerate(self.layers):
            layer.assembly()
            shape = layer.K.shape
            r, c = (shape[0] * i, shape[1] * i)
            # Insert layer.K in appropriate position in self.S
            s[r:r + shape[0], c: c + shape[1]] += layer.K

        def M(j):
            return (self.model.Nelements + 1) * (j + 1) - 1

        # Create Z matrix
        nnodes = int(self.model.Nelements + 1)
        Z = np.zeros((num_layers - 1, num_layers * nnodes))
        zr, zc = Z.shape
        for i in range(zr):
            for j in range(zc):
                if j == M(i):
                    Z[i, j] = 1
                elif j == M(i) + 1:
                    Z[i, j] = -1

        def delta_(elem, node):
            return 4 * np.pi * node.x ** 2 / elem.G[0, 0]

        # Create Z_prime matrix
        Z_prime = np.zeros((num_layers * nnodes, num_layers - 1))
        zr, zc = Z_prime.shape

        for i in range(zr):
            for j in range(zc):
                if i == M(j):
                    layer = self.layers[i // nnodes]
                    elem = layer.elements[-1]
                    delta = delta_(elem, elem.node2)
                    Z_prime[i, j] = -delta
                elif i == M(j) + 1 and (i // nnodes) < (num_layers-2):
                    layer = self.layers[i // nnodes + 1]
                    elem = layer.elements[0]
                    delta = delta_(elem, elem.node1)
                    Z_prime[i, j] = delta

        # Create H vector
        h = np.vstack([x.Fe for x in self.layers])
        H = np.zeros((s.shape[0] + Z.shape[0], 1))
        hr, hc = h.shape
        H[0: hr, 0: hc] += h

        H_prime = np.zeros((s.shape[0] + Z.shape[0], 1))
        h_prime = np.zeros((s.shape[0], 1))

        p_internal = (self.model.endLifeInternalPressure - 0) * (self.current_step * self.step_size) / (self.total_duration - 0)

        first_elem = self.layers[0].elements[0]
        last_elem = self.layers[-1].elements[-1]
        h_prime[0, 0] = p_internal * delta_(first_elem, first_elem.node1)
        h_prime[s.shape[0]-1, 0] = self.model.ambientPressure * delta_(last_elem, last_elem.node2)

        #print("A ver: : ", self.layers[0].elements[0])
        H_prime[0:hr, 0:hc] += h_prime

        self.H = H
        self.H_prime = H_prime

        Sr = s.shape[0] + Z.shape[0]
        Sc = s.shape[1] + Z_prime.shape[1]
        self.S = np.zeros((Sr, Sc))
        # Insert s
        self.S[0: s.shape[0], 0: s.shape[1]] += s
        # Insert Z
        self.S[s.shape[0]:Sr, 0: Z.shape[1]] += Z
        # Insert Z_prime
        self.S[0: Z_prime.shape[0], s.shape[1]: Sc] += Z_prime

    # Solve element problem
    def solve( self ):
        ###############
        # FEM PROBLEM #
        ###############
        # Sum the internal and external forces
        self.H_sum = self.H + self.H_prime
        ###################
        # SOLVE # K.u = F #
        ###################
        # Solve FEM problem
        self.u = np.linalg.solve(self.S, self.H_sum)
        #for i in range(len(self.u)-1):
        #    print(self.u[i])
        
        #x = range(len(self.u)-3) 
        #y = self.u[x] 
        #plt.title("Matplotlib demo") 
        #plt.xlabel("x axis caption") 
        #plt.ylabel("y axis caption") 
        #plt.plot(x,y) 
        #plt.show()
        ###################
        # LOOP # ELEMENTS #
        ###################
        # Loop over elements
        return
        for i in range( self.Nelements.astype( int ) ):
            ###############################
            # ALLOCATE RESULTS # ELEMENTS #
            ###############################
            # Allocation results
            self.elements[ i ].ue[ 0 , 0 ] = self.u[ i + 0 , 0 ]
            self.elements[ i ].ue[ 1 , 0 ] = self.u[ i + 1 , 0 ]

    # Calculate stress and strain
    def calculate_stress_strain(self):
        # we have access to self.u
        nnodes = self.model.Nelements + 1

        for layer_index, layer in enumerate(self.layers):
            start = int(layer_index * nnodes)
            end = int(start + nnodes)
            layer_u = self.u[start: end, 0:1]
            layer.calculate_stress_strain(layer_u, debug=layer_index > 0)

    # Update function for updating values of phi, er and et for each step
    def update(self):
        self.current_step += 1

        #Update phi
        self.phi = self.current_step * self.del_phi

        # Update elements er and et, i.e set er_prev and et_prev
        for layer in self.layers:
            layer.update(self.phi, 'b')

    # Solving for multiple steps
    def solve_all_steps(self):
        while True:
            if self.current_step >= self.steps:
                break

            print("step", self.current_step)
            print('------------------')
            self.assembly()
            self.solve()
            self.update()
            self.calculate_stress_strain()