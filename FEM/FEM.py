import os
import sys

from numpy.core.fromnumeric import shape

sys.path.append(os.path.dirname(__file__))
# Import dataclass library
from dataclasses import dataclass

# Import list library
import typing

# Import numpy library
import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
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

    step_size: float = 0.01

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

        self.positions = []

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
            for l in range(len(layer.nodes)):
                self.positions.append(layer.nodes[l].x)
        

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
            return 4 * np.pi * node.x ** 2 / node.G[0, 0]

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

        p_internal = 0
        
        first_elem = self.layers[0].elements[0]
        last_elem = self.layers[-1].elements[-1]
        h_prime[0, 0] = p_internal * delta_(first_elem, first_elem.node1)
        h_prime[s.shape[0]-1, 0] = -self.model.ambientPressure * delta_(last_elem, last_elem.node2)

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
        #print(self.S)
        #print(self.u)
        ru, cu = shape(self.u)
        U = np.zeros((ru,cu))

        for i in range(ru):
            for j in range(cu):
                U[i,j] = self.u[i,j]
        
        #x = self.positions
        #xx = np.zeros((8,1))
        #for i in range(8):
        #    xx[i,0] = x[i]

        #lo = len(x)
        #xo = range(lo)
        #y = self.u[xo]


        #u_r = np.zeros((lo,1))
        #for i in range(lo):
        #    u_r[i,0] = y[i,0]/xx[i,0]

        var =    np.allclose(self.S.dot(self.u), self.H_sum)
        print(self.S.dot(self.u)-self.H_sum)
        print("sol..........")
        print(self.u)
        print("sol..........")
        print("sol:: ",var)

        lN = len(self.layers)
        uo = 0
        for layersN in range(lN):
            les = len(self.layers[layersN].elements)
            for elementsN in range(les):
                print("limits",uo,uo+1)
                print("lay ", layersN,"element ", elementsN) 
                UU = np.zeros((2,1))
                #UU[0:2,0] = U[ layersN*les:(layersN+1)*les,0]
                UU[0:2,0] = U[ uo:uo+1+1,0]
                #print("IOE:.",U[ layersN*les:(layersN+1)*les,0])
                #print("IOE2:.",UU)
                F = self.layers[layersN].elements[elementsN].Fee
                z = np.zeros((2,1))
                for i in range(2):
                    z[i,0] = -F[i,0]
                    for j in range(2):
                        #print("   <<>>",i,j)
                        #print("U>>>>",UU[j,0],"j>>>",j)
                        z[i,0]+=self.layers[layersN].elements[elementsN].Ke[i,j]*UU[j,0]
                    if i == 0:
                        z[i,0] = -z[i,0]/(4.0*np.pi*self.layers[layersN].elements[elementsN].node1.x**2)
                    if i == 1:
                        z[i,0] = +z[i,0]/(4.0*np.pi*self.layers[layersN].elements[elementsN].node2.x**2)
                uo+=1

                print("z:::: ",z)


    def solve_all_steps(self):
        while True:
            if self.current_step >= self.steps:
                break

            print("step", self.current_step)
            print('------------------')
            self.assembly()
            self.solve()
            self.current_step+=1
