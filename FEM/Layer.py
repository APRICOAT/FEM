import os
import sys

sys.path.append(os.path.dirname(__file__))
# Import dataclass library
from dataclasses import FrozenInstanceError, dataclass

# Import list library
import typing

# Import numpy library
import numpy as np

# Import functions
from ModelInput import ModelInput
from Material import Material
from Node import Node
from Element import Element

@dataclass
class Layer:
    Li: float
    Lf: float
    Nnodes: int
    elements_list: 'typing.Any'

    def __init__ (self, model: ModelInput, material: Material, layer_name: str, initial_boundary: float, final_boundary:float):
        self.model = model
        self.material = material        #object 
        self.layer_name = layer_name    
        self.Li = initial_boundary      
        self.Lf = final_boundary
        # CH:= computed here
        self.L  = self.Lf - self.Li
        #data extracted from the model, an object
        self.Nelements = model.Nelements
        self.beta       = model.beta

        # CH <<<---
        self.Nnodes = self.Nelements + 1
        self.K = self.material.K
        self.E = self.material.E
        self.nu = self.material.v
        self.mu = self.material.vc
        self.er = self.material.er
        self.et = self.material.et
        self.nodes = []
        self.elements = []

        #loop over the nodes
        #radial step
        dl = self.L/self.Nelements
        for i in range(self.Nnodes.astype(int)):
            if i == 0:
                print("Layer Name:: >>> ", self.layer_name )
            #compute the cordinate
            r = self.Li + dl*i
            print('Coordinate:',r)
            #node
            #numbering as in the paper from 1 to M
            #Here M=number of nodes
            node = Node(i+1,r)
            #append an element to the list
            self.nodes.append(node)
            self.nodes[i].initial_state(self.er, self.et, sr=0, st=0)
            self.nodes[i].initial_params(self.E, self.K, self.mu,self.nu)
            self.nodes[i].setC()
            self.nodes[i].setA()
            self.nodes[i].setX(1,1)
            self.nodes[i].setG(1,1)
            self.nodes[i].setB(1,1)
            self.nodes[i].setzetas()
            self.nodes[i].setlambdas_coef()


            #print(i,self.nodes[i].x,self.nodes[i].ID, self.nodes[i].e[0,0], self.nodes[i].e[1,0])
            #if i==0:
            #    Mtoinv = self.nodes[i].C + self.nodes[i].A
            #    Minv = np.linalg.inv(Mtoinv)
            #    print(i,"G>>",self.nodes[i].G)
            #    print("TestG>>", Minv)
            #    print("C>>", self.nodes[i].C)
            #    print("A>>", self.nodes[i].A)


        #loop over elements
        
        for i in range(self.Nelements.astype(int)):
            #if i == 0:
            #    print("Layer Name:: >>> ", self.layer_name )
            # here i is the index from the list createad above
            node1 = self.nodes[i]
            node2 = self.nodes[i+1]
            element = Element(i+1,node1,node2)
            self.elements.append(element)
            self.elements[i].setElementsOne()
            #print(i,self.elements[i].node1.x,self.elements[i].node2.x,self.elements[i].ID,self.elements[i].node1.e[0,0])
            #if i==3:
            #    Mtoinv = self.elements[i].node1.C + self.elements[i].node1.A
            #    Minv = np.linalg.inv(Mtoinv)
                #print(i,"G>>",self.elements[i].node1.G)
                #print("TestG>>", Minv)
                #print("C>>", self.elements[i].node1.C)
                #print("A>>", self.elements[i].node1.A)

        #assembly

    def assembly(self):
        self.p_L = 0  
        self.p_R = self.model.ambientPressure
        self.K  = np.zeros((self.Nnodes.astype(int),self.Nnodes.astype(int)))
        self.Fi = np.zeros((self.Nnodes.astype(int),1))
        self.Fe = np.zeros((self.Nnodes.astype(int),1))
            # Loop over elements
        for i in range( self.Nelements.astype( int ) ):
                
            # Allocation on global stiffness
            self.K[ i + 0 , i + 0 ] += self.elements[ i ].Ke[ 0 , 0 ]
            self.K[ i + 0 , i + 1 ] += self.elements[ i ].Ke[ 0 , 1 ]
            self.K[ i + 1 , i + 0 ] += self.elements[ i ].Ke[ 1 , 0 ]
            self.K[ i + 1 , i + 1 ] += self.elements[ i ].Ke[ 1 , 1 ]

            # Allocation global internal forces
            self.Fi[ i + 0 , 0 ] = self.elements[ i ].Fie[ 0 , 0 ]
            self.Fi[ i + 1 , 0 ] = self.elements[ i ].Fie[ 1 , 0 ]

            # Allocation global external forces
            self.Fe[ i + 0 , 0 ] += self.elements[ i ].Fee[ 0 , 0 ]
            self.Fe[ i + 1 , 0 ] += self.elements[ i ].Fee[ 1 , 0 ]            

        i_here = 0
        
        zeta3here = self.nodes[i_here].G[0,1]/self.nodes[i_here].G[0,0]
        zeta4here = self.nodes[i_here].B[0,0]/self.nodes[i_here].G[0,0]
        zeta5here = self.nodes[i_here].B[0,1]/self.nodes[i_here].G[0,0]
        fprime0 = self.nodes[i_here].lc1 + zeta3here*self.nodes[i_here].lc2
        fprime0 += -zeta4here*self.nodes[i_here].s_p[0,0] - zeta5here*self.nodes[i_here].s_p[1,0]

        self.K[i_here,i_here] -= 4.0*np.pi*self.nodes[i_here].x * zeta3here
        self.Fe[i_here,0] -= fprime0
        deltaj = 4.0*np.pi*self.nodes[i_here].x**2/self.nodes[i_here].G[0,0]
        self.Fi[i_here,0] += deltaj * self.p_L

        i_here = -1#self.Nelements.astype(int)

        zeta3here = self.nodes[i_here].G[0,1]/self.nodes[i_here].G[0,0]
        zeta4here = self.nodes[i_here].B[0,0]/self.nodes[i_here].G[0,0]
        zeta5here = self.nodes[i_here].B[0,1]/self.nodes[i_here].G[0,0]
        fprimeM = self.nodes[i_here].lc1 + zeta3here*self.nodes[i_here].lc2
        fprimeM += -zeta4here*self.nodes[i_here].s_p[0,0] - zeta5here*self.nodes[i_here].s_p[1,0]
        
        self.K[i_here,i_here] += 4.0*np.pi*self.nodes[i_here].x * zeta3here
        self.Fe[i_here,0] += fprimeM
        deltaj = 4.0*np.pi*self.nodes[i_here].x**2/self.nodes[i_here].G[0,0]
        self.Fi[i_here,0] -= deltaj * self.p_R