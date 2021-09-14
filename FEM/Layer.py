import os
import sys

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

########################
# LAYER CLASS #
########################

# Layer class initialization
@dataclass
class Layer:
    
    ##################
    # INITIAL LENGTH #
    ##################
    
    # Initial length [ um ]
    Li: float
    
    ################
    # FINAL LENGTH #
    ################
    
    # Final length [ um ]
    Lf: float
    
    ################
    # TOTAL LENGTH #
    ################
    
    # Total length [ um ]
    L: float
    
    ###################
    # NUMBER OF NODES #
    ###################
    
    # Number of nodes
    Nnodes: int
    
    ######################
    # NUMBER OF ELEMENTS #
    ######################
    
    # Number of elements
    Nelements: int 
    
    ###############################
    # NUMBER OF DEGREE OF FREEDOM #
    ###############################
    
    # Number of degree of freedom
    Ndofs: int
    
    ################
    # NODES STRUCT #
    ################
    
    # Nodes struct initialization
    nodes: 'typing.Any'
    
    ###################
    # ELEMENTS STRUCT #
    ###################
    
    # Elements struct initialization
    elements: 'typing.Any'
    
    #####################
    # METHOD MULTIPLIER #
    #####################
    
    # Method multiplier ( 0 = explicit / 0.5 = Crank Nicholson / 1.0 = implicit )
    beta: float

    # Layer material type
    material: Material
    
    ###########################
    # GLOBAL STIFFNESS MATRIX #
    ###########################
    
    # Global stiffness matrix
    K: np.ndarray
    
    ################################
    # GLOBAL INTERNAL FORCE VECTOR #
    ################################
    
    # Global internal force vector
    Fi: np.ndarray
    
    ################################
    # GLOBAL EXTERNAL FORCE VECTOR #
    ################################
    
    # Global external force vector
    Fe: np.ndarray
    
    #######################
    # GLOBAL FORCE VECTOR #
    #######################
    
    ##############################
    # GLOBAL DISPLACEMENT VECTOR #
    ##############################
    step_size: float = 0.01

    current_step: int = 0
   
    ##################
    # INITIALIZATION #
    ##################
    
    # Initialization
    def __init__ ( self, model: ModelInput, material: Material, layer_name: str, initial_length: float, final_length: float):
        # Material of the layer
        self.material = material
        self.layer_name = layer_name

        # Initial length ( Kernel radius )
        self.Li = initial_length
        
        # Final length
        #self.Lf = self.Li + model.bufferThickness + model.IPyCThickness + model.SiCThickness + model.OPyCThickness
        self.Lf = final_length
        # Total length
        self.L = self.Lf - self.Li
        
        # Number of elements
        self.Nelements = model.Nelements
        
        # Number of nodes
        self.Nnodes = self.Nelements + 1
        
        # Number of degree of freedom
        self.Ndofs = self.Nnodes
        
        # Method multiplier
        self.beta = model.beta

        # duration
        self.total_duration = model.EFDP

        # model
        self.model = model

        # steps
        self.steps = 20

        self.del_phi = self.model.endLifeFluence / self.steps

        self.phi = 0.0

        ################
        # LOOP # NODES #
        ################
        #Here we are going to generate a 
        # list of nodes
        # node[0], node[1], ... node[Nnodes-1]
        # Nodes structure initialization
        self.nodes = [ ]
        
        # Loop over nodes
        for i in range( self.Nnodes.astype( int ) ):
            
            ##############
            # COORDINATE #
            ##############
            
            # Coordinate definition
            r = self.Li + ( ( self.L / self.Nelements ) * i )
            print("r:::",r, "i::::", i )
            ########
            # NODE #
            ########
            
            # Node struct
            node = Node( i + 1 , r )
            #node.print()
            
            ##############
            # ALLOCATION #
            ##############
            
            # Allocation on nodes structure
            self.nodes.append( node )

        ###################
        # LOOP # ELEMENTS #
        ###################
        #Here we are going to generate a list of Nelements elements
        #elements[0], elements[1], ... , elements[Nelements-1]

        # Elements structure initialization
        self.elements = [ ]
        
        # Loop over elements
        for i in range( self.Nelements.astype( int ) ):
            
            ##########
            # NODE 1 #
            ##########
            
            # Node 1 struct definition
            node1 = self.nodes[ i ]
            
            ##########
            # NODE 2 #
            ##########
            
            # Node 2 struct definition
            node2 = self.nodes[ i + 1 ]
            
            element = Element(i+1, node1, node2, self.layer_name, self.material)
            self.elements.append(element)
            
            ###########################
            # SET CONSTITUTIVE MATRIX #
            ###########################
                    
            # Set constitutive matrix
            self.elements[ i ].setC()
            
            ##########################
            # SET IRRIDIATION MATRIX #
            ##########################
                    
            # Set irridiation matrix
            self.elements[ i ].setA()   
            
            ################
            # SET G MATRIX #
            ################
            
            # Set G matrix
            self.elements[ i ].setG( self.beta )

            # Set B matrix
            self.elements[i].setB(self.beta, self.del_phi)
            
            ##########################
            # SET INITIAL CONDITIONS #
            ##########################
            
            # Set initial conditions
            self.elements[ i ].setInitialConditions()
        #as in the paper Mcapital is the number of degrees of freedom or number of nodes.
        Mcapital = self.Ndofs.astype(int)
        #nelem is the number of elements
        nelem = len(self.elements)
        #I print only for checking 
        print("M and n::",Mcapital,nelem,self.Nelements.astype( int ))

    ############
    # ASSEMBLY #
    ############
    
    # Assembly finite element problem
    def assembly( self ):
        
        ##################
        # INITIALIZATION #
        ##################
        #Mcapital is tthe number of degrees of freedom
        Mcapital = self.Ndofs.astype(int)
        #Mcapital = nelements + 1
        nelements = len(self.elements)

        # Global stiffness initialization
        self.K = np.zeros( ( Mcapital, Mcapital) )
        
        # Global internal forces initialization
        self.Fi = np.zeros( (Mcapital , 1 ) )
        
        # Global external forces initialization
        self.Fe = np.zeros( ( Mcapital , 1 ) )
        
        ###################
        # LOOP # ELEMENTS #
        ###################
        
        # Loop over elements
        for i in range(nelements):
            ##########################
            # SET ELEMENT PARAMETERS #
            ##########################
            
            # Set element parameters
            self.elements[ i ].setElementParameters(self.current_step)
            
            #################################
            # ALLOCATION # GLOBAL STIFFNESS #
            #################################
            
            # Allocation on global stiffness
            self.K[ i + 0 , i + 0 ] += self.elements[ i ].Ke[ 0 , 0 ]
            self.K[ i + 0 , i + 1 ] += self.elements[ i ].Ke[ 0 , 1 ]
            self.K[ i + 1 , i + 0 ] += self.elements[ i ].Ke[ 1 , 0 ]
            self.K[ i + 1 , i + 1 ] += self.elements[ i ].Ke[ 1 , 1 ]
            
            #######################################
            # ALLOCATION # GLOBAL INTERNAL FORCES #
            #######################################
            
            # Allocation global internal forces
            self.Fi[ i + 0 , 0 ] = self.elements[ i ].Fei[ 0 , 0 ]
            self.Fi[ i + 1 , 0 ] = self.elements[ i ].Fei[ 1 , 0 ]
            
            #######################################
            # ALLOCATION # GLOBAL EXTERNAL FORCES #
            #######################################
            
            # Allocation global external forces
            self.Fe[ i + 0 , 0 ] += self.elements[ i ].Fee[ 0 , 0 ]
            self.Fe[ i + 1 , 0 ] += self.elements[ i ].Fee[ 1 , 0 ]

        #Outside the assembly we need to correcto some non-trivial sumations.
        # At the frist element Fe(0,0) we need to add -self.elements[0].f_prime[0]
    
        #last_elem_index = nelements-1 #coz indexing is from 0 to n-1, where n is the dimension
        #print("last element index ", last_elem_index,"Mcapital",Mcapital)
        self.Fe[0,0]        -=  self.elements[0].f_prime[0]
        self.Fe[nelements,0]  +=  self.elements[nelements-1].f_prime[1]

        r1 = self.nodes[0].x
        zeta31 = self.elements[0].get_zeta3()
        self.K[0,0]             -= 4.0*np.pi*r1*zeta31
        rM = self.nodes[nelements].x
        zeta3M = self.elements[nelements-1].get_zeta3()
        self.K[nelements,nelements] += 4.0*np.pi*rM*zeta3M

        Fi_updated = np.zeros((Mcapital, 1))
        p_L = 0.0
        p_R = self.model.ambientPressure
        Fi_updated[0, 0]  =   (4 * np.pi * r1**2 / self.elements[0].G[0, 0]) * p_L
        Fi_updated[nelements, 0] = - (4 * np.pi * rM**2 / self.elements[nelements-1].G[0, 0]) * p_R

    # Calculate stress and strain
    def calculate_stress_strain(self, u: np.ndarray, debug=False):
        def print_(*args, **kwargs):
            if debug:
                print(*args, **kwargs)

        for i in range(self.Nelements.astype(int)):
            elem = self.elements[i]

            if i == 0:
                elem_stress = elem.stress[0]
                print_(self.layer_name, 'Radial Stress', elem_stress[0,0])
                print_(self.layer_name, 'Tangential Stress', elem_stress[1,0])
                print_()

            # Assign previous strain and stress values
            elem.strain_prev = elem.strain
            elem.stress_prev = elem.stress

            # Evaluate du and dr for ratial strain calculation
            #du = u[i+1, 0] - u[i, 0]
            #dr = elem.node2.x - elem.node1.x

            # Radial and tangential strain for first nodes of given element
            strain1 = np.zeros((2, 1))
            strain1[1, 0] = u[i] / elem.node1.x
            #strain1[0, 0] = du / dr
            L = elem.node2.x - elem.node1.x
            dN1dx = - 1.0/L 
            dN2dx = + 1.0/L
            strain1[0,0] = u[i]*dN1dx  +  u[i+1]*dN2dx

            # Radial and tangential strain for first nodes of given element
            strain2 = np.zeros((2, 1))
            strain2[1, 0] = u[i+1] / elem.node2.x
            #strain2[0, 0] = du / dr
            strain2[0, 0] = u[i]*dN1dx  +  u[i+1]*dN2dx
            elem.strain = (strain1, strain2)
            elem.strain_history.append(elem.strain)

            # TODO: MAYBE NEED TO UPDATE THIS
            thermal_strain_delta = np.zeros((2, 1))

            # CALCULATE STRESS

            # B(n) times sigma(n-1) term calculation: The last term in the paper in equation 9
            b_sigma1 = elem.B * elem.stress_prev[0]
            b_sigma2 = elem.B * elem.stress_prev[1]

            # Change in irradiation-induced strain
            er_et_diff_vector = np.zeros((2, 1))
            er_et_diff_vector[0, 0] = elem.er - elem.er_prev
            er_et_diff_vector[1, 0] = elem.et - elem.et_prev

            # Experssion for the middle term in the paper in equation 9
            strain_expression1 = elem.strain_prev[0] + thermal_strain_delta + er_et_diff_vector
            strain_expression2 = elem.strain_prev[1] + thermal_strain_delta + er_et_diff_vector

            # Stress calculation at two nodes of the given element
            stress1 = elem.G * strain1 - elem.G * strain_expression1 + b_sigma1
            stress2 = elem.G * strain2 - elem.G * strain_expression2 + b_sigma2

            elem.stress = (stress1, stress2)
            elem.stress_history.append(elem.stress)

    # Update function for updating values of phi, er and et for each step
    def update(self, phi: float, case: str):
        self.current_step += 1

        # Update elements er and et, i.e set er_prev and et_prev
        for element in self.elements:
           element.update(phi, case)