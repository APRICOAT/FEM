#Import dataclass library
from dataclasses import dataclass
import typing

#Import numpy library
import numpy as np

#Import inverse library
from numpy.linalg import inv

#Import functions
from Material import Material
from FEM.Node import Node

#################################
### Element Class             ###
#################################

#Element class init
@dataclass
class Element:
   ######
    # ID #
    ######
    
    # Element ID
    ID: int
    
    #################
    # NODE 1 STRUCT #
    #################
    
    # Node 1 struct
    node1: Node
    
    #################
    # NODE 2 STRUCT #
    #################
    
    # Node 2 struct
    node2: Node
    
    ##########
    # REGION #
    ##########
    
    # Region name
    region: str
    
    ###################
    # MATERIAL STRUCT #
    ###################
    
    # Material struct
    material: Material
    
    #######################
    # CONSTITUTIVE MATRIX #
    #######################
    
    # Constitutive matrix
    C: np.ndarray
    
    ######################
    # IRRIDIATION MATRIX #
    ######################
    
    # Irridiation matrix
    A: np.ndarray
    
    ############
    # G MATRIX #
    ############
    
    # G matrix
    G: np.ndarray

    # B matrix
    B: np.ndarray

    stress_history: typing.List
    strain_history: typing.List

    ####################
    # STIFFNESS MATRIX #
    ####################
    
    # Stiffness matrix
    Ke: np.ndarray
    
    #########################
    # INTERNAL FORCE VECTOR #
    #########################
    
    # Internal force vector
    Fei: np.ndarray
    
    #########################
    # EXTERNAL FORCE VECTOR #
    #########################
    
    # External force vector
    Fee: np.ndarray
    
    #######################
    # DISPLACEMENT VECTOR #
    #######################
    
    # Displacement vector
    ue: np.ndarray

    ################
    # D PARAMETERS #
    ################
    
    # D parameters
    d11: float = 0.0
    d12: float = 0.0
    
    ##########################################
    # RADIAL IRRIDIATION INDUCED CHANGE RATE #
    ##########################################
    
    # Radial irridiation induced change rate ( - )
    er: float = 0.0
    er_prev: float = 0.0

    
    ##############################################
    # TANGENTIAL IRRIDIATION INDUCED CHANGE RATE #
    ##############################################
    
    # Tangential irridiation induced change rate ( - )
    et: float = 0.0
    et_prev: float = 0.0

    f_prime = (0, 0)
    
    
    ##################
    # INITIALIZATION #
    ##################
    
    # Initialization
    def __init__ ( self, ID: int, node1: Node, node2: Node, region:str, material:Material ):
        
        # ID allocation
        self.ID = ID
        
        # Node 1 struct allocation
        self.node1 = node1
        
        # Node 2 struct allocation
        self.node2 = node2
        
        # Region allocation
        self.region = region
        
        # Material struct allocation
        self.material = material

        #Creates a matriz (2x2) full of zeros
        self.Ke = np.zeros( ( 2 , 2 ) )
        #Creates a matrix (2x1) / colum vector dim-2
        self.Fei = np.zeros( ( 2 , 1 ) )
        #Creates a matrix (2x1) / colum vector dim-2
        self.Fee = np.zeros( ( 2 , 1 ) )
        #Creates a matrix (2x1) / colum vector dim-2
        self.ue = np.zeros( ( 2 , 1 ) )
        #Creates a multidimensional array 4-elments: [0]*(0,0)^T and [1]*(0,0)^T
        self.strain = (np.zeros((2, 1)), np.zeros((2, 1)))
        #Creates a multidimensional array 4-elments: [0]*(0,0)^T and [1]*(0,0)^T
        self.strain_prev = (np.zeros((2, 1)), np.zeros((2, 1)))
        #Creates a multidimensional array 4-elments: [0]*(0,0)^T and [1]*(0,0)^T
        self.stress = (np.zeros((2, 1)), np.zeros((2, 1)))
        #Creates a multidimensional array 4-elments: [0]*(0,0)^T and [1]*(0,0)^T
        self.stress_prev = (np.zeros((2, 1)), np.zeros((2, 1)))
        #creates an array of seflt.stress elements type
        self.stress_history = [self.stress]
        #creates an array of self.strain elements type
        self.strain_history = [self.strain] 

    #############################
    # SET # CONSTITUTIVE MATRIX #
    #############################
        
    # Set consitutive matrix - According to equation (2)

    def setC(self):
        #Initialize the matrix, defined before in the class declaration
        self.C=np.zeros((2,2))
        #Allocating
        self.C[0,0] =  1.0/self.material.E
        self.C[0,0] = -2.0*self.material.v/self.material.E
        self.C[0,0] = -1.0*self.material.v/self.material.E
        self.C[1,1] =  (1.0-self.material.v)/self.material.E

    ############################
    # SET # IRRIDIATION MATRIX #
    ############################
     # Set irridiation matrix - According to equation (4)
    def setA(self):
        # Initialize irridiation matrix
        self.A = np.zeros( ( 2 , 2 ) )
        # Allocate irridiation matrix
        self.A[ 0 , 0 ] =   1.0 * self.material.K
        self.A[ 0 , 1 ] = - 2.0 * self.material.vc * self.material.K 
        self.A[ 1 , 0 ] = - 1.0 * self.material.vc * self.material.K 
        self.A[ 1 , 1 ] = ( 1.0 - self.material.vc ) * self.material.K

    ##################
    # SET # G MATRIX #
    ##################
        
    # Set G matrix - According to equation (8)
    def setG( self, beta:float ):
        
        ###############
        # CALCULATION #
        ###############
        
        # Calculate G matrix
        self.G = inv( self.C + ( beta * self.material.phi * self.A ) )  

    def setB(self, beta: float, del_phi:float):
        self.B = self.G * (self.C - (1-beta)*del_phi*self.A)
    ############################
    # SET # INITIAL CONDITIONS #
    ############################
        
    # Set initial conditions
    def setInitialConditions( self ):
         
        ################
        # D PARAMETERS #
        ################
        
        #EQ 22(f,g)
        # Set D parameters on element
        self.d11 = self.material.E * ( 1.0 - self.material.v ) / ( ( 1.0 + self.material.v ) * ( 1.0 - ( 2.0 * self.material.v ) ) )
        self.d12 = self.material.E * ( 2.0 * self.material.v ) / ( ( 1.0 + self.material.v ) * ( 1.0 - ( 2.0 * self.material.v ) ) )
        
        #This data comes from material properties
        # Set irridiation on element
        self.er = self.material.er
        self.et = self.material.et
        
        ###############
        # SET # NODES #
        ###############
        
        #Now we set d11 and d12 for each node
        # Set D parameters on node 1
        self.node1.setd11( self.d11 )
        self.node1.setd12( self.d12 )
        
        # Set D parameters on node 2
        self.node2.setd11( self.d11 )
        self.node2.setd12( self.d12 )
        
        #As an initial condition, we set er and et as it was declared in the material properties
        # Set irridiation on node 1
        self.node1.seter( self.material.er )
        self.node1.setet( self.material.et )
        
        # Set irridiation on node 2
        self.node2.seter( self.material.er )
        self.node2.setet( self.material.et )

    #This will work only for it eq 0
    #######################
    # SET # ELEMENT MODEL #
    #######################
        
    # Set element parameters
    # Element stiffness matrix - Ke
    # Element internal force vector - Fei
    # Element external force vector - Fee
    # Element model matrix - Ve

    def setElementParameters(self, step:int):
        if step != 0:
            exit
        else:
            #INITIALIZATION 
            #This will be the linear algebraic equation in (12)
            #stiffness matrix initialization
            self.Ke = np.zeros((2,2))
            #internal force vector
            self.Fei = np.zeros((2,1))
            #external force vector
            self.Fee = np.zeros((2,1))
            #displacement vector
            self.ue = np.zeros((2,1))

            #In order to avoid the exact integral 
            #we approximate with gaussi quadrature
            #Wee need to evaluate in the gauss points

            #number of gauss points
            Ngauss = 2

            #Initialization of gauss quadrature points
            pointsGauss  = np.zeros( (Ngauss,1) )
            #gauss quadrature points definition
            pointsGauss[0,0] =  1.0/np.sqrt(3.0)
            pointsGauss[1,0] = -1.0/np.sqrt(3.0)
            #Initialization of gauss quadrature weights
            weightsGauss  = np.zeros( (Ngauss,1) )
            #gauss weight quadrature definition
            weightsGauss[ 0 , 0 ] = 1.0
            weightsGauss[ 1 , 0 ] = 1.0 

            #JACOBIAN

            #Element length definition
            L=self.node2.x - self.node1.x
            #Jacobian
            J=L/2.0

            for i in range(Ngauss):
                eta=pointsGauss[i,0]
                w=weightsGauss[i,0]
                N1eta = (1.0-eta)/2.0
                N2eta = (1.0+eta)/2.0
                DN1eta = -1.0/2.0
                DN2eta =  1.0/2.0
                #getting the functions evaluated at the nodes 
                #here is the main difference tha I've found respecto to the old colde.
                #Instead of using finite difference approximation
                #we use the FEM interpolant for the derivative
                Dd11 = DN1eta*self.node1.d11+DN2eta*self.node2.d11
                Dd12 = DN1eta*self.node1.d12+DN2eta*self.node2.d12
                #d11 on gauss points d11(x(eta)) using the interpolant again...
                d11 = N1eta*self.node1.d11 + N2eta*self.node2.d11
                #radial coordinate is at gaussian point r = x_i*(N_1(eta)) + x_j*(N_2(eta))
                r = N1eta*self.node1.x + N2eta*self.node2.x
                zeta1 = 2.0 + r*Dd11/d11
                zeta2 =-2.0 + r*Dd12/d11

                Der = DN1eta*self.node1.er + DN2eta*self.node2.er
                Det = DN1eta*self.node1.et + DN2eta*self.node2.et

                er = N1eta*self.node1.er + N2eta*self.node2.er
                et = N1eta*self.node1.et + N2eta*self.node2.et

                lambda1 = (Der-Det)+Der*(1.0+self.material.v)/(1.0-self.material.v)
                lambda2 = 2.0*(er-et)*(1.0-2.0*self.material.v)/(1.0-self.material.v)
                lambda2 += r*Dd11*er/d11 + r*Dd12*et/d11

                #With all the elements evaluated at the nodes
                #We can generate the matrices
                #rembeber we are inside of the loop over the gaussian points
                self.Ke[0,0] += 4.0*np.pi*w*(r*r*DN1eta*DN1eta + (2.0-zeta1)*r*N1eta*DN1eta - zeta2*J*N1eta*N1eta)
                self.Ke[0,1] += 4.0*np.pi*w*(r*r*DN1eta*DN2eta + (2.0-zeta1)*r*N1eta*DN2eta - zeta2*J*N1eta*N2eta)
                self.Ke[1,0] += 4.0*np.pi*w*(r*r*DN2eta*DN1eta + (2.0-zeta1)*r*N2eta*DN1eta - zeta2*J*N2eta*N1eta)
                self.Ke[1,1] += 4.0*np.pi*w*(r*r*DN2eta*DN2eta + (2.0-zeta1)*r*N2eta*DN2eta - zeta2*J*N2eta*N2eta)

                rminus = self.node1.x
                rplus  = self.node2.x

                self.Fee[0,0] += -4.0*w*np.pi*r*J*(r*lambda1 + lambda2)*N1eta 
                self.Fee[1,0] += -4.0*w*np.pi*r*J*(r*lambda1 + lambda2)*N2eta 

    def get_zeta3(self):
        return self.G[0,1]/self.G[0,0]       

    # UPdate er and et
    def update(self, phi: float, case:str):
        self.er_prev = self.er
        self.et_prev = self.et
        self.setIrradiationCase(case, phi)

    # COPIED FROM Material.py
    def setIrradiationCase( self, irrCase: str, phi: float ):
        
        # Correlation case allocation
        self.irrCase = irrCase

        # Correlation case B
        if( self.irrCase == 'b' ):
            
            # Radial irridiation induced change calculation
            self.er = + ( 1.36334e-3 * ( phi**3.0 ) ) - ( 7.77024e-3 * ( phi**2.0 ) ) + ( 2.00861e-2 * ( phi**1.0 ) ) - 2.22642e-2
            
            # Tangent irridiation induced change calculation
            self.et = - ( 3.53804e-4 * ( phi**3.0 ) ) + ( 1.67251e-3 * ( phi**2.0 ) ) + ( 2.63307e-3 * ( phi**1.0 ) ) - 1.91253e-2
            
        # Correlation case C
        if ( self.irrCase == 'c' ):
            
            # Radial irridiation induced change calculation
            self.er = + ( 4.03266e-4 * ( phi**3.0 ) ) - ( 2.25937e-3 * ( phi**2.0 ) ) + ( 9.82884e-3 * ( phi**1.0 ) ) - 1.80613e-2
            
            # Tangent irridiation induced change calculation
            self.et = - ( 4.91648e-4 * ( phi**3.0 ) ) + ( 2.32979e-3 * ( phi**2.0 ) ) + ( 1.71315e-3 * ( phi**1.0 ) ) - 1.78392e-2
            
        # Correlation case C
        if ( self.irrCase == 'd' ):
            
            # Check phi condition
            if( phi <= 6.08 ):
                
                # Radial irridiation induced change calculation
                self.er = + ( 4.52013e-4 * ( phi**5.0 ) ) - ( 8.36313e-3 * ( phi**4.0 ) ) + ( 5.67549e-2 * ( phi**3.0 ) ) - ( 1.74247e-1 * ( phi**2.0 ) ) + ( 2.62692e-1 * ( phi**1.0 ) ) - 1.43234e-1
            
                # Tangent irridiation induced change calculation
                self.et = + ( 1.30457e-4 * ( phi**3.0 ) ) - ( 2.10029e-3 * ( phi**2.0 ) ) + ( 9.07826e-3 * ( phi**1.0 ) ) - 3.24737e-2
                
            else:
                
                # Radial irridiation induced change calculation
                self.er = + 0.0954
            
                # Tangent irridiation induced change calculation
                self.et = - 0.0249