# Import dataclass library
from dataclasses import dataclass

# Import list library
import typing

# Import numpy library
import numpy as np

@dataclass
class Node:
    ID:     int #Label of each node
    x :     float # radial coordinate

    #current vector at (n) iteration
    e :     np.ndarray #strain
    e_th :     np.ndarray #strain
    e_sw :     np.ndarray #strain
    s :     np.ndarray #stress

    #vector at (n-1) iteration
    e_p :   np.ndarray #strain at past
    e_p_th :     np.ndarray #strain
    e_p_sw :     np.ndarray #strain
    s_p :   np.ndarray #stress at past

    E:    float 
    K:    float 
    mu:   float 
    nu:   float 

    E_p:    float 
    K_p:    float 
    mu_p:   float 
    nu_p:   float 

    #as a matrix field we need the matrices evaluated at each node
    d_ic:   np.ndarray #matrix reduced at initial conditions
    X :     float
    A :     np.ndarray #matrix at any time
    B :     np.ndarray #matrix at any time
    C :     np.ndarray #matrix at any time
    G :     np.ndarray #matrix at any time

    lam: np.ndarray
    zet: np.ndarray

    ue      = np.ndarray
    ue_p    = np.ndarray
    f_prime = np.ndarray


    
    #INITIALIZE ONLY WITH THE ID AND THE POSITION
    def __init__ (self, ID:int, x: float):
        self.ID = ID
        self.x = x


    def initial_state(self, er:float, et:float, sr:float, st:float):
        #we are goint to get e and s at time (n)
        self.e = np.zeros((2,1))
        self.e_th = np.zeros((2,1))
        self.e_sw = np.zeros((2,1))
        self.s = np.zeros((2,1))
        #as initial condition we are goint to take t=0 as (n-1) iteration
        self.e_p = np.zeros((2,1))
        self.e_p_th = np.zeros((2,1))
        self.e_p_sw = np.zeros((2,1))
        self.s_p = np.zeros((2,1))

        self.e[0,0] = er #er and et are the initial data
        self.e[1,0] = et

        self.s[0,0] = sr #sr and st are the initial data
        self.s[1,0] = st


    def initial_params(self, E:float, K:float, mu: float, nu:float):
        self.E_p = E
        self.K_p = K
        self.mu_p = mu
        self.nu_p = nu

        self.E = E
        self.K = K
        self.mu = mu
        self.nu = nu

    def upgrade_node_state(self):
        self.e_p = self.e
        self.s_p = self.s

    def upgrade_params_state(self):
        self.E_p = self.E
        self.K_p = self.K
        self.mu_p = self.mu
        self.nu_p = self.nu


    #Matrix Eq (2)
    def setC(self):
        self.C = np.zeros((2,2))

        self.C[0,0] = 1.0/self.E
        self.C[0,1] = -2.0*self.nu/self.E
        self.C[1,0] = -self.nu/self.E
        self.C[1,1] = (1.0 - self.nu)/self.E

    #Matrix Eq (4)
    def setA(self):

        self.A = np.zeros((2,2))

        self.A[0,0] = self.K 
        self.A[0,1] = self.K*(-2.0*self.mu) 
        self.A[1,0] = self.K*(-self.mu) 
        self.A[1,1] = self.K*(1.0-self.mu) 


    def setX(self, beta:float, del_phi:float):
        self.X = ((1.0/self.E) + beta * del_phi * self.K )
        self.X *= (((1.0-self.nu)/self.E) + beta * del_phi * self.K *(1.0 - self.mu))
        self.X -= 2.0*((self.nu/self.E) +  beta * del_phi * self.K * self.mu)**2

    def setG(self, beta: float, del_phi: float):

        self.G = np.zeros((2,2))

        self.G[0,0] = (((1.0 - self.nu)/self.E) + beta*del_phi*self.K*(1.0-self.mu))/self.X
        self.G[0,1] = 2.0*((self.nu/self.E) + beta*del_phi*self.K*self.mu)/self.X
        self.G[1,0] = ((self.nu/self.E) + beta*del_phi*self.K*self.mu)/self.X
        self.G[1,1] = ((1.0/self.E) + beta*del_phi*self.K)/self.X


    def setB(self, beta:float, del_phi: float):

        g11 = self.G[0,0]
        g12 = self.G[0,1]
        g21 = self.G[1,0]
        g22 = self.G[1,1]

        self.B = np.zeros((2,2))

        c1 = (1.0/self.E_p) - (1.0-beta) * del_phi*self.K_p 
        c2 = - (self.nu_p/self.E_p) + (1.0-beta) * del_phi * self.K_p * self.mu_p
        self.B[0,0] = g11*(c1)+g12*(c2)
        c3 = c2 
        c4 = ((1.0-self.nu_p)/self.E_p) - (1.0-beta) * del_phi*self.K_p*(1.0-self.nu_p)
        self.B[0,1] = 2.0*g11*(c3) +g12*(c4)
        c5 = c1
        c6 = c2
        self.B[1,0] = g21*(c5) + g22*(c6)
        c7 = c2
        c8 = c4
        self.B[1,1] = 2.0*g21*(c7) + g22*(c8)

    def setlambdas_coef(self):
        deltaer_th = self.e_th[0,0] - self.e_p_th[0,0]
        deltaer_sw = self.e_sw[0,0] - self.e_p_sw[0,0]
        deltaet_th = self.e_th[1,0] - self.e_p_th[1,0]
        deltaet_sw = self.e_sw[1,0] - self.e_p_sw[1,0]
        self.lc1 = self.e_p[0,0] + deltaer_th + deltaer_sw    
        self.lc2 = self.e_p[1,0] + deltaet_th + deltaet_sw




    def setzetas(self):

        g11, g12, g21, g22 = self.G[0,0], self.G[0,1], self.G[1,0], self.G[1,1]
        b11, b12, b21, b22 = self.B[0,0], self.B[0,1], self.G[1,0], self.G[1,1]

        self.zeta = np.zeros((9,1))

        self.zeta[0,0] = (g12 + 2.0 *g11 - 2.0*g21)/g11
        self.zeta[1,0] = (g12 - 2.0 *g22)/ g11
        self.zeta[2,0] = g12/g11
        self.zeta[3,0] = b11/g11
        self.zeta[4,0] = b12/g11
        self.zeta[5,0] = 2.0*((g11-g21)/g11)
        self.zeta[6,0] = 2.0*((g12-g22)/g11)
        self.zeta[7,0] = 2.0*((b11-b21)/g11)
        self.zeta[8,0] = 2.0*((b12-b22)/g11)

        



