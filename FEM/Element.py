# Import dataclass library
from dataclasses import dataclass
import typing

# Import numpy library
import numpy as np

# Import inverse library
from numpy.linalg import inv

# Import functions
from Material import Material
from FEM.Node import Node

@dataclass
class Element:
    #inputs/definitions
    ID:     int
    node1 : Node
    node2 : Node

    zetaGauss: np.array
    lambdaGauss: np.array


    def __init__(self, ID:int, node1: Node, node2: Node):
        self.ID = ID
        self.node1 = node1
        self.node2 = node2

    def setElementsOne(self):
        L = self.node2.x - self.node1.x
        Ngauss = 2
        #We are going to need the gauss ponits
        eta = np.zeros((2,1))
        eta[0,0] = -1.0/np.sqrt(3.0)
        eta[1,0] = +1.0/np.sqrt(3.0)
        #We need also the weight gauss
        w = np.zeros((2,1))
        w[0,0] = 1.0
        w[1,0] = 1.0

        self.node1.setlambdas_coef()
        self.node1.setzetas()

        self.node2.setlambdas_coef()
        self.node2.setzetas()

        self.zetaGauss = np.zeros((9,2))
        self.lambdaGauss = np.zeros((2,2))

        for k in range(Ngauss):
            N1etak = (1.0 - eta[k,0])/2.0
            N2etak = (1.0 + eta[k,0])/2.0
            DN1etak = (- 1.0)/2.0
            DN2etak = (+ 1.0)/2.0
            #radial coordinate
            rk = 0.5*(1.0-eta[k,0])*self.node1.x + 0.5*(1.0+eta[k,0])*self.node2.x
            rk_g11 = N1etak*(self.node1.x/self.node1.G[0,0]) + N2etak*(self.node2.x/self.node2.G[0,0])
            #We need to compute zetas
            for l in range(9):
                self.zetaGauss[l,k] = (N1etak*(self.node1.zeta[l,0]) + N2etak*(self.node2.zeta[l,0]))



            self.zetaGauss[0,k]  += rk_g11 * (DN1etak*self.node1.G[0,0] + DN2etak*self.node2.G[0,0])
            self.zetaGauss[1,k]  += rk_g11 * (DN1etak*self.node1.G[0,1] + DN2etak*self.node2.G[0,1])
            self.zetaGauss[2,k]  += 0.0
            self.zetaGauss[3,k]  += 0.0 
            self.zetaGauss[4,k]  += 0.0
            self.zetaGauss[5,k]  += rk_g11 * (DN1etak*self.node1.G[0,0] + DN2etak*self.node2.G[0,0])
            self.zetaGauss[6,k]  += rk_g11 * (DN1etak*self.node1.G[0,1] + DN2etak*self.node2.G[0,1])
            self.zetaGauss[7,k]  += rk_g11 * (DN1etak*self.node1.B[0,0] + DN2etak*self.node2.B[0,0])
            self.zetaGauss[8,k]  += rk_g11 * (DN1etak*self.node1.B[0,1] + DN2etak*self.node2.B[0,1])


            self.lambdaGauss[0,k] += (DN1etak*self.node1.lc1 + DN2etak*self.node2.lc1) 
            self.lambdaGauss[0,k] += self.zetaGauss[2,k] * (DN1etak*self.node1.lc2 + DN2etak*self.node2.lc2)
            self.lambdaGauss[0,k] -= self.zetaGauss[3,k] * (DN1etak*self.node1.s_p[0,0] + DN2etak*self.node2.s_p[0,0])
            self.lambdaGauss[0,k] -= self.zetaGauss[4,k] * (DN1etak*self.node1.s_p[1,0] + DN2etak*self.node2.s_p[1,0])
            
            self.lambdaGauss[1,k] += self.zetaGauss[5,k] * (DN1etak*self.node1.lc1 + DN2etak*self.node1.lc1) 
            self.lambdaGauss[1,k] += self.zetaGauss[6,k] * (DN1etak*self.node1.lc2 + DN2etak*self.node2.lc2)
            self.lambdaGauss[1,k] -= self.zetaGauss[7,k] * (N1etak*self.node1.s_p[0,0] + N2etak*self.node2.s_p[0,0])
            self.lambdaGauss[1,k] -= self.zetaGauss[8,k] * (N1etak*self.node1.s_p[1,0] + N2etak*self.node2.s_p[1,0])

        print("zetaGauss: ", self.zetaGauss)
        print("lambdaGauss: ", self.lambdaGauss)

        
        self.Ke = np.zeros((2,2))
        self.Fee = np.zeros((2,1))
        self.Fie = np.zeros((2,1))
        r, c = self.Ke.shape
        for i in range(r):
            for j in range(c):
                for k in range(Ngauss):

                    N1etak = (1.0 - eta[k,0])/2.0
                    N2etak = (1.0 + eta[k,0])/2.0
                    DN1etak = (- 1.0)/2.0
                    DN2etak = (+ 1.0)/2.0
                    #radial coordinate
                    

                    if i == 0:
                        Nikp= DN1etak

                    if i == 1:
                        Nikp= DN2etak


                    if j == 0:
                        Njkp= DN1etak

                    if j == 1:
                        Njkp= DN2etak


                    if i == 0:
                        Nik= N1etak

                    if i == 1:
                        Nik= N2etak


                    if j == 0:
                        Njk= N1etak

                    if j == 1:
                        Njk= N2etak

                    rk = 0.5*(1.0-eta[k,0])*self.node1.x + 0.5*(1.0+eta[k,0])*self.node2.x

                    Jk = L/2.0
                    self.Ke[i,j] += 4.0*np.pi*w[k,0]*(((rk**2)*Nikp*Njkp/Jk) + (2.0-self.zetaGauss[0,k])*rk*Nik*Njkp - self.zetaGauss[1,k]*Jk*Nik*Njk) 
                    #if i == 0 & j==1:
                    #    print("inside:: ",(2.0-self.zetaGauss[0,k])*rk*Nik*Njkp, self.zetaGauss[0,k])
                    #if i == 1 & j==0:
                    #    print("inside:: ",(2.0-self.zetaGauss[0,k])*rk*Nik*Njkp, self.zetaGauss[0,k])

                    self.Fee[i] -= 4.0*np.pi*w[k,0]*rk*Jk*(rk*self.lambdaGauss[0,k] + self.lambdaGauss[1,k])*Nik
