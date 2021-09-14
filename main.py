
import sys
import os
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from FEM.Layer import ModelInput, Material, Node, Element, Layer
from FEM.FEM import FEM

def solve():
    
    pyc_temperature = 1000
    pyc_phi = 1
    pyc_E = 3.96e4
    pyc_v = 0.33
    pyc_p = 1.9e-9
    pyc_a = 5.5e-6
    pyc_S = 200
    pyc_k = 5

    sic_temperature = 1000
    sic_phi = 1
    sic_E = 3.7e5
    sic_v = 0.13
    sic_p = 3.2e-9
    sic_a = 4.9e-6
    sic_S = 873
    sic_k = 8.02


    model = ModelInput.ModelInput('./InputData.xlsx')
    model.print()

    print("@::-->>This is onlye the beginning<<--::@")

    #Once the model parameters have been substracted
    #We cand define the layer lenghts of the CFP.
    l0 = model.kernelDiameter/2.0 #Radius
    l1 = l0 + model.bufferThickness
    l2 = l1 + model.IPyCThickness
    l3 = l2 + model.SiCThickness
    l4 = l3 + model.OPyCThickness
    print("_________________________________")
    print("The boundaries of each layer:")
    print("Kernel radius ::","0",l0)
    print("Buffer domain ::",l0,l1)
    print("IPyC domain ::",l1,l2)
    print("SiC domain ::",l2,l3)
    print("OPyC domain ::",l3,l4)
    print("_________________________________")   

    pyc = Material.Material(1, 'pyc', pyc_temperature, pyc_phi, pyc_E, pyc_v, pyc_p, pyc_a, pyc_S, pyc_k)
    sic = Material.Material(2, 'sic', sic_temperature, sic_phi, sic_E, sic_v, sic_p, sic_a, sic_S, sic_k)



    buffer_layer = Layer(model, pyc, 'Buffer', l0, l1)
    print("-----------D3\n")
    ipyc_layer = Layer(model, pyc, 'IPyC', l1, l2)
    print("-----------D4\n")
    sic_layer = Layer(model, sic, 'SiC', l2, l3)
    print("-----------D5\n")
    opyc_layer = Layer(model, pyc, 'OPyC', l3, l4)
    print("-----------D6\n")

    fem = FEM(model, layers=[buffer_layer, ipyc_layer, sic_layer, opyc_layer], steps=1)
    fem.solve_all_steps()
    return fem



if __name__ == '__main__':
    solve()
