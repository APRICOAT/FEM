# Import dataclass library
from dataclasses import dataclass

# Import pandas library
import pandas as pd

###############
# MODEL CLASS #
###############

# Model class initialization
@dataclass
class ModelInput:
    
    ###################
    # KERNEL DIAMETER # 1
    ###################
    
    # Kernel diameter [ um ]
    kernelDiameter: float
    
    ####################
    # BUFFER THICKNESS # 2
    ####################
    
    # Buffer thickness [ um ]
    bufferThickness: float
        
    ##################
    # IPYC THICKNESS # 3
    ##################
    
    # IPyC thickness [ um ]
    IPyCThickness: float
    
    #################
    # SIC THICKNESS # 4
    #################
    
    # SiC thickness [ um ]
    SiCThickness: float
    
    ##################
    # OPYC THICKNESS #5
    ##################
    
    # OPyC thickness [ um ]
    OPyCThickness: float
    
    ##################
    # KERNEL DENSITY #6
    ##################
    
    # Kernel density [ g / cm³ ]
    kernelDensity: float
    
    ##################
    # BUFFER DENSITY #7
    ##################
    
    # Buffer density [ g / cm³ ]
    bufferDensity: float
        
    ################
    # IPYC DENSITY #8
    ################
    
    # IPyC density [ g / cm³ ]
    IPyCDensity: float
    
    ###############
    # SIC DENSITY #9
    ###############
    
    # SiC density [ g / cm³ ]
    SiCDensity: float
    
    ################
    # OPYC DENSITY #10
    ################
    
    # OPyC density [ g / cm³ ]
    OPyCDensity: float
    
    ############
    # IPYC BAF #11
    ############
    
    # IPyC BAF [ - ]
    IPyCBAF: float
    
    ############
    # OPYC BAF #12
    ############
    
    # OPyC BAF [ - ]
    OPyCBAF: float
    
    ########################
    # IRRIDIATION DURATION #13
    ########################
    
    # Irridiation duration [ - ]
    EFDP: float
    
    #####################
    # END OF LIFE BUMUP #14
    #####################
    
    # End of life bumup [ % ]
    endLifeBumup: float
    
    #######################
    # END OF LIFE FLUENCE #15
    #######################
    
    # End of life fluence [ MeV ]
    endLifeFluence: float
    
    ###########################
    # IRRIDIATION TEMPERATURE # 16
    ###########################
    
    # Irridation temperature [ ºC ]
    irridiationTemperature: float
    
    ########################
    # END OF LIFE PRESSURE # 17
    ########################
    
    # End of life internal pressure [ MPa ]
    endLifeInternalPressure: float
    
    ####################
    # AMBIENT PRESSURE # 18
    ####################
    
    # Ambient pressure [ MPa ]
    ambientPressure: float    
    
    ######################
    # NUMBER OF ELEMENTS # 19
    ######################
    
    # Number of elements
    Nelements: int 
    
    #####################
    # METHOD MULTIPLIER # 20
    #####################
    
    # Method multiplier
    beta: float 
    
    ##################
    # INITIALIZATION # 
    ##################
    
    # Initialization
    def __init__ ( self, input_xlsx_path: str ):
        
        #############
        # READ FILE #
        #############
        
        # Read excel file  ## pandas is used here
        excel = pd.read_excel( input_xlsx_path )
        
        ############
        # ALLOCATE #
        ############
        
        # Kernel diameter [ um ] 1
        self.kernelDiameter = excel.iloc[ 0 ][ 1 ]
        
        # Buffer thickness [ um ] 2
        self.bufferThickness = excel.iloc[ 1 ][ 1 ]
        
        # IPyC thickness [ um ] 3
        self.IPyCThickness = excel.iloc[ 2 ][ 1 ]
        
        # SiC thickness [ um ] 4
        self.SiCThickness = excel.iloc[ 3 ][ 1 ]
        
        # OPyC thickness [ um ] 5
        self.OPyCThickness = excel.iloc[ 4 ][ 1 ]
        
        # Kernel density [ g / cm³ ] 6
        self.kernelDensity = excel.iloc[ 5 ][ 1 ]
        
        # Buffer density [ g / cm³ ] 7
        self.bufferDensity = excel.iloc[ 6 ][ 1 ]
        
        # IPyC density [ g / cm³ ] 8
        self.IPyCDensity = excel.iloc[ 7 ][ 1 ]
        
        # SiC density [ g / cm³ ] 9
        self.SiCDensity = excel.iloc[ 8 ][ 1 ]
        
        # OPyC density [ g / cm³ ] 10
        self.OPyCDensity = excel.iloc[ 9 ][ 1 ]
        
        # IPyC BAF [ - ] 11
        self.IPyCBAF = excel.iloc[ 10 ][ 1 ]
        
        # OPyC BAF [ - ] 12
        self.OPyCBAF = excel.iloc[ 11 ][ 1 ]
        
        # Irridiation duration [ - ] 13
        self.EFDP = excel.iloc[ 12 ][ 1 ]
        
        # End of life bumup [ % ] 14
        self.endLifeBumup = excel.iloc[ 13 ][ 1 ]
        
        # End of life fluence [ MeV ] 15
        self.endLifeFluence = excel.iloc[ 14 ][ 1 ]
        
        # Irridiation temperature [ ºC ] 16
        self.irridiationTemperature = excel.iloc[ 15 ][ 1 ]

        # End of life internal pressure [ MPa ] 17
        self.endLifeInternalPressure = excel.iloc[ 16 ][ 1 ]
        
        # Ambient pressure [ MPa ] 18
        self.ambientPressure = excel.iloc[ 17 ][ 1 ]
        
        # Number of elements 19
        self.Nelements = excel.iloc[ 18 ][ 1 ]
        
        # Method multiplier 20
        self.beta = excel.iloc[ 19 ][ 1 ]
   
    #########
    # PRINT #
    #########
        
    # Print
    def print( self ):
        
        # Print header
        print( '--------------------------------------------------------------------' )
        
        # Print kernel diameter [ um ]
        print( 'Kernel diameter [ um ]        = %.2f' % self.kernelDiameter )
        
        # Print buffer thickness [ um ]
        print( 'Buffer thickness [ um ]       = %.2f' % self.bufferThickness )
        
        # Print IPyC thickness [ um ]
        print( 'IPyC thickness [ um ]         = %.2f' % self.IPyCThickness )
        
        # Print SiC thickness [ um ]
        print( 'SiC thickness [ um ]          = %.2f' % self.SiCThickness )
        
        # Print OPyC thickness [ um ]
        print( 'OPyC thickness [ um ]         = %.2f' % self.OPyCThickness )
        
        # Print kernel density [ g / cm³ ]
        print( 'Kernel density [ g / cm³ ]    = %.2f' % self.kernelDensity )
        
        # Print buffer density [ g / cm³ ]
        print( 'Buffer density [ g / cm³ ]    = %.2f' % self.bufferDensity )
        
        # Print IPyC density [ g / cm³ ]
        print( 'IPyC density [ g / cm³ ]      = %.2f' % self.IPyCDensity )
        
        # Print SiC density [ g / cm³ ]
        print( 'SiC density [ g / cm³ ]       = %.2f' % self.SiCDensity )
        
        # Print OPyC density [ g / cm³ ]
        print( 'OPyC density [ g / cm³ ]      = %.2f' % self.OPyCDensity )
        
        # Print IPyC BAF [ - ]
        print( 'IPyC BAF [ - ]                = %.2f' % self.IPyCBAF )
        
        # Print OPyC BAF [ - ]
        print( 'OPyC BAF [ - ]                = %.2f' % self.OPyCBAF )
        
        # Print irridiation duration [ - ]
        print( 'Irridiation duration [ - ]    = %.2f' % self.EFDP )
        
        # Print end of life bumup [ - ]
        print( 'End of life bumup [ %% ]       = %.2f' % self.endLifeBumup )
        
        # Print end of life fluence [ MeV ]
        print( 'End of life fluence [ MeV ]   = %.2f' % self.endLifeFluence )
        
        # Print irridiation temperature [ ºC ]
        print( 'Irridiation temperature [ºC]  = %.2f' % self.irridiationTemperature )
        
        # Print end of life internal pressure [ MPa ]
        print( 'End of life pressure [ MPa ]  = %.2f' % self.endLifeInternalPressure )
        
        # Print ambient pressure [ MPa ]
        print( 'Ambient pressure [ MPa ]      = %.2f' % self.ambientPressure )
        
        # Print number of elements
        print( 'Number of elements per region = %d' % self.Nelements )
        
        # Print method multipliers
        print( 'Method multiplier - beta      = %.2f' % self.beta )
         
        # Print footer
        print( '--------------------------------------------------------------------' )
