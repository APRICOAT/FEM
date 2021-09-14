# Import dataclass library
from dataclasses import dataclass

##############
# NODE CLASS #
##############

# Node class initialization
@dataclass
class Node:
    
    ######
    # ID #
    ######
    
    # Node ID
    ID: int
    
    ###############
    # COORDINATES #
    ###############
    
    # Radial coordinate
    x: float    
        
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
    
    ##############################################
    # TANGENTIAL IRRIDIATION INDUCED CHANGE RATE #
    ##############################################
    
    # Tangential irridiation induced change rate ( - )
    et: float = 0.0

    ##################
    # INITIALIZATION #
    ##################
    
    # Initialization
    def __init__ ( self, ID: int, x: float ):
        
        # ID allocation
        self.ID = ID
        
        # Radial coordinate definition
        self.x = x    
        
    #############    
    # SET # D11 #
    #############
        
    # Set D11 parameter
    def setd11( self, d11: float ):
        
        ###############
        # CHECK # D11 #
        ###############
        
        # Check if d11 parameter is zero
        if( self.d11 == 0.0 ):
            
            ############
            # ALLOCATE #
            ############
            
            # Set d11 parameter
            self.d11 += d11
        
        # Check if d11 parameter different than zero
        else:
            
            ##########
            # ADJUST #
            ##########
            
            # Adjust weight
            self.d11 /= 2.0
            
            ############
            # ALLOCATE #
            ############
            
            # Set d11 parameter
            self.d11 += d11 / 2.0
            
    #############    
    # SET # D12 #
    #############
        
    # Set D12 parameter
    def setd12( self, d12: float ):
        
        ###############
        # CHECK # D12 #
        ###############
        
        # Check if d12 parameter is zero
        if( self.d12 == 0.0 ):
            
            ############
            # ALLOCATE #
            ############
            
            # Set d12 parameter
            self.d12 += d12
        
        # Check if d12 parameter different than zero
        else:
            
            ##########
            # ADJUST #
            ##########
            
            # Adjust weight
            self.d12 /= 2.0
            
            ############
            # ALLOCATE #
            ############
            
            # Set d12 parameter
            self.d12 += d12 / 2.0            
          
    ############################    
    # SET # RADIAL IRRIDIATION #
    ############################
        
    # Set radial irridiation
    def seter( self, er: float ):
        
        ##############
        # CHECK # ER #
        ##############
        
        # Check if er parameter is zero
        if( self.er == 0.0 ):
            
            ############
            # ALLOCATE #
            ############
            
            # Set er parameter
            self.er += er
        
        # Check if er parameter different than zero
        else:
            
            ##########
            # ADJUST #
            ##########
            
            # Adjust weight
            self.er /= 2.0
            
            ############
            # ALLOCATE #
            ############
            
            # Set er parameter
            self.er += er / 2.0
            
    #############################    
    # SET # TANGENT IRRIDIATION #
    #############################
        
    # Set tangent irridiation
    def setet( self, et: float ):
        
        ##############
        # CHECK # ET #
        ##############
        
        # Check if et parameter is zero
        if( self.et == 0.0 ):
            
            ############
            # ALLOCATE #
            ############
            
            # Set et parameter
            self.et += et
        
        # Check if et parameter different than zero
        else:
            
            ##########
            # ADJUST #
            ##########
            
            # Adjust weight
            self.et /= 2.0
            
            ############
            # ALLOCATE #
            ############
            
            # Set et parameter
            self.et += et / 2.0
    
    #########
    # PRINT #
    #########
        
    # Print
    def print( self ):
        
        # Print header
        print( '--------------------------------------------------------------------' )
        
        # Print ID
        print( 'Node %d' % self.ID )
        
        # Print radial coordinate
        print( 'Radial coordinate [ um ]      = %.2f' % self.x )        
                
        # Print d1 parameter
        print( 'd1 parameter [ MPa ]          = %.2f' % self.d11 )
        
        # Print d2 parameter
        print( 'd2 parameter [ MPa ]          = %.2f' % self.d12 )
        
        # Print radial irridiation induced change rate
        print( 'Radial irridiation ind. rate  = %.2e' % self.er )
        
        # Print tangential irridiation induced change rate
        print( 'Tangent irridiation ind. rate = %.2e' % self.et )
        
        # Print footer
        print( '--------------------------------------------------------------------' )  
