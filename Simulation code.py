#start of code, import some stuff
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import scipy.sparse.linalg as lin

#input data Rod and clad dimensions from: https://www.researchgate.net/figure/BWR-fuel-assembly-dimensions-Fensin-2004-Mueller-et-al-2013a_tbl1_323487844
#All units are at cm-level
    
Cylinder_Radius = 0.513 #Fuel and Clad are homogenized, this is in cm. 0.513 cm by default
Mesh_Points = int(input("How many mesh points should there be for the First Trial? ")) #How many mesh points are there going to be in our first trial? Changes results due to source term
Source_Strength = int(input("How Strong should the uniform fuel source be? ")) #Strength of uniform fuel source in n/(cm^3*s)
Refinement_Trials=int(input("How many refinements should be done to this mesh (max 7)? "))
while Refinement_Trials>7 or Refinement_Trials<1:
    Refinement_Trials=int(input("Invalid input, how many refinement trials (min 1, max 7)? "))
#Properties for fuel/clad and moderator. Acquired from this source: https://www.oecd-nea.org/science/wprs/eg3drtb/NEA-C5G7MOX.PDF
#Diffusion coefficients are derived from transport cross section (1/(3*(Tr_CS + Abs_CS)))

#Average all of these properties for the 1-group model,
Fuel_Clad_Abs_CS=[8.025e-3,3.717e-3,2.677e-2,9.624e-2,3.002e-2,0.1113,0.2828]
Fuel_Clad_Avg_Abs_CS = sum(Fuel_Clad_Abs_CS)/len(Fuel_Clad_Abs_CS)
Fuel_Clad_Total_CS=[2.12450E-01,3.55470E-01,4.85540E-01,5.59400E-01,3.18030E-01,4.01460E-01,5.70610E-01 ]
Fuel_Clad_Avg_Total_CS = sum(Fuel_Clad_Total_CS)/len(Fuel_Clad_Total_CS)
Fuel_Clad_Fission_CS=[7.212e-3,8.193e-4,6.453e-3,1.8565e-2,1.781e-2,8.303e-2,0.216]
Fuel_Clad_Avg_Fission_CS=sum(Fuel_Clad_Fission_CS)/len(Fuel_Clad_Fission_CS)
Nu_Values=[2.7815,2.4744,2.4338,2.4338,2.4338,2.4338,2.4338] #Neutrons per fission
Avg_Nu=sum(Nu_Values)/len(Nu_Values)
Fuel_Clad_Diff_Coeffs=[1.7924,0.9994,0.6573,0.5123,0.9752,0.6582,0.3935]
Fuel_Clad_Avg_Diff_Coeff=sum(Fuel_Clad_Diff_Coeffs)/len(Fuel_Clad_Diff_Coeffs)

    #No fission in the moderator (duh)
Mod_Diff_Coeffs =[2.0858,0.8071,0.5644,0.5686,0.4606,0.2627,0.1240] 
Mod_Avg_Diff_Coeff = sum(Mod_Diff_Coeffs)/len(Mod_Diff_Coeffs)
Mod_Abs_CS= [6.0105e-4,1.5793e-5,3.37160E-04,1.94060E-03,5.74160E-03,1.50010E-02,3.72390E-02]
Mod_Avg_Abs_CS= sum(Mod_Abs_CS)/len(Mod_Abs_CS)
Mod_Total_CS= [2.30070E-01,7.76460E-01,1.4842,1.5052,1.5592,2.0254,3.3057]
Mod_Avg_Total_CS =  sum(Mod_Total_CS)/len(Mod_Total_CS)

Greatest_Extent=Cylinder_Radius+2*(Mod_Avg_Diff_Coeff)
Step_Width=Greatest_Extent/(Mesh_Points-1)

R_Points=np.linspace(0,Greatest_Extent,Mesh_Points) #Points between centerline and greatest extent
Greatest_Extent=Cylinder_Radius+2*(Mod_Avg_Diff_Coeff)
Flux_Array=np.zeros((len(R_Points)))#Create empty array to store flux values at each point

def Initialize(): #Reset everything for repeated trials
    global Mesh_Points,Greatest_Extent,Step_Width,R_Points,Flux_Array
    Step_Width=Greatest_Extent/(Mesh_Points-1)
    R_Points=np.linspace(0,Greatest_Extent,Mesh_Points) #Points between centerline and greatest extent
    Flux_Array=np.zeros((len(R_Points)))#Create empty array to store flux values at each point

def Get_Properties(mesh_position): #Find properties at a point, maybe between two groups
    global Cur_R,R_minhalf,R_plushalf,Next_R,Prev_R,Flux_plusone,Flux_minone,Cur_Diff_Coeff,Diff_Cur_Flux,D_plushalf,D_minhalf,Fission_Coeff,Removal_Coeff
    #First define everything here as global
    Cur_R = R_Points[mesh_position] #Position of current cell, all that matters for everything but diffusion
    R_minhalf = Cur_R - Step_Width/2 #Position i-1/2
    R_plushalf= Cur_R + Step_Width/2 #Position i+1/2
    Next_R = R_Points[mesh_position] + Step_Width #Position i+1
    Prev_R = R_Points[mesh_position] - Step_Width #Position i-1
    #Finding properties to use, define 

    if Next_R<Cylinder_Radius: #if next R is still in cylinder, it's all in cylinder
        D_plushalf=Fuel_Clad_Avg_Diff_Coeff 
        D_minhalf=Fuel_Clad_Avg_Diff_Coeff
        Cur_Diff_Coeff = Fuel_Clad_Avg_Diff_Coeff
        Fission_Coeff = Fuel_Clad_Avg_Fission_CS*Avg_Nu
        Removal_Coeff= Fuel_Clad_Avg_Abs_CS + Fuel_Clad_Avg_Fission_CS #No in or out scattering in 1-group
 #Because current R is in cylinder
       
    elif Prev_R>Cylinder_Radius: #If Previous R is outside cylinder, all outside
        D_plushalf=Mod_Avg_Diff_Coeff
        D_minhalf=Mod_Avg_Diff_Coeff
        Fission_Coeff = 0 #No fission outside cylinder
        Cur_Diff_Coeff=Mod_Avg_Diff_Coeff
        Removal_Coeff= Mod_Avg_Abs_CS #No in or out scattering in 1-group

    elif Prev_R<Cylinder_Radius: #If previous R is inside cylinder, and Next R isn't, 
        if Cur_R<Cylinder_Radius: #If current R is in cylinder, next R must not be
            D_plushalf=((Mod_Avg_Diff_Coeff**-1 + Fuel_Clad_Avg_Diff_Coeff**-1)/2)**-1  #Harmonic averaging
            D_minhalf=Mod_Avg_Diff_Coeff
            """Make this function use average properties, everything below here should be made into average values"""
            Cur_Diff_Coeff = Fuel_Clad_Avg_Diff_Coeff
            Fission_Coeff = Fuel_Clad_Avg_Fission_CS*Avg_Nu
            Removal_Coeff = Fuel_Clad_Avg_Abs_CS + Fuel_Clad_Avg_Fission_CS

        else: #Current and next R are outside cylinder
            D_plushalf= Fuel_Clad_Avg_Diff_Coeff
            D_minhalf= ((Mod_Avg_Diff_Coeff**-1 + Fuel_Clad_Avg_Diff_Coeff**-1)/2)**-1 
            Fission_Coeff=0
            Cur_Diff_Coeff=Mod_Avg_Diff_Coeff
            Removal_Coeff= Mod_Avg_Abs_CS
    
    if mesh_position==0:
        Diff_Cur_Flux=0 #Removal and fission portion added elsewhere
        Flux_plusone=0 #Zero'd like the diffusion term
        Flux_minone=0
    else:
        Diff_Cur_Flux = -((R_plushalf*D_plushalf + R_minhalf*D_minhalf) / (Cur_R*Step_Width**2))#This one is all negative
        Flux_plusone = R_plushalf*D_plushalf/(Cur_R*Step_Width**2) #Diffusion term coefficient
        Flux_minone=R_minhalf*D_minhalf /(Cur_R*Step_Width**2)


def Numerical_Soln(): #This function was heavily edited by ChatGPT
    Flux_Matr = np.zeros((Mesh_Points, Mesh_Points)) #Create empty matrix for flux values (A matrix)
    Source_Matr = np.zeros(Mesh_Points) #Create empty matrix for source strength (b matrix)
    
    for i in range(Mesh_Points):
        Get_Properties(i) #Acquire material properties for current point

        # -------------------
        # CENTERLINE (symmetry condition)
        # -------------------
        if i == 0:
            
            Flux_Matr[i, i]   = -2 * (Cur_Diff_Coeff / Step_Width**2) + Fission_Coeff - Removal_Coeff #Apply symmetry condition if at r=0
            Flux_Matr[i, i+1] =  2 * Cur_Diff_Coeff / Step_Width**2 #Same as above, but only diffusion portion

            Source_Matr[i] = -Source_Strength #Define source array at this point
            continue

        # -------------------
        # OUTER BOUNDARY (Assumed zero)
        # -------------------
        if i == Mesh_Points - 1:
            Flux_Matr[i, i] = 1.0 #Setting these values ensures flux will always be at zero at the last point (Extrapolation distance)
            Source_Matr[i] = 0.0
            continue

        # -------------------
        # INTERIOR POINTS (Anything between origin and zero bound)
        # -------------------
        #Work is largely done in Get_Properties for these
        Flux_Matr[i, i-1] = Flux_minone
        Flux_Matr[i, i]   = Diff_Cur_Flux + Fission_Coeff - Removal_Coeff
        Flux_Matr[i, i+1] = Flux_plusone

        if Cur_R < Cylinder_Radius: #Only apply source if we're in the fuel cylinder, coolant has no source
            Source_Matr[i] = -Source_Strength
        else:
            Source_Matr[i] = 0.0

    Sparse_A = sparse.csr_matrix(Flux_Matr) #Solve, then return flux matrix
    return lin.spsolve(Sparse_A, Source_Matr)

#Run flux calculation questions
def Find_L2_Error(Prev_Flux,Cur_Flux): #Least mean square error between consecutive trials
    Eo=0 #Difference between current values
    Ep=0 #Difference between next value
    Total_Error=0 #Summation portion of error
    max_error=0 #Greatest point of error

    for z in range(len(Prev_Flux) - 1): #runs up to second to last value in previous array
        #start at 0, iterate to end
        cur_prev_flux = Prev_Flux[z] #Pull current flux value from array
        next_prev_flux = Prev_Flux[z+1] #Pull next flux value from array
        #calculate error, then add to total
        Eo = abs(Cur_Flux[2*z] - cur_prev_flux)
        
        if Eo>max_error:
            max_error=Eo #Find and print out the greatest deviation
        Ep = (Cur_Flux[2*z+2] - next_prev_flux)
        
        Total_Error += Step_Width/2 * (Eo**2 + Ep**2)
    #solve for and print
    L_square_error=np.sqrt(Total_Error) #Get final answer, then print out results
    print("For", len(Prev_Flux), "Points to", len(Cur_Flux),"Points:")
    print("Least squares error = ", L_square_error)
    print("Greatest error =", max_error)

def Plot_Fluxes(): #Print out flux values
    plt.axvline(x=Cylinder_Radius, color='black',linestyle='--', label='Edge of Fuel/Clad') #Indicate edge of the cylinder
    plt.xlabel("Radial distance from center (cm)")
    plt.ylabel("Total Flux (n / (cm\u00b2 * s)")
    plt.title("Radial Distance vs Total Flux at Steady State (S=" + str(Source_Strength) +" n/cm\u00b3*s)")
    plt.legend()
    plt.grid()
    plt.show()
    
#Actually perform convergence test
Color_Array=["Red","Orange","Yellow","Green","Blue","Indigo","Violet"] #Pretty colors that each trial will be printed as
Prev_Array=np.zeros((len(R_Points)))
for N in range(Refinement_Trials):
    Initialize()
    Flux_Array = Numerical_Soln()
    plt.plot(R_Points,Flux_Array,color=str(Color_Array[N]),label="Flux for " + str(Mesh_Points) + " Mesh Points")
    if N!=0: #If first trial, no previous trial to compare to
        Find_L2_Error(Prev_Array,Flux_Array)
    Prev_Array=Flux_Array #Define current array as previous for error comparison
    Mesh_Points=Mesh_Points * 2 #Double number of points for next trial
    
Plot_Fluxes() #Run plotting code to create graph

#End of code

