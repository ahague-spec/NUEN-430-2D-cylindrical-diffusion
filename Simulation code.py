"""Fuel_Height = 381 #Height of overall rod, doesn't matter due to 1D assumption
#. This is handled at the function for initialization of arrays
#Everything past the cladding is assumed to be coolant, fuel and clad are homogenized
#Properties for fuel/clad and moderator. Acquired from this source: https://www.oecd-nea.org/science/wprs/eg3drtb/NEA-C5G7MOX.PDF"""

""""#Diffusion coefficients are derived from transport cross section (1/(3*(Tr_CS + Abs_CS)))"""
#start of code, import some stuff

"""Current problems:
    Graph changes with mesh size, but shape does not. Source term is most likely to blame
    Need Analytical solutions
    Likely gonna make this a 1-group solution with total flux, averaging cross sections
    We likely won't have time for the multigroup
    Average out material properties
    """

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import sparse
import scipy.sparse.linalg as lin
from scipy.special import i0

#input data Rod and clad dimensions from: https://www.researchgate.net/figure/BWR-fuel-assembly-dimensions-Fensin-2004-Mueller-et-al-2013a_tbl1_323487844
#All units are at cm-level

Test_Convergence=True
    
Cylinder_Radius = 0.513 #Fuel and Clad are homogenized, this is in cm
Mesh_Points = 50 #How many mesh points are there going to be in our first trial? Changes results due to source term
Source_Strength=100 #Strength of uniform fuel source
Total_Source = Source_Strength * np.pi * Cylinder_Radius**2 #Spread source throughout the cylinder. AI generated

  #Fuel/Clad data

#Average all of these properties for the 1-group model, now port them to the rest of the code
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

R_Points=np.linspace(0,Greatest_Extent,Mesh_Points) #Points between centerline and greatest extent. Not all groups will reach this far
#Singularity at zero requires special treatment

Flux_Array=np.zeros((len(R_Points)))#Create empty array to store flux values at each point
Anal_Flux=np.zeros((len(R_Points))) #Create 1D array to store flux values for 1-group analytical solution

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
    
    
def Analytical_Soln():
    print("yeah")
    Anal_Fluxes=np.empty(Mesh_Points)
    #Finish this once it gets fixed
    for i in range(Mesh_Points):
        Get_Properties(i) #Get properties for current point. Always in first group
        B_Term = (Removal_Coeff-Fission_Coeff)/Cur_Diff_Coeff #This is the B2 term
        Extrap_Bessel=i0(math.sqrt(B_Term)*Greatest_Extent) #I0 function for Rex
        Cur_Bessel=i0(math.sqrt(B_Term)*R_Points[i]) #I0 function for r
        if R_Points[i] > Cylinder_Radius: #If we're outside the cylinder
            Part_Soln = Source_Strength/(Removal_Coeff-Fission_Coeff)
            """Fission needs to be taken into account here"""
        else: #We must be in the fuel/clad
            print("In fuel, particular solution is zero")
            Part_Soln=0
    return Anal_Fluxes

def Numerical_Soln(): #This function was heavily edited by ChatGPT
    Group_Matr = np.zeros((Mesh_Points, Mesh_Points))
    Source_Matr = np.zeros(Mesh_Points)
    
    for i in range(Mesh_Points):
        Get_Properties(i)
        r = R_Points[i]

        # -------------------
        # CENTERLINE (symmetry)
        # -------------------
        if i == 0:
            
            Group_Matr[i, i]   = -2 * Cur_Diff_Coeff / Step_Width**2 \
                                 + Fission_Coeff - Removal_Coeff
            Group_Matr[i, i+1] =  2 * Cur_Diff_Coeff / Step_Width**2

            Source_Matr[i] = -Source_Strength
            continue

        # -------------------
        # OUTER BOUNDARY (vacuum)
        # -------------------
        if i == Mesh_Points - 1:
            Group_Matr[i, i] = 1.0
            Source_Matr[i] = 0.0
            continue

        # -------------------
        # INTERIOR POINTS
        # -------------------
        

        Group_Matr[i, i-1] = Flux_minone
        Group_Matr[i, i]   = Diff_Cur_Flux + Fission_Coeff - Removal_Coeff
        Group_Matr[i, i+1] = Flux_plusone

        if r < Cylinder_Radius:
            Source_Matr[i] = -Source_Strength
        else:
            Source_Matr[i] = 0.0

    Sparse_A = sparse.csr_matrix(Group_Matr)
    return lin.spsolve(Sparse_A, Source_Matr)


Flux_Array = Numerical_Soln()
Anal_Flux_Array=Analytical_Soln()

def Find_L2_Error(): #For analytical solution, finds Least mean squared error between it and numerical solution
    Eo=0 #Difference between current values
    Ep=0 #Difference between next value
    Total_Error=0 #Summation portion of error
    max_error=0
#L2 code can be found in HW2 Q2C file
    for z in range(len(Flux_Array) - 1): #runs up to second to last value in both analytical and Numerical data arrays
        #start at 0, iterate to end
        cur_anal_flux = Anal_Flux[z] #Pull current flux value from array
        next_anal_flux = Anal_Flux[z+1] #Pull next flux value from array
        #calculate error, then add to total
        Eo = abs(Flux_Array[0][z] - cur_anal_flux)
        if Eo>max_error:
            max_error=Eo #Find and print out the greatest deviation
        Ep = (Flux_Array[0][z+1] - next_anal_flux)
        Total_Error += Step_Width/2 * (Eo**2 + Ep**2)
    #solve for and print
    L_square_error=math.sqrt(Total_Error) #Get final answer, then print out results
    print("Least squares error = ", L_square_error)
    print("Greatest error =", max_error)

def Plot_Fluxes(): #Print out flux values for each group
    for i in range(Mesh_Points):
        print("Flux at point", R_Points[i], "cm =", Flux_Array[i])
    plt.plot(R_Points,Flux_Array,color="red",label="Numerical Flux")
    #plt.plot(R_Points,Anal_Flux,color="blue",linestyle='-',label='Analytical Flux')
    plt.axvline(x=Cylinder_Radius, color='black',linestyle='--', label='Edge of Fuel/Clad') #Indicate edge of the cylinder
    plt.xlabel("Radial distance from center (cm)")
    plt.ylabel("Group Flux (n / (cm\u00b2 * s)")
    plt.title("Radial Distance vs Total Flux at Steady State")
    plt.legend()
    plt.show()
    
Plot_Fluxes() #Run plotting code
#Find_L2_Error()#Find and print out L2 error if we're doing analytical solution. Finish that
#End of code
