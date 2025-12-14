#start of code, import some stuff
import numpy as np
import matplotlib.pyplot as plt
import math

#input data Rod and clad dimensions from: https://www.researchgate.net/figure/BWR-fuel-assembly-dimensions-Fensin-2004-Mueller-et-al-2013a_tbl1_323487844
#All units are at cm-level
One_Group_Toggle=True #Toggle whether or not we're using one group or seven, for analytical testing
Test_Convergence=True
if One_Group_Toggle==True: #Set number of groups. Code could theoreticall accomodate any number, but we're only using 1 and 7
    Num_Groups= 1
else:
    Num_Groups = 7
    
Cylinder_Radius = 0.513 #Fuel and Clad are homogenized, this is in cm
Mesh_Points = 50 #How many mesh points are there going to be in our first trial?
Cur_Time=-1 #Define value for current time, which will be used to print graph titles and simplify formulas
Time_Steps=[0.01,0.1,1,2,3,4,5,6,7,8,9,10] #Time steps for finite element transient
Source_Strength=100 #Strength of point centreline source
"""Fuel_Height = 381 #Height of overall rod, doesn't matter due to 1D assumption
#. This is handled at the function for initialization of arrays
#Everything past the cladding is assumed to be coolant, fuel and clad are homogenized
#Properties for fuel/clad and moderator. Acquired from this source: https://www.oecd-nea.org/science/wprs/eg3drtb/NEA-C5G7MOX.PDF"""

""""#Diffusion coefficients are derived from transport cross section (1/(3*(Tr_CS + Abs_CS)))"""
  #Fuel/Clad data
Fuel_Clad_Scat_Matr=[[1.27537E-01, 4.23780E-02, 9.43740E-06, 5.51630E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
[0.00000E+00, 3.24456E-01, 1.63140E-03, 3.14270E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
[0.00000E+00, 0.00000E+00, 4.50940E-01, 2.67920E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 4.52565E-01, 5.56640E-03, 0.00000E+00, 0.00000E+00],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 1.25250E-04, 2.71401E-01, 1.02550E-02, 1.00210E-08],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.29680E-03, 2.65802E-01, 1.68090E-02],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.54580E-03, 2.73080E-01]]

Fuel_Clad_Fission_Spectrum=[0.5879,0.4118,3.391e-4,1.176e-7,0,0,0]
Fuel_Clad_Abs_CS=[8.025e-3,3.717e-3,2.677e-2,9.624e-2,3.002e-2,0.1113,0.2828]
Fuel_Clad_Total_CS=[2.12450E-01,3.55470E-01,4.85540E-01,5.59400E-01,3.18030E-01,4.01460E-01,5.70610E-01 ]
Fuel_Clad_Fission_CS=[7.212e-3,8.193e-4,6.453e-3,1.8565e-2,1.781e-2,8.303e-2,0.216]
Nu_Values=[2.7815,2.4744,2.4338,2.4338,2.4338,2.4338,2.4338] #Neutrons per fission
Fuel_Clad_Diff_Coeffs=[1.7924,0.9994,0.6573,0.5123,0.9752,0.6582,0.3935]

  #Moderator data
Mod_Scat_Matr = [[4.44777E-02, 1.13400E-01, 7.23470E-04, 3.74990E-06, 5.31840E-08, 0.00000E+00, 0.00000E+00],
[0.00000E+00, 2.82334E-01, 1.29940E-01, 6.23400E-04, 4.80020E-05, 7.44860E-06, 1.04550E-06],
[0.00000E+00, 0.00000E+00, 3.45256E-01, 2.24570E-01, 1.69990E-02, 2.64430E-03, 5.03440E-04],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 9.10284E-02, 4.15510E-01, 6.37320E-02, 1.21390E-02],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 7.14370E-05, 1.39138E-01, 5.11820E-01, 6.12290E-02],
[00.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.32440E-01, 2.48070E+00]]
    #No fission in the moderator (duh)
Mod_Diff_Coeffs =[2.0858,0.8071,0.5644,0.5686,0.4606,0.2627,0.1240] 
Mod_Abs_CS= [6.0105e-4,1.5793e-5,3.37160E-04,1.94060E-03,5.74160E-03,1.50010E-02,3.72390E-02]
Mod_Total_CS= [2.30070E-01,7.76460E-01,1.4842,1.5052,1.5592,2.0254,3.3057]

"""#Scattering is from column group to row group 
#(i.e 0.2823 is  group 2 -> group 2, 0.1324 is 7->6 etc.)
#Coordinates in calls are [rows down][columns right] from top left (Which is 0,0)
This is just reference material for me, cuz I tend to forget python does x,y in the opposite way it usually is"""
#def Initialize(Mesh_Points): #Set all of the empty/zeroes arrays, in case multiple trials need to be run. Might be a pain because global variables

Greatest_Extents=np.empty(Num_Groups) #based on zero-edge assumptuion, find the zero-point for each energy group
for i in range(Num_Groups): #Loop populates extrapolation distance array. Contains positions for zero-flux assumption. Constant
    Greatest_Extents[i]=Cylinder_Radius+2*(Mod_Diff_Coeffs[i]) 
    
R_Points=np.linspace(0,max(Greatest_Extents),Mesh_Points) #Points between centerline and greatest extent. Not all groups will reach this far
Step_Width=R_Points[1]-R_Points[0] #Find the width of a step between two slices
Flux_Array=np.zeros((Num_Groups,len(R_Points)))#Create empty array to store flux values at each point
Anal_Flux=np.zeros((len(R_Points))) #Create 1D array to store flux values for 1-group analytical solution
#Loop just defines values for graph testing,

for i in range(Num_Groups): #Run this loop to define placeholder values for Flux_Array
    for j in range(Mesh_Points):
        if R_Points[j]>Greatest_Extents[i]:
            Flux_Array[i][j]=0
            Anal_Flux[j]=0
        else:
            Flux_Array[i][j]=i**2 + j**2
            Anal_Flux[j]=i**2 + j**2 + 1


def Get_Properties(mesh_position,cur_group): #Find properties at a point, maybe between two groups
    global Cur_R,R_minhalf,R_plushalf,Next_R,Prev_R,InScattering_Coeff,Flux_plusone,Flux_minone,Diff_Cur_Flux,D_plushalf,D_minhalf,Fission_Coeff,Inscattering_Coeff,Removal_Coeff
    #First define everything here as global
    Cur_R = R_Points[mesh_position] #Position of current cell, all that matters for everything but diffusion
    R_minhalf = Cur_R - Step_Width/2 #Position i-1/2
    R_plushalf= Cur_R + Step_Width/2 #Position i+1/2
    Next_R = R_Points[mesh_position] + Step_Width #Position i+1
    Prev_R = R_Points[mesh_position] - Step_Width #Position i-1
    #Finding properties to use, define 
    if Next_R<Cylinder_Radius: #if next R is still in cylinder, it's all in cylinder
        D_plushalf=Fuel_Clad_Diff_Coeffs[cur_group] 
        D_minhalf=Fuel_Clad_Diff_Coeffs[cur_group] 
        if One_Group_Toggle==True:
            Fission_Coeff = Fuel_Clad_Fission_CS[0]*Nu_Values[0] 
            InScattering_Coeff=0 #No inscattering with only 1 group
        else:
            Fission_Coeff = -1 #Placeholder, this and inscattering are based on other groups
            InScattering_Coeff= -1 #Placeholder
        Removal_Coeff = Fuel_Clad_Total_CS[cur_group] - Fuel_Clad_Scat_Matr[cur_group][cur_group] #Because current R is in cylinder
        """Fission Coeff and scattering gain are dependant on other energy groups"""
    elif Prev_R>Cylinder_Radius: #If Previous R is outside cylinder, all outside
        D_plushalf=Mod_Diff_Coeffs[cur_group]
        D_minhalf=Mod_Diff_Coeffs[cur_group]
        Fission_Coeff = 0 #No fission outside cylinder
        Removal_Coeff = Mod_Total_CS[cur_group] - Mod_Scat_Matr[cur_group][cur_group] #Removal = total-inscattering
    elif Prev_R<Cylinder_Radius: #If previous R is inside cylinder, and Next R isn't, 
        if Cur_R<Cylinder_Radius: #If current R is in cylinder, next R must not be
            D_plushalf=((Mod_Diff_Coeffs[cur_group]**-1 + Fuel_Clad_Diff_Coeffs[cur_group]**-1)/2)**-1  #Harmonic averaging
            D_minhalf=Mod_Diff_Coeffs[cur_group]
            if One_Group_Toggle==True:
                Fission_Coeff = Fuel_Clad_Fission_CS[0]*Nu_Values[0] 
            else: 
                Fission_Coeff = -1 #Placeholder
            Removal_Coeff = Fuel_Clad_Total_CS[cur_group] - Fuel_Clad_Scat_Matr[cur_group][cur_group] #Current R in cylinder
        else: #Current and next R are outside cylinder
            D_plushalf= Fuel_Clad_Diff_Coeffs[cur_group]
            D_minhalf= ((Mod_Diff_Coeffs[cur_group]**-1 + Fuel_Clad_Diff_Coeffs[cur_group]**-1)/2)**-1 
            Fission_Coeff=0
            Removal_Coeff = Mod_Total_CS[cur_group] - Mod_Scat_Matr[cur_group][cur_group] #Removal = total-inscattering
    Flux_plusone=R_plushalf*D_plushalf/(Cur_R*Step_Width**2) #Diffusion term coefficient
    Diff_Cur_Flux=-((R_plushalf*D_plushalf + R_minhalf*D_minhalf) / (Cur_R*Step_Width**2)) #This one is all negative
    Flux_minone=R_minhalf*D_minhalf /(Cur_R*Step_Width**2)
    
def Analytical_Soln():
    print("yeah")
    
def Numerical_Soln():
    #Spatial discretization (for each group). Currently doing 1-group
    #NOTE: Inscattering 
    for a in range(Num_Groups):
        #Group_v= 1.38e6 * math.sqrt(Group_Avg_Energy[a]) #Averaged velocity of a given group
        for i in range(1, len(R_Points)):#Iterate through mesh points. Start at 1 because 0 results in a singularity
            #DEfining positions    
            Get_Properties(i,a)
            #Diffusion term
            if One_Group_Toggle==True:
                print(Cur_R)
                print("gaming")
                #Now do the numerical solution, add values to empty matrix A. 
                #Source strength is defined at point 0, but is zero everywhere else (Matrix b)
                #Then solve
                
                #Diffusion term, get coefficients that will be used in matrix
                #Fission Term (Already calculated)
                #in-Scattering term (N/A for one group)
                
                #Removal term (Again, already found)
            else:
                print("Multigroup system coming soon")
Numerical_Soln()







def Find_L2_Error(): #For analytical solution, finds Least mean squared error between it and numerical solution
    Eo=0 #Difference between current values
    Ep=0 #Difference between next value
    Total_Error=0 #Summation portion of error
    max_error=0
#L2 code can be found in HW2 Q2C file
    for z in range(len(Flux_Array[0]) - 1): #runs up to second to last value in both analytical and Numerical data arrays
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
    if (One_Group_Toggle == False): #Print out all 7 groups in pretty rainbow colors
        plt.plot(R_Points,Flux_Array[0],color="red",label="Group 1 Flux")
        plt.plot(R_Points,Flux_Array[1],color="orange",label="Group 2 Flux")
        plt.plot(R_Points,Flux_Array[2],color="yellow",label="Group 3 Flux")
        plt.plot(R_Points,Flux_Array[3],color="green",label="Group 4 Flux")
        plt.plot(R_Points,Flux_Array[4],color="blue",label="Group 5 Flux")
        plt.plot(R_Points,Flux_Array[5],color="indigo",label="Group 6 Flux")
        plt.plot(R_Points,Flux_Array[6],color="violet",label="Group 7 Flux")
    else: #If only one group, we only print out the first group and compare it to the analytical solution
        plt.plot(R_Points,Flux_Array[0],color="red",label="Numerical Flux")
        plt.plot(R_Points,Anal_Flux,color="blue",linestyle='-',label='Analytical Flux')
    plt.axvline(x=Cylinder_Radius, color='black',linestyle='--', label='Edge of Fuel/Clad') #Indicate edge of the cylinder
    plt.xlabel("Radial distance from center (cm)")
    plt.ylabel("Group Flux (n / (cm\u00b2 * s)")
    plt.title("Radial Distance vs Flux at t=" + str(Cur_Time) + " Seconds")
    plt.legend()
    plt.show()
    
Plot_Fluxes() #Run plotting code
Find_L2_Error()#Find and print out L2 error if we're doing analytical solution
#End of code

