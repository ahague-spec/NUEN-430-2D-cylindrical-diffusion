#start of code, import some stuff
import numpy as np

#input data that is meant to vary
Pitch_Width = float(input("What is the lattice pitch (in cm)? ")) #1.5 is at 9.375 cm
Pu239_Frac = float(input("What is the Pu-239 fraction (in percent)? ")) / 100
Element_Length = float(input("How long is the fuel element (in cm)? "))
Core_Diameter=float(input("How wide is the core (in cm)? "))

Bg2 = (2.405/(20 + Core_Diameter/2))**2 + (3.1415926/(Element_Length+40))**2  #calculate geometric buckling, same as hex reactor because full system is cylinder

#fuel is cylindrical and its diameter is fixed. 30.288=pi*(6.21/2)^2, to simplify calculation. This is in cm^3 by the way

Fuel_Volume = 30.2882 * Element_Length #30.6796 = pi * (6.25/2)^2
Cladding_Volume = 30.6796 * Element_Length - Fuel_Volume 
Cell_Volume= 3.1415926*(Pitch_Width/2)**2 * Element_Length #pi*r^2 * h
#Cell_Volume= Pitch_Width**2 * Element_Length #This one's for square

#find coolant volume
Coolant_Volume= float(Cell_Volume - Cladding_Volume - Fuel_Volume) 
#volume factors for homogenization
Fuel_Frac=float(Fuel_Volume/Cell_Volume)
Clad_Frac=float(Cladding_Volume/Cell_Volume)
Cool_Frac=float(Coolant_Volume/Cell_Volume)
Core_Volume=3.1415926*(Core_Diameter/2)**2 * Element_Length #Volume of whole core


#calculate number of cells, cast to integer to round down
Num_Cells=int(Core_Volume/Cell_Volume) #Do I need to compensate for voids here? Guessing no.
Fuel_Load=0.019*Fuel_Volume*Num_Cells
#Find load of plutonium from enrichment
Pu239_Load = Pu239_Frac * Fuel_Load

#calculate fuel load in kg, using the density of 19.0 g/cm^3
Fuel_Load=0.019*Fuel_Volume*Num_Cells
#Find the load of plutonium from enrichment
Pu239_Load = Pu239_Frac * Fuel_Load
#Volume fractions, used to homogenize properties later
Fuel_Frac=float(Fuel_Volume/Cell_Volume)
Clad_Frac=float(Cladding_Volume/Cell_Volume)
Cool_Frac=float(Coolant_Volume/Cell_Volume)

Flux_Dis_Toggle="n"
#Makes flux disadv factors always equal 1 if toggle isn't enabled
if Flux_Dis_Toggle == "y":
    Fuel_Flux_Dis=np.array([1.60,1.60,1.60,1.45,1.45,1.34,1.34,1.34])
    Cool_Flux_Dis=np.array([0.9,0.9,0.91,0.91,0.91,0.94,0.94,0.94])
    Clad_Flux_Dis=np.array([1.01,1.01,1.01,1.02,1.02,1.03,1.03,1.03])
    
else: 
    #set them all to 1s, effectively ignoring them
    Fuel_Flux_Dis=np.array([1,1,1,1,1,1,1,1])
    Cool_Flux_Dis=np.array([1,1,1,1,1,1,1,1])
    Clad_Flux_Dis=np.array([1,1,1,1,1,1,1,1])

#define constants

Chi_Matr=np.array([0.365,0.396,0.173,0.05,0.012,0.003,0.001,0]) #groups 1-8 fission cross sections for all fissions
Chi1=0.365
Chi2=0.396
Chi3=0.173
Chi4=0.05
Chi5=0.012
Chi6=0.003
Chi7=0.001
Chi8=0
#define and solve for lead data
"""D=1/(3*macro transp. C.S)"""
#removal is fission + capture + scat removal
#done this way instead of in a matrix because python whines at me about depreciation when I pull values from one
#Lead macro comp. factor = 3.2987e22 (barns to cm^-1)
#calculates transport cross section, removal cross section, and diffusion coefficient for each energy group
Lead_Trans_CS1=1.5e-24 * 3.2987e22
Lead_Removal_CS1=0.628e-24 * 3.2987e22 #add radiative capture to removal CONFIRMED RIGHT
Lead_D1= 1/(3*Lead_Trans_CS1)
Lead_Trans_CS2=2.2e-24 * 3.2987e22
Lead_Removal_CS2=0.691e-24 * 3.2987e22
Lead_D2= 1/(3*Lead_Trans_CS2)
Lead_Trans_CS3=3.6e-24 * 3.2987e22
Lead_Removal_CS3=0.4462e-24 * 3.2987e22
Lead_D3= 1/(3*Lead_Trans_CS3)
Lead_Trans_CS4=3.5e-24 * 3.2987e22
Lead_Removal_CS4=0.291e-24 * 3.2987e22
Lead_D4= 1/(3*Lead_Trans_CS4)
Lead_Trans_CS5=4e-24 * 3.2987e22
Lead_Removal_CS5=0.351e-24 * 3.2987e22
Lead_D5= 1/(3*Lead_Trans_CS5)
Lead_Trans_CS6=3.9e-24 * 3.2987e22
Lead_Removal_CS6=0.301e-24 * 3.2987e22
Lead_D6= 1/(3*Lead_Trans_CS6)
Lead_Trans_CS7=7.3e-24 * 3.2987e22
Lead_Removal_CS7=0.049e-24 * 3.2987e22
Lead_D7= 1/(3*Lead_Trans_CS7) 
Lead_Trans_CS8=3.2e-24 * 3.2987e22
Lead_Removal_CS8=0.008e-24 * 3.2987e22
Lead_D8= 1/(3*Lead_Trans_CS8)
#create LGxG matrix for lead  (its fission matrix is all zeroes)
Lead_Scat_Matr=np.array([[Lead_D1*Bg2+Lead_Removal_CS1,0,0,0,0,0,0,0],
                [-0.01715,Lead_D2*Bg2+Lead_Removal_CS2,0,0,0,0,0,0],
                [-0.0029688,-0.02276,Lead_D3*Bg2+Lead_Removal_CS3,0,0,0,0,0],
                [-0.000098961,0,-0.01451428,Lead_D4*Bg2+Lead_Removal_CS4,0,0,0,0],
                [-0.000296883,-0.0000131948,-0.000164935,-0.00956623,Lead_D5*Bg2+Lead_Removal_CS5,0,0,0],
                [-0.000032987,-0.0000131948,-0.0000263896,0,-0.01154545,Lead_D6*Bg2+Lead_Removal_CS6,0,0],
                [0,0,0,0,0,-0.00230909,Lead_D7*Bg2+Lead_Removal_CS7,0],
                [0,0,0,0,0,0,-0.0000131948,Lead_D8*Bg2+Lead_Removal_CS8]])

#Solve for Iron's properties and matrix (cladding material)
#micro to macro factor is 8.475e22 atoms/cm^3
Iron_Trans_CS1= 2.2e-24 * 8.475e22
Iron_Removal_CS1= 1.0308e-24 * 8.475e22
Iron_D1= 1/(3*Iron_Trans_CS1)
Iron_Trans_CS2= 2.1e-24* 8.475e22#do removals
Iron_Removal_CS2= 0.463e-24 * 8.475e22
Iron_D2= 1/(3*Iron_Trans_CS2)
Iron_Trans_CS3= 2.4e-24 * 8.475e22
Iron_Removal_CS3= 0.125e-24 * 8.475e22
Iron_D3= 1/(3*Iron_Trans_CS3)
Iron_Trans_CS4= 3.1e-24 * 8.475e22
Iron_Removal_CS4= 0.146e-24 * 8.475e22
Iron_D4= 1/(3*Iron_Trans_CS4)
Iron_Trans_CS5= 4.5e-24 * 8.475e22
Iron_Removal_CS5= 0.288e-24 * 8.475e22
Iron_D5= 1/(3*Iron_Trans_CS5)
Iron_Trans_CS6= 6.1e-24 * 8.475e22
Iron_Removal_CS6= 0.082e-24 * 8.475e22
Iron_D6= 1/(3*Iron_Trans_CS6)
Iron_Trans_CS7= 6.9e-24 * 8.475e22
Iron_Removal_CS7= 0.072e-24 * 8.475e22
Iron_D7= 1/(3*Iron_Trans_CS7)
Iron_Trans_CS8= 10.4e-24 * 8.475e22
Iron_Removal_CS8= 0.02e-24 * 8.475e22
Iron_D8= 1/(3*Iron_Trans_CS8)
#now construct the matrix for Iron
Iron_Scat_Matr=np.array([[Iron_D1*Bg2+Iron_Removal_CS1,0,0,0,0,0,0,0],
                [-0.0635625,Iron_D2*Bg2+Iron_Removal_CS2,0,0,0,0,0,0],
                [-0.01695,-0.0279675,Iron_D3*Bg2+Iron_Removal_CS3,0,0,0,0,0],
                [-0.042375,-0.008475,-0.01017,Iron_D4*Bg2+Iron_Removal_CS4,0,0,0,0],
                [-0.0008475,-0.001695,0,-0.011865,Iron_D5*Bg2+Iron_Removal_CS5,0,0,0],
                [-0.0000678,-0.0008475,0,0,-0.02373,Iron_D6*Bg2+Iron_Removal_CS6,0,0],
                [0,0,0,0,0,-0.0059325,Iron_D7*Bg2+Iron_Removal_CS7,0],
                [0,0,0,0,0,0,-0.00339,Iron_D8*Bg2+Iron_Removal_CS8]])

#Fuel matrix, sums U and Pu cross section
#U comp factor=4.8065e22
#Pu comp factor= 4.78636e22

Found_Keff = 0
Flux_Matr=np.empty(8)
Core_Fission_Matr=np.empty(8)
def Enrichment_Update():
    Fuel_Trans_CS1= (1-Pu239_Frac)*4.3e-24*4.8065e22 + Pu239_Frac*4.5e-24*4.78636e22
    Fuel_Removal_CS1= (1-Pu239_Frac)*2.883e-24*4.8065e22 + Pu239_Frac*3.355e-24*4.78636e22
    Fuel_Fission_vCS1=(1-Pu239_Frac)*2.91*0.58e-24*4.8065e22 + Pu239_Frac*3.4*1.85e-24*4.78636e22
    Fuel_D1=1/(3*Fuel_Trans_CS1)
    #For group 2  FIXED
    Fuel_Trans_CS2= (1-Pu239_Frac)*4.8e-24*4.8065e22 + Pu239_Frac*5.1e-24*4.78636e22
    Fuel_Removal_CS2= (1-Pu239_Frac)*1.78e-24*4.8065e22 + Pu239_Frac*2.676e-24*4.78636e22
    Fuel_Fission_vCS2=(1-Pu239_Frac)*2.58*0.2e-24*4.8065e22 + Pu239_Frac*3.07*1.82e-24*4.78636e22
    Fuel_D2=1/(3*Fuel_Trans_CS2)
    #Below this there's no U portion for fission vCS as its fission CS is zero. Group 3 FIXED
    Fuel_Trans_CS3= (1-Pu239_Frac)*6.3e-24*4.8065e22 + Pu239_Frac*6.3e-24*4.78636e22
    Fuel_Removal_CS3= (1-Pu239_Frac)*0.4859e-24*4.8065e22 + Pu239_Frac*2.0809e-24*4.78636e22
    Fuel_Fission_vCS3=Pu239_Frac*2.95*1.6e-24*4.78636e22
    Fuel_D3=1/(3*Fuel_Trans_CS3)
    #Group 4
    Fuel_Trans_CS4= (1-Pu239_Frac)*9.3e-24*4.8065e22 + Pu239_Frac*8.6e-24*4.78636e22
    Fuel_Removal_CS4= (1-Pu239_Frac)*0.4435e-24*4.8065e22 + Pu239_Frac*1.9005e-24*4.78636e22
    Fuel_Fission_vCS4= Pu239_Frac*2.9 *1.58e-24*4.78636e22
    Fuel_D4=1/(3*Fuel_Trans_CS4)
    #Group 5
    Fuel_Trans_CS5= (1-Pu239_Frac)*11.7e-24*4.8065e22 + Pu239_Frac*11.3e-24*4.78636e22
    Fuel_Removal_CS5= (1-Pu239_Frac)*0.46e-24*4.8065e22 + Pu239_Frac*2.1e-24*4.78636e22
    Fuel_Fission_vCS5= Pu239_Frac*2.88 * 1.6e-24*4.78636e22
    Fuel_D5=1/(3*Fuel_Trans_CS5)
    #Group 6
    Fuel_Trans_CS6= (1-Pu239_Frac)*12.7e-24*4.8065e22 + Pu239_Frac*13.1e-24*4.78636e22
    Fuel_Removal_CS6= (1-Pu239_Frac)*0.56e-24*4.8065e22 + Pu239_Frac*2.35e-24*4.78636e22
    Fuel_Fission_vCS6= Pu239_Frac*2.88 *1.67e-24*4.78636e22
    Fuel_D6=1/(3*Fuel_Trans_CS6)
    #Group 7
    Fuel_Trans_CS7= (1-Pu239_Frac)*13.1e-24*4.8065e22 + Pu239_Frac*16.5e-24*4.78636e22
    Fuel_Removal_CS7= (1-Pu239_Frac)*0.85e-24*4.8065e22 + Pu239_Frac*4.77e-24*4.78636e22
    Fuel_Fission_vCS7= Pu239_Frac*2.87 * 2.78e-24*4.78636e22
    Fuel_D7=1/(3*Fuel_Trans_CS7)
    #for group 8
    Fuel_Trans_CS8= (1-Pu239_Frac)*11e-24*4.8065e22 + Pu239_Frac*31.8e-24*4.78636e22
    Fuel_Removal_CS8= (1-Pu239_Frac)*1.47e-24*4.8065e22 + Pu239_Frac*19.17e-24*4.78636e22
    Fuel_Fission_vCS8= Pu239_Frac*2.87 * 10.63e-24*4.78636e22
    Fuel_D8=1/(3*Fuel_Trans_CS8)
    
    #now create matrix, needs macro values, fraction*macro*compensation, for both U and Pu in fuel
    Fuel_Scat_Matr=np.array([[Fuel_D1*Bg2+Fuel_Removal_CS1,0,0,0,0,0,0,0],
                    [-((1-Pu239_Frac)*1.28e-24*4.8065e22+Pu239_Frac*0.66e-24*4.78636e22),Fuel_D2*Bg2+Fuel_Removal_CS2,0,0,0,0,0,0],
                    [-((1-Pu239_Frac)*0.78e-24*4.8065e22+Pu239_Frac*0.6e-24*4.78636e22),-((1-Pu239_Frac)*1.05e-24*4.8065e22+Pu239_Frac*0.64e-24*4.78636e22),Fuel_D3*Bg2+Fuel_Removal_CS3,0,0,0,0,0],
                    [-((1-Pu239_Frac)*0.2e-24*4.8065e22+Pu239_Frac*0.19e-24*4.78636e22),-((1-Pu239_Frac)*0.42e-24*4.8065e22+Pu239_Frac*0.15e-24*4.78636e22),-((1-Pu239_Frac)*0.33e-24*4.8065e22+Pu239_Frac*0.31e-24*4.78636e22),Fuel_D4*Bg2+Fuel_Removal_CS4,0,0,0,0],
                    [-((1-Pu239_Frac)*0.03e-24*4.8065e22+Pu239_Frac*0.04e-24*4.78636e22),-((1-Pu239_Frac)*0.01e-24*4.8065e22+Pu239_Frac*0.03e-24*4.78636e22),-((1-Pu239_Frac)*0.04e-24*4.8065e22+Pu239_Frac*0.05e-24*4.78636e22),-((1-Pu239_Frac)*0.29e-24*4.8065e22+Pu239_Frac*0.18e-24*4.78636e22),Fuel_D5*Bg2+Fuel_Removal_CS5,0,0,0],
                    [-((1-Pu239_Frac)*0.003e-24*4.8065e22+Pu239_Frac*0.005e-24*4.78636e22),-((1-Pu239_Frac)*0.01e-24*4.8065e22+Pu239_Frac*0.006e-24*4.78636e22),-((1-Pu239_Frac)*0.005e-24*4.8065e22+Pu239_Frac*0.01e-24*4.78636e22),-((1-Pu239_Frac)*0.003e-24*4.8065e22+Pu239_Frac*0.01e-24*4.78636e22),-((1-Pu239_Frac)*0.18e-24*4.8065e22+Pu239_Frac*0.13e-24*4.78636e22),Fuel_D6*Bg2+Fuel_Removal_CS6,0,0],
                    [0,0,-((1-Pu239_Frac)*0.0009e-24*4.8065e22+Pu239_Frac*0.0009e-24*4.78636e22),-((1-Pu239_Frac)*0.0005e-24*4.8065e22+Pu239_Frac*0.0005e-24*4.78636e22),-((1-Pu239_Frac)*0.02e-24*4.8065e22+Pu239_Frac*0.02e-24*4.78636e22),-((1-Pu239_Frac)*0.09e-24*4.8065e22+Pu239_Frac*0.09e-24*4.78636e22),Fuel_D7*Bg2+Fuel_Removal_CS7,0],
                    [0,0,0,0,0,0,-((1-Pu239_Frac)*0.01e-24*4.8065e22+Pu239_Frac*0.01e-24*4.78636e22),Fuel_D8*Bg2+Fuel_Removal_CS8]])
    
    Fuel_vCS_Matr=np.array([Fuel_Fission_vCS1,Fuel_Fission_vCS2,Fuel_Fission_vCS3,Fuel_Fission_vCS4,Fuel_Fission_vCS5,Fuel_Fission_vCS6,Fuel_Fission_vCS7,Fuel_Fission_vCS8])
    #empty core L matrix so it can be populated by loop
    Core_Scat_Matr=np.empty((8,8))
    Core_Fission_Matr=np.array([Fuel_Fission_vCS1,Fuel_Fission_vCS2,Fuel_Fission_vCS3,Fuel_Fission_vCS4,Fuel_Fission_vCS5,Fuel_Fission_vCS6,Fuel_Fission_vCS7,Fuel_Fission_vCS8])
    #Now homogenize properties for all the 8x8 matrices
    for i in range(8):
        #I moves down, j moves right. Sweeps through matrices and homogenizes properties
        Core_Fission_Matr[i] = Core_Fission_Matr[i] * Fuel_Frac * Fuel_Flux_Dis[i]
        for j in range(8):
            #averages values across all three section for volume
            Core_Scat_Matr[i][j] = Lead_Scat_Matr[i][j]*Cool_Frac*Cool_Flux_Dis[i] + Iron_Scat_Matr[i][j]*Clad_Frac*Clad_Flux_Dis[i] + Fuel_Scat_Matr[i][j] * Fuel_Frac*Fuel_Flux_Dis[i]
            
            
    #Got all that, now solve for keff and flux
    Flux_Matr=np.linalg.solve(Core_Scat_Matr,Chi_Matr)
    #Find Keff value
    Found_Keff = 0
    for k in range (8):
       Found_Keff += Flux_Matr[k]*Core_Fission_Matr[k]
    return(Found_Keff)

Keff = Enrichment_Update()
#find and print out Keff at current set enrichment
print ("keff = " + str(Keff))
print("Number of cells: ", Num_Cells)
print("Fuel load =", Fuel_Load, "kg")
print("Pu-239 load =", Pu239_Load, "kg")
print(Flux_Matr)
