# Input parameter generation script for 2D CFD Problems
''' 
Written by Matt Blomquist
Last Update: 2017-12-04 (YYYY-MM-DD)

This script generates the input parameter file for 2D 
CFD problems.
'''

# Geometry Parameters
# Problem geometry :: x = length, y = width, z = depth
x = 1.0
y = 1.0
z = 1.0

# Material Properties
# Fluid Thermal Conductivity
k = 0.6

# Fluid Viscosity
mu = 0.89

# Fluid Density
rho = 1.0

# Reynolds Number of Flow
Re = 100.0

# Characteristic Length
l0 = 1.0

# Characteristic Velocity
v0 = 89.0

# Solution Properties
# Solution Tolerance
tol = 1.0e-6

# Max Iterations for SIMPLER Loop
itrmax = 10

# Boundary Conditions
#
# West Wall
bc_F_west_type = 'velocity'

# East Wall
bc_F_east_type = 'pressure'

# South Wall
bc_F_south_type = 'no-slip'

# North Wall
bc_F_north_type = 'no-slip' 

# Output paramters to :: input2d.txt 
# ---------------------------------------
# --------- DO NOT MODIFY BELOW ---------
# ---------------------------------------

# Generate File Name Number
case = '0001'
file_name = 'src/input2d_{}.txt'.format(case)

f = open(file_name,'w')

f.write('! 2D CFD Parameter Input File\n')
f.write('! Geometry parameters :: height, width, length\n')
f.write('{}, {}, {}\n'.format(x,y,z))
f.write('! Media parameters :: thermal conductivity, viscosity, density, Reynolds number\n')
f.write('{}, {}, {}, {}\n'.format(k, mu, rho, Re))
f.write('! Equation parameters :: length scale, velocity scale\n')
f.write('{}, {}\n'.format(l0, v0))
f.write('! Convergence parameters :: solution tolerance, maximum iterations\n')
f.write('{}, {}\n'.format(tol, itrmax))
f.close()
