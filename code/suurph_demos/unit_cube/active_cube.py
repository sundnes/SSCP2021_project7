import matplotlib.pyplot as plt
import numpy as np
from fenics import *

from guccionematerial import GuccioneMaterial

"""
This version of the simple contracting cube demo sallows two different boundary
conditions on the right side of the cube:
- A constant pressure load (zero or non-zero)
- A Robin boundary condition, i.e., a linear spring attached to the entire
right boundary. Spring constant can be adjusted to allow free contraction or
restrict the shortening of the cube.
Switching between the versions require editing the code and adding the lines
that are now "commented out"
"""


# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"] = 4

# Setup the mesh and the function space for the solutions
mesh = UnitCubeMesh(4,4,4)
V = VectorFunctionSpace(mesh, "Lagrange", 2)


# Define functions
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration

# Mark boundary subdomains
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)

boundary_markers = MeshFunction("size_t", mesh,mesh.topology().dim() - 1)
boundary_markers.set_all(0)
left.mark(boundary_markers, 1)
right.mark(boundary_markers, 2)

# Redefine boundary measure
ds = Measure('ds',domain=mesh,subdomain_data=boundary_markers)

# Define Dirichlet boundary (x = 0 or x = 1)
clamp = Constant((0.0, 0.0, 0.0))
bc = DirichletBC(V, clamp, left)
bcs = [bc]

# Kinematics
d = len(u)
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
#C = F.T*F                   # the right Cauchy-Green tensor
#E = 0.5*(C - I)             # the Green-Lagrange strain tensor

F = variable(F)

# Tissue microstructure
f0 = as_vector([ 1.0, 0.0, 0.0 ])
s0 = as_vector([ 0.0, 1.0, 0.0 ])
n0 = as_vector([ 0.0, 0.0, 1.0 ])

mat = GuccioneMaterial(e1=f0,e2=s0,e3=n0,kappa=1e3,Tactive=0.0)
psi = mat.strain_energy(F)

P = diff(psi,F) # first Piola-Kirchhoff stress tensor

#p_right = Constant(0.0) #the pressure load (zero for now)

#for Robin BC:
bc_spring = Constant(0.0)   #0.0 means free contraction, set to, e.g
#bc_spring = Constant(20.0) #Cube contracts against a linear spring

# Definition of the weak form:
N = FacetNormal(mesh)
#G_neumann = p_right * inner(v, det(F)*inv(F)*N) * ds(2) #ds(2) = left boundary
G_robin = inner(bc_spring * u, v) * ds(2)
R = inner(P,grad(v))*dx + G_robin #G_neumann

# The middle point on the right boundary
point0 = np.array([1.0,0.5,0.5])

# Step-wise loading (for plotting and convergence)
load_steps = 6
target_load = 5.0
loads = np.linspace(0,target_load,load_steps)

d0 = np.zeros(3)                #displacement at point0
disp = np.zeros(load_steps) #array to store displacement for all steps

disp_file = XDMFFile("results/guccione_active/u.xdmf")
for step in range(load_steps):
    # Stretch is a negative pressure
    mat.set_active_stress(loads[step])
    #p_right.assign(-loads[step])

    #solve the problem:
    solve(R == 0, u, bcs)

    #evaluate displacement at point defined above
    u.eval(d0,point0)
    disp[step] = d0[0]
    disp_file.write_checkpoint(u, "Displacement", step, append=True)

disp_file.close()

#displacement on x-axis, load on y-axis
plt.plot(loads,disp)
plt.xlabel('Active tension in $x$-direction')
plt.ylabel('Displacement of point (1.0,0.5,0.5)')

plt.show()
