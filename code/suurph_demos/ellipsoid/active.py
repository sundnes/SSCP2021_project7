import matplotlib.pyplot as plt
import numpy as np
from dolfin import *
from guccionematerial import *

parameters["form_compiler"]["quadrature_degree"] = 4



def load_ellipsoid_data():
    """Returns 4-tuple:
    mesh - the mesh,
    mf - MeshFunction defining boundary markers,
    numbering - dict of marking numbers,
    fibers - list of functions defining microstructure"""

    mesh = Mesh(MPI.comm_world, "data/mesh.xml")
    mf = MeshFunction("size_t", mesh, "data/facet_function.xml")

    numbering = {
        "BASE": 10,
        "ENDO": 30,
        "EPI": 40
    }

    # load fibers, sheet, cross_sheet data
    fiber_element = VectorElement(family="Quadrature",
                                     cell=mesh.ufl_cell(),
                                     degree=4,
                                     quad_scheme="default")
    fiber_space = FunctionSpace(mesh, fiber_element)
    fiber = Function(fiber_space, "data/fiber.xml")
    sheet = Function(fiber_space, "data/sheet.xml")
    cross_sheet = Function(fiber_space, "data/cross_sheet.xml")

    fibers = [fiber,sheet, cross_sheet]

    fiber_space = VectorFunctionSpace(mesh, "CG", 1)

    fiber_int = project(fiber, fiber_space)

    #fiber_file = File("fibers.pvd")
    #fiber_file << fiber_int

    return mesh, mf, numbering, fibers


def compute_cavity_volume(mesh,mf, numbering,u=None):
    X = SpatialCoordinate(mesh)
    N = FacetNormal(mesh)

    if u is not None:
        I = Identity(3) # the identity matrix
        F = I + grad(u) # the deformation gradient
        J = det(F)
        vol_form = (-1.0/3.0) * dot(X + u,
                                           J * inv(F).T * N)
    else:
        vol_form = (-1.0/3.0) * dot(X, N)


    ds = Measure('ds',domain=mesh,subdomain_data=mf)

    return assemble(vol_form*ds(numbering["ENDO"]))


parameters["form_compiler"]["cpp_optimize"] = True

mesh, boundary_markers, numbering, fibers = load_ellipsoid_data()

V = VectorFunctionSpace(mesh,'P',1)

# Redefine boundary measure
ds = Measure('ds',domain=mesh,subdomain_data=boundary_markers)

clamp = Constant((0.0, 0.0, 0.0))
bc = DirichletBC(V, clamp, boundary_markers, numbering["BASE"])
bcs = [bc]

# Define solution u and test function v
u = Function(V)
v = TestFunction(V)

# Define strain measures
I = Identity(3) # the identity matrix
F = I + grad(u) # the deformation gradient
F = variable(F)

mat = GuccioneMaterial(e1=fibers[0],e2=fibers[1],e3=fibers[2],kappa=1e3,Tactive=0.0)

psi = mat.strain_energy(F)

#S = diff(W, E) # the second Piola-Kirchoff stress tensor
P = diff(psi,F)        # the first Piola-Kirchoff stress tensor

p_endo = Constant(0.0)

# Define nonlinear problem
N = FacetNormal(mesh)
Gext = p_endo * inner(v, det(F)*inv(F)*N)  * ds(numbering["ENDO"]) #_bottom
R = inner(P,grad(v))*dx + Gext #dot(T,v)*ds(2)
#J = derivative(R,u,du)

# Step-wise loading (for plotting and convergence)
pressure_steps = 20
active_steps = 20
target_pressure = 10.0
target_active = 20.0

#first ramp up pressure, then keep constant
filling_pressure = np.linspace(0,target_pressure,pressure_steps)
const_pressure = np.ones(active_steps)*target_pressure
pressures = np.concatenate((filling_pressure,const_pressure))

#zero active tension during filling, then increase linearly
active1 = np.zeros_like(filling_pressure)
active2 = np.linspace(0,target_active, active_steps)
active =  np.concatenate((active1, active2))

volumes = np.zeros_like(pressures)

disp_file = XDMFFile("results/inflate_contract/u.xdmf")
for step in range(pressure_steps+active_steps):
    p_endo.assign(pressures[step])
    mat.set_active_stress(active[step])
    solve(R == 0, u, bcs)
    volumes[step] = compute_cavity_volume(mesh,boundary_markers,numbering,u)
    disp_file.write_checkpoint(u, "Displacement", step, append=True)

disp_file.close()


plt.plot(volumes,pressures)
plt.xlabel('Volume')
plt.ylabel('Pressure')
plt.axis([2.0,5.5,0.0,11.0])
plt.show()
