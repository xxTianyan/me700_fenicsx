from dolfinx import mesh, fem, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import ufl
import numpy as np

# Define beam geometry and material parameters
L = 1         # Beam length
W = 0.2       # Beam width/height
mu = 1        # Shear modulus
rho = 1       # Density
delta = W / L
gamma = 0.4 * delta**2  # Gravity-like term
beta = 1.25
lambda_ = beta           # First Lamé parameter
g = gamma                # Gravitational acceleration

# Create 3D hexahedral mesh for the domain: [0, L] x [0, W] x [0, W]
domain = mesh.create_box(MPI.COMM_WORLD,
                         [np.array([0, 0, 0]), np.array([L, W, W])],
                         [20, 6, 6],
                         cell_type=mesh.CellType.hexahedron)

# Define vector-valued Lagrange finite element space (P1 elements)
V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim, )))

# Function to mark clamped boundary (x = 0)
def clamped_boundary(x):
    return np.isclose(x[0], 0)

# Locate boundary facets and apply zero displacement (clamped)
fdim = domain.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(domain, fdim, clamped_boundary)
u_D = np.array([0, 0, 0], dtype=default_scalar_type)
bc = fem.dirichletbc(u_D,
                     fem.locate_dofs_topological(V, fdim, boundary_facets),
                     V)

# Define zero traction (no Neumann BC)
T = fem.Constant(domain, default_scalar_type((0, 0, 0)))
ds = ufl.Measure("ds", domain=domain)

# Define strain and stress tensors
def epsilon(u):
    return ufl.sym(ufl.grad(u))

def sigma(u):
    return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

# Define trial and test functions
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

# Define body force (gravity in -z direction)
f = fem.Constant(domain, default_scalar_type((0, 0, -rho * g)))

# Define bilinear and linear forms
a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
L = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds

# Solve linear elasticity problem
problem = LinearProblem(a, L, bcs=[bc],
                        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

# Evaluate solution at points along the beam centerline (y=z=0.1)
from dolfinx import geometry
bb_tree = geometry.bb_tree(domain, domain.topology.dim)
points = np.zeros((3, 101))
tol = 0.001
x = np.linspace(tol, 1 - tol, 101)
center = np.full((101,), 0.1)
points[0] = x
points[1] = center
points[2] = center

cells = []
points_on_proc = []
cell_candidates = geometry.compute_collisions_points(bb_tree, points.T)
colliding_cells = geometry.compute_colliding_cells(domain, cell_candidates, points.T)
for i, point in enumerate(points.T):
    if len(colliding_cells.links(i)) > 0:
        points_on_proc.append(point)
        cells.append(colliding_cells.links(i)[0])

points_on_proc = np.array(points_on_proc, dtype=np.float64)
u_values = uh.eval(points_on_proc, cells)

# Plot z-displacement vs x-position along the beam centerline
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(12,8))
ax.plot(points_on_proc[:, 0], u_values[:,2])
ax.grid()
ax.set_xlabel('length')
ax.set_ylabel('deflection')
plt.savefig('result.png', dpi=300)
