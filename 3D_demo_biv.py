import dolfin as df
import ldrb
import argparse


def get_physical_region(fname):
    #generate physical groups to xdmf
    import pathlib as pl
    tetra_mesh_name = pl.Path(f"mesh_{fname}.xdmf")
    auxmesh = df.Mesh()
    with df.XDMFFile(tetra_mesh_name.as_posix()) as infile:
        infile.read(auxmesh)

    tecido = df.MeshFunction("size_t", auxmesh, 3)
    with df.XDMFFile(pl.Path(tetra_mesh_name).as_posix()) as f:
        f.read(tecido, "name_to_read")

    tecido.array()[:] = tecido.array()==2 
    tecido.rename("tissue", "tissue")

    import os 

    os.system(f"rm mesh_{fname}.*")
    os.system(f"rm triangle_mesh_{fname}.*")

    return tecido

def solve_laplace(mesh, boundary_markers, boundary_values, ldrb_markers):
    V = df.FunctionSpace(mesh, 'P', 1)

    u_rv, u_lv, u_epi = boundary_values

    bc1 = df.DirichletBC(V, u_rv, boundary_markers, ldrb_markers["rv"]) 
    bc2 = df.DirichletBC(V, u_lv, boundary_markers, ldrb_markers["lv"])
    bc3 = df.DirichletBC(V, u_epi, boundary_markers, ldrb_markers["epi"])

    bcs=[bc1, bc2 ,bc3]

    ds = df.Measure('ds', domain=mesh, subdomain_data=boundary_markers)
    dx = df.Measure('dx', domain=mesh)

    # Define variational problem
    u = df.TrialFunction(V)
    v = df.TestFunction(V)
    f = df.Constant(0.0)   
    a = df.dot(df.grad(u), df.grad(v))*dx  
    L = f*v*dx

    # Compute solution
    u = df.Function(V)
    df.solve(a == L, u, bcs, solver_parameters=dict(linear_solver='gmres', preconditioner='hypre_amg')) 

    return u


#mesh, ffun, markers = ldrb.gmsh2dolfin("malha3.msh")

#for m in markers:
#    print(m)

# Update the markers which are stored within the mesh


parser = argparse.ArgumentParser() #argumentos para cada segmentação
parser.add_argument('-meshname', type=str, default='patient', help='Gmsh file name without extension (.msh)')
argumentos = parser.parse_args()

#convert mesh to fenics format
meshname = argumentos.meshname
#mesh = df.Mesh(meshname + '.xml')

#markers = df.MeshFunction("size_t", mesh, meshname + '_facet_region.xml')

mesh, ffun, markers = ldrb.gmsh2dolfin(meshname + ".msh")

materials = df.MeshFunction("size_t", mesh, meshname + '_physical_region.xml')

V0 = df.FunctionSpace(mesh, "DG", 0) 
tecido = df.Function(V0) 
tecido.vector()[:] = materials.array()==2
tecido.rename("tecido", "tecido")


for m in markers:
    print(m)


print(markers)


ldrb_markers = {
    "base": markers["base"][0],
    "lv": markers["ve"][0],
    "epi": markers["epi"][0],
    "rv": markers["vd"][0]
}

#print('passou aqui')

# Choose space for the fiber fields
# This is a string on the form {family}_{degree}
fiber_space = "DG_0"

#--------
# Create field to define the action potential fenotype
# Solve Laplace problems with different boundary conditions
# u=1 on epicardium
phi_epi = solve_laplace(mesh, ffun, [0, 0, 1], ldrb_markers)
# u=1 on LV endocardium
phi_lv = solve_laplace(mesh, ffun, [0, 1, 0], ldrb_markers)
# u=1 on RV endocardium
phi_rv = solve_laplace(mesh, ffun, [1, 0, 0], ldrb_markers)

# Compute field with Laplace solutions
V = df.FunctionSpace(mesh, 'Lagrange', 1)
u = df.Function(V)
u.interpolate(df.Expression('-(epi + 2*rv*lv/(rv+lv) ) + 1', epi=phi_epi, rv=phi_rv, lv=phi_lv, degree=1))

bc3 = df.DirichletBC(V, 0, ffun, ldrb_markers["epi"])
bc3.apply(u.vector())

u.rename("fenotipo","fenotipo")
#--------

# Compute the microstructure
fiber, sheet, sheet_normal = ldrb.dolfin_ldrb(
    mesh=mesh,
    fiber_space=fiber_space,
    ffun=ffun,
    markers=ldrb_markers,
    alpha_endo_lv=30,  # Fiber angle on the LV endocardium
    alpha_epi_lv=-30,  # Fiber angle on the LV epicardium
    beta_endo_lv=0,  # Sheet angle on the LV endocardium
    beta_epi_lv=0,  # Sheet angle on the LV epicardium
    alpha_endo_sept=60,  # Fiber angle on the Septum endocardium
    alpha_epi_sept=-60,  # Fiber angle on the Septum epicardium
    beta_endo_sept=0,  # Sheet angle on the Septum endocardium
    beta_epi_sept=0,  # Sheet angle on the Septum epicardium
    alpha_endo_rv=80,  # Fiber angle on the RV endocardium
    alpha_epi_rv=-80,  # Fiber angle on the RV epicardium
    beta_endo_rv=0,  # Sheet angle on the RV endocardium
    beta_epi_rv=0,
)


fiber.rename("f_0","f_0")
sheet.rename("s_0","s_0")
sheet_normal.rename("n_0","n_0")


print("Salvando...")

with df.XDMFFile(mesh.mpi_comm(), meshname + ".xdmf") as xdmf:
    xdmf.parameters.update(
    {
        "functions_share_mesh": True,
        "rewrite_function_mesh": False
    })
    xdmf.write(mesh)
    xdmf.write(fiber, 0)
    xdmf.write(sheet, 0)
    xdmf.write(sheet_normal,0)
    xdmf.write(tecido, 0)
    xdmf.write(u, 0)

print("Done.")

# Store the results
#with df.HDF5File(mesh.mpi_comm(), "biv.h5", "w") as h5file:
#    h5file.write(fiber, "/fiber")
#    h5file.write(sheet, "/sheet")
#    h5file.write(sheet_normal, "/sheet_normal")


# You can also store files in XDMF which will also compute the fiber angle as scalars on the glyph to be visualised in Paraview. Note that these functions don't work (yet) using mpirun

# (These function are not tested in parallel)
#ldrb.fiber_to_xdmf(fiber, "biv_fiber_new1")
#ldrb.fiber_to_xdmf(sheet, "biv_sheet")
#ldrb.fiber_to_xdmf(sheet_normal, "biv_sheet_normal")

# ![_](_static/figures/biv_fiber.png)
# [Link to source code](https://github.com/finsberg/ldrb/blob/main/demos/demo_biv.py)
