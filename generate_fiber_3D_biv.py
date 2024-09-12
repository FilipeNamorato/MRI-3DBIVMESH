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

def convert_xdmf_to_vtu(meshname):
    filename = meshname+".xdmf"

    t, point_data, cell_data = None, None, None

    #FAZER ALTERAÇÃO DE TEMPO PARA O 3D
    with meshio.xdmf.TimeSeriesReader(filename) as reader:
        points, cells = reader.read_points_cells()
        t, point_data, cell_data = reader.read_data(0)

    mesh = meshio.Mesh(points, cells, point_data=point_data, cell_data=cell_data,)
    mesh.write(meshname+".vtu", file_format="vtu",  )
        