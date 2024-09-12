import argparse
from generate_mesh import generate_mesh_from_points
from generate_fiber2D_biv import generate_fiber2D
from generate_mesh import *

import os

parser = argparse.ArgumentParser()
parser.add_argument('-epi', type=str, default='epi.txt', help='File with segmentation epicardium points')
parser.add_argument('-vd', type=str, default='vd.txt', help='File with segmentation right ventricle points')
parser.add_argument('-ve', type=str, default='ve.txt', help='File with segmentation left ventricle points')
parser.add_argument('-numfib', type=int, default=0, help='Number of fibroses')
parser.add_argument('-fibbase', type=str, default='f', help='Base name for fibrosis files. Files must start with 0 id.')
parser.add_argument('-o', type=str, default='output', help='Output file name')
parser.add_argument('-t', type=int, default=0, help='Choose input type: txt or matlab')
parser.add_argument('-m', type=str, default='0', help='matlab file name')
parser.add_argument('-s', type=int, default=5, help='Slice number')
parser.add_argument('-dx', type=float, default=0.1, help='dx')
parser.add_argument('-dy', type=float, default=0.1, help='dy')
parser.add_argument('-dz', type=float, default=0.1, help='dz')

args = parser.parse_args()

if args.t == 0:
    if args.fibbase == "-1":
        generate_mesh_from_points(args.epi,args.vd,args.ve,0, 0, args.o)
    else:
        generate_mesh_from_points(args.epi,args.vd,args.ve,args.fibbase, args.numfib, args.o)
    generate_fiber2D(args.o, args.numfib)

    #os.chdir('hexa-mesh-from-VTK/')
    #os.system(f'./bin/HexaMeshFromVTK -i "../outputs_other/{args.o}.vtu" --dx {args.dx} --dy {args.dy} --dz {args.dz} -r 1000 -o "../outputs_alg/{args.o}.alg" -c ../config_file.ini --2d')

elif args.t == 1:
    numfib = generate_mesh_from_matlab(args.m, args.o, args.s)
    generate_fiber2D(args.o, numfib)

    #os.chdir('hexa-mesh-from-VTK/')
    #os.system(f'./bin/HexaMeshFromVTK -i "../outputs_other/{args.o}.vtu" --dx {args.dx} --dy {args.dy} --dz {args.dz} -r 1000 -o "../outputs_alg/{args.o}.alg" -c ../config_file.ini --2d')

else:
    
    zSize = verify_zSize(args.m)
    slices = verify_slices(args.m)

    for i, slice in enumerate(slices):

        numfib = generate_mesh_from_matlab(args.m, args.o+"_slice_"+str(slice), slice)
        generate_fiber2D(args.o+"_slice_"+str(slice), numfib)
        os.chdir('hexa-mesh-from-VTK/')
        os.system(f'./bin/HexaMeshFromVTK -i "../outputs_other/{args.o}_slice_{slice}.vtu" --dx {args.dx} --dy {args.dy} --dz {args.dz} -r 1000 -o "../outputs_alg/{args.o}_slice_{slice}.alg" -c ../config_file.ini --2d')
        os.chdir('../')