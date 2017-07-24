import numpy as np
import collections
import sys
from optparse import OptionParser
usage = "usage: %prog  output_file"
parser = OptionParser(usage=usage)

(options, args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(1)

filename = str(args[0])

######################################################################

def write_to_VTK(filename, origin, spacing, dimensions, variables):
    fout = open("%s"%filename, 'wb')
    fout.write(bytes("# vtk DataFile Version 3.0\n", 'utf-8'))
    fout.write(bytes("vtk output\n", 'utf-8'))
    fout.write(bytes("ASCII\n", 'utf-8'))
    fout.write(bytes("DATASET STRUCTURED_POINTS\n", 'utf-8'))
    fout.write(bytes("DIMENSIONS {0} {1} {2}\n".format(*dimensions), 'utf-8'))
    fout.write(bytes("ORIGIN {0} {1} {2}\n".format(*origin), 'utf-8'))
    fout.write(bytes("SPACING {0} {1} {2}\n".format(*spacing), 'utf-8'))
    fout.write(bytes("POINT_DATA {0}\n".format(np.prod(dimensions)), 'utf-8'))
    for k,v in variables.items():
        fout.write(bytes("\nSCALARS %s float\n" %k, 'utf-8'))
        fout.write(bytes("LOOKUP_TABLE default\n", 'utf-8'))
        for data_slice in v.T:
            np.savetxt(fout, data_slice, fmt='%-15.7f')
    fout.close()

######################################################################
#                              BEGIN                                 #
######################################################################

# Grid characteristics
origin = [0, 0, 0]

Lx = 1.
scale = 1
nx = 3000

dx = float(Lx/(nx-1))
ny = int(scale*(nx-1)+1)
Ly = scale * Lx
nz = 1

spacing = [dx, dx, dx]
dimensions = [nx, ny, nz]

Pref = 101325.
P = np.zeros((nx, ny, nz))
ux  = np.zeros((nx, ny, nz))
uy  = np.zeros((nx, ny, nz))
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            x = i*spacing[0]
            y = j*spacing[1]

# Double Shear layer (Minion et Brown)
            if (y <= 0.5):
                P[i,j,k]  = Pref
                ux[i,j,k]   = 10.*np.tanh(80.*(y-0.25))
                uy[i,j,k]   = 10.*0.05*np.sin(2.*np.pi*(x+0.25))
            else:
                P[i,j,k]  = Pref
                ux[i,j,k]   = 10.*np.tanh(80.*(0.75-y))
                uy[i,j,k]   = 10.*0.05*np.sin(2.*np.pi*(x+0.25))

variables = collections.OrderedDict()
variables['Pressure']   = P
variables['velocity_X'] = ux
variables['velocity_Y'] = uy

write_to_VTK(filename, origin, spacing, dimensions, variables)


