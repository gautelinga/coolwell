import dolfin as df
import argparse
#import matplotlib.pyplot as plt
import numpy as np
import math
from mesh_box import square_cavity_mesh
import os
from bcs import PeriodicBC, Wall, TopWall, BtmWall, NotWall
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
folder = "data_square"


# Form compiler options
df.parameters["form_compiler"]["optimize"] = True
df.parameters["form_compiler"]["cpp_optimize"] = True

parser = argparse.ArgumentParser(description="Taylor dispersion")
parser.add_argument("-res", type=int, default=32, help="Resolution")
parser.add_argument("-Lx", type=float, default=2.0, help="Length x")
parser.add_argument("-b", type=float, default=1.0, help="Roughness b")
parser.add_argument("-H", type=float, default=0.5, help="Height h")
parser.add_argument("-fx", type=float, default=1.0, help="Force x")
parser.add_argument("-Ttop", type=float, default=1.0, help="T top")
parser.add_argument("-Tbtm", type=float, default=0.0, help="T btm")
parser.add_argument("-alpha", type=float, default=0.1,
                    help="gravity strength")
parser.add_argument("-nutop", type=float, default=1.0, help="nu top")
parser.add_argument("-nubtm", type=float, default=2.0, help="nu btm")
parser.add_argument("-kappa", type=float, default=1.0,
                    help="thermal diffusivity")

#parser.add_argument("--onlyflow", action="store_true", help="only flow")
args = parser.parse_args()

N = args.res
Lx = args.Lx
b = args.b
h = args.H

dx = b/N

Ttop = df.Constant(args.Ttop)
Tbtm = df.Constant(args.Tbtm)
T0 = 0.5*(Ttop + Tbtm)
dT = Ttop-Tbtm
nutop = df.Constant(args.nutop)
nubtm = df.Constant(args.nubtm)
alpha = df.Constant(args.alpha)
kappa = df.Constant(args.kappa)

# Change this
mesh = square_cavity_mesh(Lx, h, b, dx)
coords = mesh.coordinates()[:]

Eu = df.VectorElement("Lagrange", mesh.ufl_cell(), 2)
Ep = df.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
ET = df.FiniteElement("Lagrange", mesh.ufl_cell(), 1)

pbc = PeriodicBC(Lx, b, h)
wall = Wall(Lx, b, h)
notwall = NotWall(Lx, b, h)
topwall = TopWall(Lx, b, h)
btmwall = BtmWall(Lx, b, h)

subd = df.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
subd.set_all(0)
wall.mark(subd, 1)
topwall.mark(subd, 2)
btmwall.mark(subd, 3)
notwall.mark(subd, 0)

if rank == 0 and not os.path.exists(folder):
    os.makedirs(folder)

with df.XDMFFile(mesh.mpi_comm(), "{}/subd_b{}.xdmf".format(
        folder, args.b)) as xdmff:
    xdmff.write(subd)

W = df.FunctionSpace(mesh, df.MixedElement([Eu, Ep, ET]),
                     constrained_domain=pbc)
#S = df.FunctionSpace(mesh, ET, constrained_domain=pbc)

w_ = df.Function(W)
u_, p_, T_ = df.split(w_)
w = df.TrialFunction(W)
u, p, T = df.split(w)
v, q, psi = df.TestFunctions(W)

nu_ = (nutop*(T_-Tbtm) + nubtm * (Ttop-T_))/dT

f = df.Constant((args.fx, 0.)) + alpha * (T_ - T0) * df.Constant((0., 1.))

F = (
    df.inner(df.grad(u_)*u_, v)*df.dx
    + 2*nu_*df.inner(df.sym(df.grad(u_)), df.sym(df.grad(v)))*df.dx
    - df.div(v)*p_*df.dx - df.div(u_)*q*df.dx
    - df.dot(f, v)*df.dx
    - T_ * df.dot(u_, df.grad(psi)) * df.dx
    + kappa * df.dot(df.grad(T_), df.grad(psi)) * df.dx
)

J = df.derivative(F, w_, du=w)

bcuwall = df.DirichletBC(W.sub(0), df.Constant((0., 0.)),
                         subd, 1)
bcutop = df.DirichletBC(W.sub(0), df.Constant((0., 0.)),
                        subd, 2)
bcubtm = df.DirichletBC(W.sub(0), df.Constant((0., 0.)),
                        subd, 3)
bcTtop = df.DirichletBC(W.sub(2), Ttop, subd, 2)
bcTbtm = df.DirichletBC(W.sub(2), Tbtm, subd, 3)

x0, y0 = coords[0, 0], coords[0, 1]
x0 = comm.bcast(x0, root=0)
y0 = comm.bcast(y0, root=0)
# distribute!

bcp = df.DirichletBC(W.sub(1), df.Constant(0.),
                     ("abs(x[0]-({x0})) < DOLFIN_EPS && "
                      "abs(x[1]-({y0})) < DOLFIN_EPS").format(x0=x0, y0=y0),
                     "pointwise")

bcs = [bcuwall, bcutop, bcubtm, bcp, bcTtop, bcTbtm]

problem = df.NonlinearVariationalProblem(F, w_, bcs=bcs, J=J)
solver = df.NonlinearVariationalSolver(problem)

solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-14
solver.parameters["newton_solver"]["krylov_solver"]["absolute_tolerance"] = 1e-14
solver.parameters["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = True

solver.solve()

one_well = df.interpolate(df.Expression('x[1] < DOLFIN_EPS ? 1. : 0.', degree=0), W.sub(2).collapse())
V_Omega_well = df.assemble(one_well*df.dx)
one_channel = df.interpolate(df.Expression('x[1] < DOLFIN_EPS ? 0. : 1.', degree=0), W.sub(2).collapse()) 
V_Omega_channel = df.assemble(one_channel*df.dx)

u_b = df.assemble(df.sqrt(df.dot(u_,u_))*one_well*df.dx)/V_Omega_well
u_0 = df.assemble(df.sqrt(df.dot(u_,u_))*one_channel*df.dx)/V_Omega_channel

if rank == 0:
    print("u_b = {}".format(u_b))
    print("u_0 = {}".format(u_0))
    if u_b != 0:
       print("u_0/u_b = {}".format(u_0/u_b))
   

nu_mean = (args.nutop+args.nubtm)*0.5
Re = args.b*u_b/nu_mean
Re_0 = args.b*u_0/nu_mean
Pr = nu_mean/args.kappa
Pe = Re*Pr
Ra = args.alpha*(args.Ttop-args.Tbtm)*b**3/(nu_mean*args.kappa)
RT = np.inf
if args.fx != 0:
    RT = args.alpha*(args.Ttop-args.Tbtm)/args.fx

ds = df.Measure("ds", domain=mesh, subdomain_data=subd)

qy = df.assemble(kappa*T_.dx(1)*ds(3))
qy_ref = args.kappa * (args.Ttop-args.Tbtm) / (args.b+args.H) * args.b
Nu = qy/qy_ref

if rank == 0:
    print("Re = {}".format(Re))
    print("Re_0 = {}".format(Re_0))
    print("Ra = {}".format(Ra))
    print("Pe = {}".format(Pe))
    print("RT = {}".format(RT))
    print("Nu = {}".format(Nu))
    print("qy = {}".format(qy))
    print("Ri = {}".format(Ri))
    print("Ri*Re = {}".format(Ri*Re))

U_, P_, To_ = w_.split(deepcopy=True)
U_.rename("u", "tmp")
#U_.vector()[:] /= u_b
To_.rename("T", "tmp")

subfolder = "{}/RT{:.2e}_Re{:.2e}_Pe{:.2e}_b{}_H{}".format(folder, RT, Re, Pe, args.b, args.H)
if rank == 0 and not os.path.exists(subfolder):
    os.makedirs(subfolder)

with df.HDF5File(mesh.mpi_comm(),
                 "{}/flow.h5".format(subfolder), "w") as h5f:
    h5f.write(w_, "w")

with df.XDMFFile(mesh.mpi_comm(),
                 "{}/U.xdmf".format(subfolder)) as xdmff:
    xdmff.write(U_)

with df.XDMFFile(mesh.mpi_comm(),
                 "{}/T.xdmf".format(subfolder)) as xdmff:
    xdmff.write(To_)

if rank == 0:
   with open("{}/parameters.dat".format(subfolder),"w") as file:
      file.write("u_b = {}\nu_0 = {}\nRe_b = {}\nRe_0 = {}\nRa = {}\nPe = {}\nRT = {}\nNu = {}\nqy = {}\nRi = {}\nRi*Re = {}".format(u_b,u_0,Re,Re_0,Ra,Pe,RT,Nu,qy, Ri, Ri*Re))
      file.write("\nmpiexec -n 4 python3 steady.py -res {} -Lx {} -b {} -H {} -fx {} -Ttop {} -Tbtm {} -alpha {} -nutop {} -nubtm {} -kappa {}".format(args.res, args.Lx, args.b, args.H, args.fx, args.Ttop, args.Tbtm, args.alpha, args.nutop, args.nubtm, args.kappa))
   
