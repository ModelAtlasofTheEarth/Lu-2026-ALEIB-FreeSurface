#!/usr/bin/env python
# coding: utf-8

# In[1]:


import petsc4py
from petsc4py import PETSc

import underworld3 as uw
from underworld3 import function
from underworld3.cython.petsc_discretisation import petsc_dm_find_labeled_points_local

import numpy as np
import sympy
import os
from datetime import datetime
import sys
import matplotlib.pyplot as plt

tRatio =  int(sys.argv[1])
use_fssa = sys.argv[2].lower() == 'true' 

u = uw.scaling.units
ndim = uw.scaling.non_dimensionalise
dim = uw.scaling.dimensionalise

# scaling 3: vel
half_rate = 1.0 * u.centimeter / u.year
model_length = 500. * u.kilometer
gravity = 9.81 * u.meter / u.second**2
bodyforce = 3300 * u.kilogram / u.metre**3 *gravity 

KL = model_length
Kt = KL / half_rate
KM = bodyforce * KL**2 * Kt**2

scaling_coefficients                    = uw.scaling.get_coefficients()
scaling_coefficients["[length]"] = KL
scaling_coefficients["[time]"] = Kt
scaling_coefficients["[mass]"]= KM


xmin, xmax = ndim(-250 * u.kilometer), ndim(250 * u.kilometer)
ymin, ymax = ndim(-500 * u.kilometer), ndim(0 * u.kilometer)

xres, yres = 50,50
dy = (ymax-ymin)/yres
dx = (xmax-xmin)/xres

if use_fssa:
    outputPath = "op_Ex_2DTopoRelax_Cartesain_FreeSurf_ALE_withFSSA0.5_yres{:n}_tRatio{:n}_noSwarm/".format(yres,tRatio)
else:
    outputPath = "op_Ex_2DTopoRelax_Cartesain_FreeSurf_ALE_noFSSA_yres{:n}_tRatio{:n}_noSwarm/".format(yres,tRatio)

if uw.mpi.rank == 0:
    if os.path.exists(outputPath):
        for i in os.listdir(outputPath):
            os.remove(outputPath+ i)
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

mesh = uw.meshing.StructuredQuadBox(elementRes=(int(xres), int(yres)), minCoords=(xmin, ymin), maxCoords=(xmax, ymax))      
mesh0 = uw.meshing.StructuredQuadBox(elementRes=(int(xres), int(yres)), minCoords=(xmin, ymin), maxCoords=(xmax, ymax)) 

# # dq2dq1
# v = uw.discretisation.MeshVariable("V", mesh, mesh.dim, degree=2)
# p = uw.discretisation.MeshVariable("P", mesh, 1, degree=1)

# q1dq0
v = uw.discretisation.MeshVariable("V", mesh, mesh.dim, degree=1,continuous=True)
p = uw.discretisation.MeshVariable("P", mesh, 1, degree=0,continuous=False)
timeField     = uw.discretisation.MeshVariable("time", mesh, 1, degree=1)

# botwall = petsc_dm_find_labeled_points_local(mesh.dm,"Bottom")
# topwall = petsc_dm_find_labeled_points_local(mesh.dm,"Top")

ratio = 1
D = np.abs(ymin)
Lambda = D/ratio
k = 2.0 * np.pi / Lambda
mu0 = ndim(1e21  * u.pascal * u.second)
g = ndim(gravity)
rho0 = ndim(3300* u.kilogram / u.metre**3)
drho = rho0-0.
w_m = ndim(10*u.kilometer)

tau0 = 2*k*mu0/drho/g
tau = (D*k+np.sinh(D*k)*np.cosh(D*k))/(np.sinh(D*k)**2)*tau0

def perturbation(x):
    return w_m * np.cos(2.*np.pi*(x)/Lambda)
deform_fn = w_m * sympy.cos(2.*np.pi*(mesh.X[0])/Lambda)

# def get_analytical(x0,load_time):
#     #tau = (D*k+np.sinh(D*k)*np.cosh(D*k))/(np.sinh(D*k)**2)*tau0
#     #A = -F0/k/tau0
#     #B = -F0/k/tau
#     #C = F0/tau
#     #E = F0/tau/np.tanh(D*k)
#     #phi = np.sin(k*x)*np.exp(-tmax/tau)*(A*np.sinh(k*z)+B*np.cosh(k*z)+C*z*np.sinh(k*z)+E*z*np.cosh(k*z))
#     w = w_m*np.exp(-load_time/tau)
#     return w

max_time = tau*4
dt_set = tau/tRatio
save_every = 1

Dz = uw.discretisation.MeshVariable("Dz", mesh, 1, degree=1)
diffuser = uw.systems.Poisson(mesh, Dz)
diffuser.constitutive_model = uw.constitutive_models.DiffusionModel
diffuser.constitutive_model.Parameters.diffusivity = 1.

#diffuser.add_essential_bc((ymax,), "Top")
#diffuser.add_essential_bc((deform_fn,), "Internal")
diffuser.add_essential_bc((deform_fn,), "Top")
diffuser.add_essential_bc((0.,),"Bottom")
diffuser.solve()
#dissym = Dz.sym * mesh.CoordinateSystem.unit_e_0 
#displacement = uw.function.evaluate(Dz.sym[0], mesh.X.coords)
displacement = np.zeros((mesh.X.coords.shape[0],mesh.dim))
displacement[:,-1] = uw.function.evaluate(Dz.sym[0], mesh.X.coords)[:,0,0]
mesh._deform_mesh(mesh.X.coords + displacement)

density_fn = rho0
visc_fn = mu0
ND_gravity = g

stokes = uw.systems.Stokes(mesh, velocityField=v, pressureField=p)
stokes.constitutive_model = uw.constitutive_models.ViscousFlowModel
stokes.bodyforce = sympy.Matrix([0, -1 * ND_gravity * density_fn])
stokes.constitutive_model.Parameters.shear_viscosity_0 = visc_fn
stokes.saddle_preconditioner = 1.0 / stokes.constitutive_model.Parameters.shear_viscosity_0
stokes.add_essential_bc((0.0,None), "Left")
stokes.add_essential_bc((0.0,None), "Right")
stokes.add_essential_bc((0.0,0.0), "Bottom")

# if uw.mpi.size == 1:
#     stokes.petsc_options['pc_type'] = 'lu'

stokes.tolerance = 1.0e-6
stokes.petsc_options["ksp_rtol"] = 1.0e-6
stokes.petsc_options["ksp_atol"] = 1.0e-6
stokes.petsc_options["snes_converged_reason"] = None
stokes.petsc_options["snes_monitor_short"] = None

if use_fssa:
    theta = 0.5*density_fn*ND_gravity*dt_set
    FSSA_traction = theta*mesh.Gamma.dot(v.sym) * mesh.Gamma
    stokes.add_natural_bc(FSSA_traction, "Top")


# In[2]:


# def _adjust_time_units(val):
#     """ Adjust the units used depending on the value """
#     if isinstance(val, u.Quantity):
#         mag = val.to(u.years).magnitude
#     else:
#         val = dim(val, u.years)
#         mag = val.magnitude
#     exponent = int("{0:.3E}".format(mag).split("E")[-1])

#     if exponent >= 9:
#         units = u.gigayear
#     elif exponent >= 6:
#         units = u.megayear
#     elif exponent >= 3:
#         units = u.kiloyears
#     elif exponent >= 0:
#         units = u.years
#     elif exponent > -3:
#         units = u.days
#     elif exponent > -5:
#         units = u.hours
#     elif exponent > -7:
#         units = u.minutes
#     else:
#         units = u.seconds
#     return val.to(units)

# import underworld3 as uw
# from underworld3 import function
# from underworld3.cython.petsc_discretisation import petsc_dm_find_labeled_points_local

# from scipy.interpolate import interp1d
# from enum import Enum
# import numpy as np

# class FreeSurfType(Enum):
#     """
#     free surface method type:

#     FreeSurfType.CartesianALE     ALE in StructuredQuadBox 
#     FreeSurfType.CartesianALEIB   ALE with internal boundary in StructuredQuadBox
#     FreeSurfType.CartesianALEIBSP ALE with internal boundary and surface processes in StructuredQuadBox

#     FreeSurfType.AnnulusALE     ALE in Annulus
#     FreeSurfType.AnnulusALEIB   ALE with internal boundary in Annulus


#     FreeSurfType.RegionalSphericalALE     ALE in RegionalSphericalox 
#     FreeSurfType.RegionalSphericalALEIB   ALE with internal boundary in RegionalSphericalBox

#     """

#     CartesianALE = 0
#     CartesianALEIB = 1
#     CartesianALEIBSP = 2
#     AnnulusALE = 3
#     AnnulusALEIB = 4
#     RegionalSphericalALE = 5
#     RegionalSphericalALEIB = 6


# class FreeSurfaceProcessor_Cartesian(object): 
#     def __init__(self,init_mesh,mesh,v,type = None,):
#         """
#         Parameters
#         ----------
#         _init_mesh : the original mesh
#         mesh : the updating model mesh
#         vel : the velocity field of the model
#         dt : dt for advecting the surface
#         """

#         self.init_mesh = init_mesh
#         self.Tmesh = uw.discretisation.MeshVariable("Tmesh", self.init_mesh, 1, degree=1)
#         self.Bmesh = uw.discretisation.MeshVariable("Bmesh", self.init_mesh, 1, degree=1)
#         self.mesh_solver = uw.systems.Poisson(self.init_mesh , u_Field=self.Tmesh)
#         self.mesh_solver.constitutive_model = uw.constitutive_models.DiffusionModel
#         self.mesh_solver.constitutive_model.Parameters.diffusivity = 1. 
#         self.mesh_solver.f = 0.0
#         self.mesh_solver.add_dirichlet_bc((0.,), "Bottom")

#         if type == None:
#             type = FreeSurfType.CartesianALE
#         if not isinstance(type, FreeSurfType):
#             raise ValueError("'type' must be an instance of 'FreeSurfType'")
#         self.type = type 
#         if self.type == FreeSurfType.CartesianALEIB or type == FreeSurfType.CartesianALEIBSP:
#             self.mesh_solver.add_dirichlet_bc((0.,), "Top")
#             self.mesh_solver.add_dirichlet_bc((self.Bmesh.sym[0],), "Internal")
#             self.interface = petsc_dm_find_labeled_points_local(self.init_mesh.dm,"Internal")
#         elif self.type == FreeSurfType.CartesianALE:
#             self.mesh_solver.add_dirichlet_bc((self.Bmesh.sym[0],), "Top")
#             self.interface = petsc_dm_find_labeled_points_local(self.init_mesh.dm,"Top")

#         self.mesh = mesh
#         self.v = v

#     def _advect_surface(self): 
#         with self.init_mesh.access(self.Bmesh):
#             self.Bmesh.data[:, 0] = self.mesh.X.coords[:, -1]
#             #print("CPU.no: %d topsiez: %d \n" %(uw.mpi.rank,self.top.size))
#             if self.interface.size > 0:
#                 if self.mesh.dim == 2:         
#                     coords = self.mesh.X.coords[self.interface]
#                     vel = self.veldata[self.interface]
#                     coords2 = coords + vel * self._dt
#                     f = interp1d(coords2[:,0], coords2[:,1], kind='cubic', fill_value='extrapolate')
#                     self.Bmesh.data[self.interface, 0] = f(coords[:,0])-coords[:,-1] 
#                 else:
#                     coords = self.mesh.data[self.interface]
#                     vel = self.veldata[self.interface]
#                     new_coords = coords + vel * self._dt
#                     mesh_kdt = uw.kdtree.KDTree(coords[:,0:2].copy(order='C'))
#                     mesh_kdt.build_index()
#                     values = mesh_kdt.rbf_interpolator_local(new_coords[:,0:2].copy(order='C'),new_coords[:,-1][:, np.newaxis].copy(order='C'), self.mesh.dim+1)
#                     del mesh_kdt
#                     self.Bmesh.data[self.interface, 0] = values[:,0]-coords[:,-1] 
#         uw.mpi.barrier()
#         self.init_mesh.update_lvec()

#     def solve(self,dt):
#         self._dt = dt
#         #self.veldata = uw.function.evaluate(self.v.sym, self.mesh.data)
#         # for v type = Q1
#         with self.mesh.access():
#             self.veldata = self.v.data
#         self._advect_surface()
#         self.mesh_solver.solve() 
#         with self.init_mesh.access():
#             new_mesh_coords = np.zeros([self.Tmesh.data.shape[0],self.mesh.dim])
#             new_mesh_coords[:,-1] = self.Tmesh.data[:,0].copy()  
#         return new_mesh_coords+self.mesh.X.coords



def _adjust_time_units(val):
    """ Adjust the units used depending on the value """
    if isinstance(val, u.Quantity):
        mag = val.to(u.years).magnitude
    else:
        val = dim(val, u.years)
        mag = val.magnitude
    exponent = int("{0:.3E}".format(mag).split("E")[-1])

    if exponent >= 9:
        units = u.gigayear
    elif exponent >= 6:
        units = u.megayear
    elif exponent >= 3:
        units = u.kiloyears
    elif exponent >= 0:
        units = u.years
    elif exponent > -3:
        units = u.days
    elif exponent > -5:
        units = u.hours
    elif exponent > -7:
        units = u.minutes
    else:
        units = u.seconds
    return val.to(units)

import underworld3 as uw
from underworld3 import function
from underworld3.cython.petsc_discretisation import petsc_dm_find_labeled_points_local

from scipy.interpolate import interp1d
from enum import Enum
import numpy as np

class FreeSurfType(Enum):
    """
    free surface method type:

    FreeSurfType.CartesianALE     ALE in StructuredQuadBox 
    FreeSurfType.CartesianALEIB   ALE with internal boundary in StructuredQuadBox
    FreeSurfType.CartesianALEIBSP ALE with internal boundary and surface processes in StructuredQuadBox

    FreeSurfType.AnnulusALE     ALE in Annulus
    FreeSurfType.AnnulusALEIB   ALE with internal boundary in Annulus


    FreeSurfType.RegionalSphericalALE     ALE in RegionalSphericalox 
    FreeSurfType.RegionalSphericalALEIB   ALE with internal boundary in RegionalSphericalBox

    """

    CartesianALE = 0
    CartesianALEIB = 1
    CartesianALEIBSP = 2
    AnnulusALE = 3
    AnnulusALEIB = 4
    RegionalSphericalALE = 5
    RegionalSphericalALEIB = 6


class FreeSurfaceProcessor_Cartesian(object): 
    def __init__(self,init_mesh,mesh,v,type = None,):
        """
        Parameters
        ----------
        _init_mesh : the original mesh
        mesh : the updating model mesh
        vel : the velocity field of the model
        dt : dt for advecting the surface
        """

        self.init_mesh = init_mesh
        self.Tmesh = uw.discretisation.MeshVariable("Tmesh", self.init_mesh, 1, degree=1)
        self.Bmesh = uw.discretisation.MeshVariable("Bmesh", self.init_mesh, 1, degree=1)
        self.mesh_solver = uw.systems.Poisson(self.init_mesh , u_Field=self.Tmesh)
        self.mesh_solver.constitutive_model = uw.constitutive_models.DiffusionModel
        self.mesh_solver.constitutive_model.Parameters.diffusivity = 1. 
        self.mesh_solver.f = 0.0
        self.mesh_solver.add_dirichlet_bc((0.,), "Bottom")

        if type == None:
            type = FreeSurfType.CartesianALE
        if not isinstance(type, FreeSurfType):
            raise ValueError("'type' must be an instance of 'FreeSurfType'")
        self.type = type 
        if self.type == FreeSurfType.CartesianALEIB or type == FreeSurfType.CartesianALEIBSP:
            self.mesh_solver.add_dirichlet_bc((0.,), "Top")
            self.mesh_solver.add_dirichlet_bc((self.Bmesh.sym[0],), "Internal")
            self.interface = petsc_dm_find_labeled_points_local(self.init_mesh.dm,"Internal")
        elif self.type == FreeSurfType.CartesianALE:
            self.mesh_solver.add_dirichlet_bc((self.Bmesh.sym[0],), "Top")
            self.interface = petsc_dm_find_labeled_points_local(self.init_mesh.dm,"Top")

        self.mesh = mesh
        self.v = v

    def _advect_surface(self): 
        # with self.init_mesh.access(self.Bmesh):
        #     self.Bmesh.data[:, 0] = self.mesh.X.coords[:, -1]
            #print("CPU.no: %d topsiez: %d \n" %(uw.mpi.rank,self.top.size))
        if self.interface.size > 0:
            if self.mesh.dim == 2:         
                coords1 = self.mesh.X.coords[self.interface]
                vel = uw.function.evaluate(self.v,coords1)[:,0,:].copy()
                #vel = self.veldata[self.interface]
                coords2 = coords1 + vel * self._dt
                f = interp1d(coords2[:,0], coords2[:,1], kind='cubic', fill_value='extrapolate')
                self.Bmesh.data[self.interface, 0] = f(coords1[:,0])-coords1[:,-1] 
            else:
                coords1 = self.mesh.data[self.interface]
                vel = uw.function.evaluate(self.v,coords1)[:,0,:].copy()
                new_coords = coords + vel * self._dt
                mesh_kdt = uw.kdtree.KDTree(coords[:,0:2].copy(order='C'))
                mesh_kdt.build_index()
                values = mesh_kdt.rbf_interpolator_local(new_coords[:,0:2].copy(order='C'),new_coords[:,-1][:, np.newaxis].copy(order='C'), self.mesh.dim+1)
                del mesh_kdt
                self.Bmesh.data[self.interface, 0] = values[:,0]-coords[:,-1] 
        uw.mpi.barrier()
        self.init_mesh.update_lvec()

    def solve(self,dt):
        self._dt = dt
        #self.veldata = uw.function.evaluate(self.v.sym, self.mesh.data)
        # for v type = Q1
        #with self.mesh.access():
        #self.veldata = self.v.data
        self._advect_surface()
        self.mesh_solver.solve() 
        #with self.init_mesh.access():
        new_mesh_coords = np.zeros([self.Tmesh.data.shape[0],self.mesh.dim])
        new_mesh_coords[:,-1] = self.Tmesh.data[:,0].copy()  
        return new_mesh_coords+self.mesh.X.coords


# In[3]:


step      = 0
max_steps = 5
time      = 0
dt        = 0

#from freesurface import FreeSurfaceProcessor_Cartesian
#from freesurface import FreeSurfType
freesuface = FreeSurfaceProcessor_Cartesian(mesh0,mesh,v,type=FreeSurfType(0))

while time < max_time+dt_set:   
    if uw.mpi.rank == 0:
        string = """Step: {0:5d} Model Time: {1:6.1f} dt: {2:6.1f} ({3})\n""".format(
        step, _adjust_time_units(time),
        _adjust_time_units(dt),
        datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        sys.stdout.write(string)
        sys.stdout.flush()

    #stokes.solve(zero_init_guess=False,_force_setup=True)
    stokes.solve(zero_init_guess=False)

    if step%save_every ==0:
        if uw.mpi.rank == 0:
            print(f'\nSave data:')
        with mesh.access(timeField):
            timeField.data[:,0] = dim(time, u.megayear).m
        #mesh.petsc_save_checkpoint(meshVars=[v, p, timeField], index=step, outputPath=outputPath)
        mesh.write_timestep("TR",meshUpdates=True,meshVars=[v, p, timeField],
                            outputPath=outputPath,index=step,)


    dt_solver = stokes.estimate_dt()
    dt = min(dt_solver,dt_set)

    new_mesh_coords=freesuface.solve(dt)
    mesh._deform_mesh(new_mesh_coords)

    step += 1
    time += dt


# In[ ]:




