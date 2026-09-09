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
ymin, ymax = ndim(-500 * u.kilometer), ndim(100 * u.kilometer)
yint = ndim(0. * u.kilometer)

xres, yres = 50,60
dy = (ymax-ymin)/yres
dx = (xmax-xmin)/xres
yresa,yresb =int(np.around((ymax-yint)/dy)),int(np.around((yint-ymin)/dy))

use_fssa = False
use_diff = False

tRatio =  int(sys.argv[1])

 
if use_fssa:
    outputPath = "op_Ex_2DTopoRelax_Cartesain_FreeSurf_EulerianMC_withFSSA0.5_yres{:n}_tRatio{:n}_noSwarm/".format(yres,tRatio)
else:
    outputPath = "op_Ex_2DTopoRelax_Cartesain_FreeSurf_EulerainMC_noFSSA_yres{:n}_tRatio{:n}_noSwarm/".format(yres,tRatio)

if uw.mpi.rank == 0:
    if os.path.exists(outputPath):
        for i in os.listdir(outputPath):
            os.remove(outputPath+ i)
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

mesh = uw.meshing.StructuredQuadBox(elementRes=(int(xres), int(yres)), minCoords=(xmin, ymin), maxCoords=(xmax, ymax))      
# mesh = uw.meshing.BoxInternalBoundary(elementRes=(xres,yres),zelementRes=(yresa,yresb),minCoords=(xmin,ymin),maxCoords=(xmax, ymax),zintCoord=yint,degree=1,qdegree=2)
# mesh0 = uw.meshing.BoxInternalBoundary(elementRes=(xres,yres),zelementRes=(yresa,yresb),minCoords=(xmin,ymin),maxCoords=(xmax, ymax),zintCoord=yint,degree=1,qdegree=2)

# # dq2dq1
# v = uw.discretisation.MeshVariable("V", mesh, mesh.dim, degree=2)
# p = uw.discretisation.MeshVariable("P", mesh, 1, degree=1)

# q1dq0
v = uw.discretisation.MeshVariable("V", mesh, mesh.dim, degree=1,continuous=True)
p = uw.discretisation.MeshVariable("P", mesh, 1, degree=0,continuous=False)
timeField     = uw.discretisation.MeshVariable("time", mesh, 1, degree=1)

#material_mesh = uw.discretisation.MeshVariable("M", mesh, 1, degree=1,continuous=True)
phi0 = uw.discretisation.MeshVariable(r"\phi_0", mesh, 1, degree=1,continuous=True)
#phi = uw.discretisation.MeshVariable("phi", mesh, vtype=uw.VarType.SCALAR, degree=2,continuous=False)


wRatio = 1
D = np.abs(ymin)
Lambda = D/wRatio
k = 2.0 * np.pi / Lambda
mu0 = ndim(1e21  * u.pascal * u.second)
g = ndim(gravity)
rho0 = ndim(3300* u.kilogram / u.metre**3)
drho = rho0-0.
w_m = ndim(10*u.kilometer)
ND_gravity = g

tau0 = 2*k*mu0/drho/g
tau = (D*k+np.sinh(D*k)*np.cosh(D*k))/(np.sinh(D*k)**2)*tau0

def perturbation(x):
    return w_m * np.cos(2.*np.pi*(x)/Lambda)
deform_fn = w_m * sympy.cos(2.*np.pi*(mesh.X[0])/Lambda)

max_time = tau*4
dt_set = tau/tRatio
save_every = 1

R0 = uw.discretisation.MeshVariable("r_0", mesh, vtype=uw.VarType.SCALAR, degree=2, continuous=False)
R0.data[:,0] = uw.function.evaluate(mesh.X[1], R0.coords)[:,0,0]


# In[2]:


from scipy.spatial import distance
# https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
def points_in_polygon(pts,polygon):
    pts = np.asarray(pts,dtype='float32')
    polygon = np.asarray(polygon,dtype='float32')
    contour2 = np.vstack((polygon[1:], polygon[:1]))
    test_diff = contour2-polygon
    mask1 = (pts[:,None] == polygon).all(-1).any(-1)
    m1 = (polygon[:,1] > pts[:,None,1]) != (contour2[:,1] > pts[:,None,1])
    slope = ((pts[:,None,0]-polygon[:,0])*test_diff[:,1])-(test_diff[:,0]*(pts[:,None,1]-polygon[:,1]))
    m2 = slope == 0
    mask2 = (m1 & m2).any(-1)
    m3 = (slope < 0) != (contour2[:,1] < polygon[:,1])
    m4 = m1 & m3
    count = np.count_nonzero(m4,axis=-1)
    mask3 = ~(count%2==0)
    mask = mask1 | mask2 | mask3
    return mask

#topwall = petsc_dm_find_labeled_points_local(mesh.dm,"Top")
def init_phi():
    with interfaceSwarm.access():
        interface_coords = interfaceSwarm.data
        index = uw.kdtree.KDTree(interfaceSwarm.data)
        index.build_index()
    # index = uw.kdtree.KDTree(interface_coords)
    # index.build_index()
    indices, dist_sqr, found = index.find_closest_point(phi0.coords)
    phi_values = np.sqrt(dist_sqr)

    point_leftwall = np.array([xmin,interface_coords[0,1]])
    point_rightwall = np.array([xmax,interface_coords[0,1]])
    point_righttop = np.array([xmax,ymax+dy])
    point_lefttop = np.array([xmin,ymax+dy])
    polygon = np.vstack([point_lefttop,point_leftwall,interface_coords,point_rightwall,point_righttop,point_lefttop])
    # polygon = np.concatenate((polygon,line),axis=0)
    # polygon = np.vstack([polygon, polygon[0]])
    mask = points_in_polygon(phi0.coords, polygon)

    phi_values[~mask] = -phi_values[~mask]
    phi0.data[:,0] = phi_values
    return

# from scipy.interpolate import interp1d
# def init_phi():
#     with interfaceSwarm.access():
#         interface_coords = interfaceSwarm.data
#         index = uw.kdtree.KDTree(interfaceSwarm.data)
#         index.build_index()
#     indices, dist_sqr, found = index.find_closest_point(R0.coords)
#     phi_values = np.sqrt(dist_sqr)
#     f = interp1d(interface_coords[:,0], interface_coords[:,1], kind='cubic', fill_value='extrapolate')
#     mask = f(R0.coords[:,0]) < R0.coords[:,1]
#     phi_values[~mask] = -phi_values[~mask]
#     R0.data[:,0] = phi_values
#     return


# In[3]:


x = np.linspace(xmin,xmax,xres*2+1)
y = perturbation(x)
interface_coords = np.ascontiguousarray(np.array([x,y]).T)
interfaceSwarm = uw.swarm.Swarm(mesh)
interfaceSwarm.add_particles_with_coordinates(interface_coords) 
init_phi()


# In[4]:


# #topwall = petsc_dm_find_labeled_points_local(mesh.dm,"Top")

# #polygon = interface_coords
# import numpy as np
# #from skimage.measure import points_in_poly
# import matplotlib.pyplot as plt


# # line = mesh.X.coords[topwall]
# # x1,y1 = line[:,0],line[:,1]
# # zipxy = zip(x1,y1)
# # zipxy = sorted(zipxy,reverse=True)
# # x2,y2 = zip(*zipxy) 
# # line[:,0] = x2 
# # line[:,1] = y2
# point_leftwall = np.array([xmin,interface_coords[0,1]])
# point_rightwall = np.array([xmax,interface_coords[0,1]])
# point_righttop = np.array([xmax,ymax])
# point_lefttop = np.array([xmin,ymax])
# polygon = np.vstack([point_lefttop ,point_leftwall,interface_coords,point_rightwall,point_righttop,point_lefttop])
# #polygon = np.concatenate((polygon,line),axis=0)
# #polygon = np.vstack([polygon, polygon[0]])

# fig, ax1 = plt.subplots(nrows=1, figsize=(5,5))
# ax1.plot(*polygon.T, color="red")
# ax1.scatter(interface_coords[:,0],interface_coords[:,1],c="blue")
# # ax1.set_xlim([-5,5])
# # ax1.set_ylim([-6.6,0.])
# #plt.show()


# In[5]:


# # R0 = uw.discretisation.MeshVariable("r_0", mesh, vtype=uw.VarType.SCALAR, degree=2, continuous=False)
# # R0.data[:,0] = uw.function.evaluate(mesh.X[1], R0.coords)[:,0,0]

# viscM = ndim(1e21 * u.pascal * u.second)
# densityM = ndim(3300 * u.kilogram / u.metre**3)
# viscA = ndim(1e18 * u.pascal * u.second)
# densityA = ndim(0. * u.kilogram / u.metre**3)
# # density_fn = sympy.Piecewise((densityA,R0.sym[0]>0),
# #                              (densityM, True))
# # visc_fn = sympy.Piecewise((viscA,R0.sym[0]>0),
# #                           (viscM, True))


# use_diff = True
# alphah = 1*dy
# def material_parameter_fn(c1,c2,alphah):
#     return sympy.Piecewise((c1,R0.sym[0] <= -alphah),
#                            (c2,R0.sym[0] > alphah),
#                            ((c2-c1)*R0.sym[0]/alphah/2.+(c1+c2)/2, True))
# if use_diff:
#     visc_fn = material_parameter_fn(viscM,viscA,alphah)
#     density_fn = material_parameter_fn(densityM,densityA,alphah)


# In[6]:


AIndex = 0
MIndex = 1

viscM = ndim(1e21 * u.pascal * u.second)
densityM = ndim(3300 * u.kilogram / u.metre**3)
viscA = ndim(1e18 * u.pascal * u.second)
densityA = ndim(0. * u.kilogram / u.metre**3)


alphah = 1*dy
def material_parameter_fn(c1,c2,alphah):
    return sympy.Piecewise((c1,phi0.sym[0] <= -alphah),
                           (c2,phi0.sym[0] > alphah),
                           ((c2-c1)*phi0.sym[0]/alphah/2.+(c1+c2)/2, True))
if use_diff:
    visc_fn = material_parameter_fn(viscM,viscA,alphah)
    density_fn = material_parameter_fn(densityM,densityA,alphah)
    #bodyforce_fn = -(material_parameter_fn(densityM*ND_gravity,densityA*ND_gravity,alphah))*mesh.CoordinateSystem.unit_e_1
    print("use_diff")

else:
    visc_fn = sympy.Piecewise((viscA,phi0.sym[0] > 0.0),(viscM, True))
    density_fn = sympy.Piecewise((densityA,phi0.sym[0] > 0.0),(densityM, True))  
    #bodyforce_fn = -(sympy.Piecewise((densityA*ND_gravity,phi0.sym[0] > 0.0),(densityM*ND_gravity, True)))*mesh.CoordinateSystem.unit_e_1 
material_fn = sympy.Piecewise((MIndex,phi0.sym[0] <= 0.0,),(AIndex, True))


# In[7]:


# from matplotlib.colors import ListedColormap

# var_fn = phi.sym
# var_data = uw.function.evaluate(var_fn, mesh.X.coords)
# vmin = var_data.min()
# vmax = var_data.max()

# radio = -vmin/(vmax-vmin) 
# from matplotlib.colors import LinearSegmentedColormap
# colors = [(0, "blue"), (radio, "white"), (1, "red")]
# cmap = LinearSegmentedColormap.from_list("custom_colormap", colors)


# In[8]:


# import pyvista as pv
# import underworld3.visualisation as vis

# pvmesh = vis.mesh_to_pv_mesh(mesh)
# #pvmesh.point_data["phi"] = vis.scalar_fn_to_pv_points(pvmesh, R0.sym)
# pvmesh.point_data["phi"] = vis.scalar_fn_to_pv_points(pvmesh, visc_fn)
# #pvmesh.point_data["phi"] = vis.scalar_fn_to_pv_points(pvmesh, density_fn*ND_gravity )
# pl = pv.Plotter(window_size=(750, 750))

# pl.add_mesh(
#     pvmesh,
#     cmap='coolwarm',
#     edge_color="Black",
#     show_edges=True,
#     scalars="phi",
#     use_transparency=False,
#     opacity=1,
#     line_width = 0.0
#     )
# pl.show(cpos="xy")


# In[9]:


stokes = uw.systems.Stokes(mesh, velocityField=v, pressureField=p)
stokes.constitutive_model = uw.constitutive_models.ViscousFlowModel
stokes.bodyforce = sympy.Matrix([0, -1 * ND_gravity * density_fn])
stokes.constitutive_model.Parameters.shear_viscosity_0 = visc_fn
stokes.saddle_preconditioner = 1.0 / stokes.constitutive_model.Parameters.shear_viscosity_0
stokes.add_essential_bc((0.0,None), "Left")
stokes.add_essential_bc((0.0,None), "Right")
stokes.add_essential_bc((0.0,0.0), "Bottom")
stokes.add_essential_bc((None,0.0), "Top")

# if uw.mpi.size == 1:
#     stokes.petsc_options['pc_type'] = 'lu'

stokes.tolerance = 1.0e-6
stokes.petsc_options["ksp_rtol"] = 1.0e-6
stokes.petsc_options["ksp_atol"] = 1.0e-6
stokes.petsc_options["snes_converged_reason"] = None
stokes.petsc_options["snes_monitor_short"] = None


# In[10]:


# stokes.solve(zero_init_guess=False)
# mesh.write_timestep("TR",meshUpdates=False,meshVars=[v, p, timeField,R0],outputPath=outputPath,index=0)


# In[11]:


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


# In[12]:


step      = 0
max_steps = 5
time      = 0
dt        = 0

while time < max_time+dt_set:   
    if uw.mpi.rank == 0:
        string = """Step: {0:5d} Model Time: {1:6.1f} dt: {2:6.1f} ({3})\n""".format(
        step, _adjust_time_units(time),
        _adjust_time_units(dt),
        datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        sys.stdout.write(string)
        sys.stdout.flush()

    stokes.solve(zero_init_guess=False)

    if step%save_every ==0:
        if uw.mpi.rank == 0:
            print(f'\nSave data:')
        with mesh.access(timeField):
            timeField.data[:,0] = dim(time, u.megayear).m
        #mesh.petsc_save_checkpoint(meshVars=[v, p, timeField], index=step, outputPath=outputPath)
        mesh.write_timestep("TR",meshUpdates=False,meshVars=[v, p, timeField,phi0],outputPath=outputPath,index=step,)
        interfaceSwarm.write_timestep(filename='swarm',swarmname = 'interfaceSwarm',index = step, outputPath=outputPath)

    dt_solver = stokes.estimate_dt()
    dt = min(dt_solver,dt_set)

    interfaceSwarm.advection(V_fn=stokes.u.sym, delta_t=dt,evalf=False)
    init_phi()

    step += 1
    time += dt


# In[ ]:




