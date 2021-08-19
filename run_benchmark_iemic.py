import logging
import sys
from amuse.test.amusetest import TestWithMPI
from omuse import units
from omuse.community.iemic.interface import iemicInterface
from omuse.community.iemic.interface import iemic
import matplotlib
#~ matplotlib.use("gtk3agg")
from matplotlib import pyplot
import numpy
from visualisation import  crossection, cross_vect_iemic, plot_data
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D

instance = iemic(redirection="none")

instance.Continuation.destination_0=1.
instance.Continuation.backtracking_increase=1.2
instance.Continuation.maximum_number_of_steps=2
instance.Continuation.maximum_step_size=0.05
instance.parameters.Ocean__THCM__Depth_hdim = 5000.0

#instance.ocean__THCM.Levitus_T = 1
#instance.ocean__THCM.Levitus_S = 1
instance.parameters.Ocean__THCM__Global_Bound_ymax = 65.
instance.parameters.Ocean__THCM__Global_Bound_ymin = 10.
instance.parameters.Ocean__THCM__Starting_Parameters__Temperature_Forcing=12.5 
instance.parameters.Ocean__THCM__Starting_Parameters__Salinity_Forcing=0 #.1
instance.parameters.Ocean__THCM__Starting_Parameters__Wind_Forcing=0. #1.0



instance.parameters.Ocean__THCM__Starting_Parameters__SPL1=10.0e3
instance.parameters.Ocean__THCM__Starting_Parameters__SPL2=0.001
#instance.ocean__THCM__Starting_Parameters.P_VC=3.0
#instance.ocean__THCM.Wind_Forcing_Type = 0 #2 for standard conditions like in pop

sets=instance.parameter_set_names()
for name in sets:
  print("parameter set: {0}".format(name))
  print(getattr(instance,name))
  print()
print(instance.grid)
input()
instance.step_continuation()
input()
temp_slice=instance.grid[0,:,:].temperature
sal_slice=instance.grid[0,:,:].salinity
for i in range(1,16):
  temp_slice += instance.grid[i,:,:].temperature
  sal_slice += instance.grid[i,:,:].salinity

temp = temp_slice/16.
sal = sal_slice/16.

density_crossection=(1.035 -  0.0001 * (temp - 12.5) + 0.00076 * (sal-35))*1000

usurf=instance.grid[5,:,:].temperature
plot_data(usurf, "temperature SN single slice", "temperature_SN.png")

usurf=temp
plot_data(usurf, "temperature SN averaged", "temperature_SN_av.png")

usurf=density_crossection
plot_data(usurf, "density SN averaged", "density_SN_av.png")

usurf=sal
plot_data(usurf, "density salinity averaged", "salinity_SN_av.png")

usurf=instance.grid[10,:,:].v_velocity
plot_data(usurf, "v velocity", "v_velocity_SN.png")

vel_lat_slice=instance.grid[0,:,:].v_velocity
vel_z_slice=instance.grid[0,:,:].w_velocity

for i in range(1,15):
  vel_lat_slice = vel_lat_slice + instance.grid[i,:,:].v_velocity
  vel_z_slice = vel_z_slice + instance.grid[i,:,:].w_velocity

vel_lat=vel_lat_slice/16
vel_z=vel_z_slice/16

#Meridional streamfunction
y = numpy.linspace(0.0, 5e3, 16)
y=y[::-1]
cumsum=-vel_lat[:,0]*250.0
print(cumsum)
z_strf=vel_lat.copy()
z_strf[:,0]=cumsum

for i in range(1,len(vel_lat[0,:])):
  cumsum += -vel_lat[:,i]*250
  z_strf[:,i] = cumsum
  
usurf=z_strf
pyplot.imshow(usurf.T, origin="lower")
plot_data(usurf, "streamfunction", "streamfunction_SN.png")


#plot averaged velocity field in vertical-crossection
xs,ys = numpy.mgrid[0:16,0:16]
cross_vect_iemic(xs, ys, vel_lat, vel_z,"averaged_velocity_SN.png")

#plot velocity field in SN-crossection
vel_z=instance.grid[5,:,:].w_velocity
vel_lat=instance.grid[5,:,:].v_velocity
cross_vect_iemic(xs, ys, vel_lat, vel_z,"velocity_SN_5.png")

#plot velocity field in vertical-crossection
vel_lon=instance.grid[:,:,8].u_velocity
vel_lat=instance.grid[:,:,8].v_velocity
cross_vect_iemic(xs, ys, vel_lon, vel_lat,"velocity_top.png")

#3D plot
#xs,ys,zs = numpy.mgrid[0:16,0:16,0:16]
#pyplot.clf()
#fig=pyplot.figure()
#ax=fig.gca(projection="3d")
#ax.quiver(xs,ys,zs,instance.grid[:,:,:].u_velocity,   \
#          instance.grid[:,:,:].v_velocity,            \
#          instance.grid[:,:,:].w_velocity,            \
#          length=0.5, normalize=True)
    #pyplot.gca().invert_yaxis()

#pyplot.show()
#pyplot.clf()

instance.cleanup_code()
instance.stop() 

