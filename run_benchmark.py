import os
import numpy
import matplotlib
import math
from matplotlib import pyplot
from scipy import integrate

numpy.random.seed(123451)

from omuse.community.pop.interface import POP
from amuse.units import units
from amuse.io import write_set_to_file

from visualisation import sst_plot, crossection, cross_vectorfield


def depth_levels(N, stretch_factor=0):
  z=numpy.arange(N)/(1.*(N-1))
  if stretch_factor==0:
    return z
  else:
    return 1 - numpy.tanh(stretch_factor*(1-z))/numpy.tanh(stretch_factor)


if __name__=="__main__":
#prepare the plot stuff
  pyplot.ion()
  pyplot.show()

#size of the modeling domain, for the benchmark it is equal in all directions
  dim = 16

  p=POP(number_of_workers=2, channel_type="sockets", mode='test')
  
  cwd=os.getcwd()
  
  depth=depth_levels(dim+1, 0.0)
  dz=depth[1:]-depth[:-1]
  
  p.parameters.topography_option='internal'
  p.parameters.horiz_grid_option='amuse'
  p.parameters.vert_grid_option='amuse'     
  p.parameters.vertical_layer_thicknesses=dz * (5000 | units.m)
  #p.parameters.surface_heat_flux_forcing = 'amuse'
  #p.parameters.surface_freshwater_flux_forcing = 'amuse'
  p.parameters.windstress_forcing = 'analytic'
  
# wind forcing 
  tau_x = numpy.zeros((dim,dim))
  tau_y = numpy.zeros((dim,dim))

  forcings = p.forcings.copy()
  channel=forcings.new_channel_to(p.forcings)

  f0 = 8.4e-5
  beta = 1.85e-11
  tau0 = 0.1
  
  for i in range(0,dim):
    for j in range(0,dim):
      cor_par = (f0 + beta*j*4e5)/f0
      tau_x[i][j] = tau0*cor_par*math.sin(2*math.pi*j/dim)

  forcings.tau_x = tau_x | units.Pa
  forcings.tau_y = tau_y | units.Pa
  channel.copy_attributes(["tau_x", "tau_y"])
  
  input()

  print (p.elements.lat.min().in_(units.deg),p.elements.lat.max().in_(units.deg))
  print (p.elements.lon.min().in_(units.deg),p.elements.lon.max().in_(units.deg))
  print ()
  print (p.nodes.lat.min().in_(units.deg),p.nodes.lat.max().in_(units.deg))
  print (p.nodes.lon.min().in_(units.deg),p.nodes.lon.max().in_(units.deg))
  
  x=p.elements.lon.flat
  y=p.elements.lat.flat
  
  input()

  tnow=p.model_time
  dt = 100*365 | units.day
  tend=tnow+(365 * 10000 | units.day)


  while tnow< tend-dt/2:
      print ("evolving to:", tnow+dt)
      p.evolve_model(tnow+dt)
      tnow=p.model_time
      t = tnow.value_in(units.day)
      t = int(t/(10*365))
     
#Preparing for plotting different parameters in SN direction
#zonal averaged 
#coordinates and variables calculated in nodes
      lat = p.nodes3d.lat.value_in(units.deg)
      z = p.nodes3d.z.value_in(units.m) 

      vel_lat_slice = p.nodes3d.yvel[0][:][:]
      vel_z_slice = p.nodes3d.zvel[0][:][:]

      for i in range(1,dim):
        vel_lat_slice = vel_lat_slice + p.nodes3d.yvel[i][:][:]
        vel_z_slice = vel_z_slice + p.nodes3d.zvel[i][:][:]

      vel_lat = vel_lat_slice.value_in(units.cm / units.s)/dim
      vel_z = vel_z_slice.value_in(units.cm / units.s)/dim

#Meridional overturning streamfunction
#based on iemic/matlab/mstream.m
#corrected v_z to v_lat
      y = numpy.linspace(0.0, 5e5, 16)
      z_strf = integrate.cumtrapz(vel_lat, y, initial = 0)
      z_strf = z_strf/1e4

#coordinates and variables calcuated in elements
      lat_el = p.elements3d.lat.value_in(units.deg)
      z_el = p.elements3d.z.value_in(units.m)
      
      density_slice = p.elements3d.rho[0][:][:]
      salinity_slice = p.elements3d.salinity[0][:][:]
      temp_slice = p.elements3d.temperature[0][:][:]

      for i in range(1,dim):
        temp_slice = temp_slice + p.elements3d.temperature[i][:][:]
        density_slice = density_slice + p.elements3d.rho[i][:][:]
        salinity_slice = salinity_slice + p.elements3d.salinity[i][:][:]

      temp_av = temp_slice.value_in(units.C)/dim
      density_av = density_slice.value_in(units.g / units.cm**3)/dim
      salinity = salinity_slice.value_in(units.g / units.kg)/dim
   


#Creating all plots and writing data to files, [] - range for colorbar  
#Plot crossections
      crossection(lat,z,vel_lat,[-6.0,6.0],"y_velocity, cm/s","velocity_lat_"+ \
                  str(t)+".png")
      crossection(lat,z,salinity,[0.34,0.36],"salinity, g/kg","salinity_av_"+  \
                  str(t)+".png")
      crossection(lat_el,z_el,temp_av,[2.0,15.0],"temperature, C", \
                  "temperature_av_"+str(t)+".png")
      crossection(lat_el,z_el,density_av,[0.0,0.04],"density, g/cm3", \
                  "density_av_"+str(t)+".png")
      crossection(lat,z,z_strf,[-6.0,6.0],"meridional overturning streamfunction, Sv", \
                  "merid_streamf"+str(t)+".png")
      cross_vectorfield(lat,z,vel_lat,vel_z,"vel_total"+str(t)+".png")

#Plot surface data
      bar_x = p.nodes.vx_barotropic.value_in(units.cm / units.s ) 
      bar_y = p.nodes.vy_barotropic.value_in(units.cm / units.s ) 
      #sst_plot(p.nodes,tau_y,"psi"+str(t)+".png") 

      tau_x = p.forcings.tau_x.value_in(units.Pa)
      tau_y = p.forcings.tau_y.value_in(units.Pa)

      #sst_plot(p.nodes,tau_x,"tau_x_"+str(t)+".png") 
      #sst_plot(p.nodes,tau_y,"tau_y_"+str(t)+".png") 
      sst_plot(p.nodes,bar_x,"barotropic_x_vel"+str(t)+".png")
      sst_plot(p.nodes,bar_y,"barotropic_y_vel"+str(t)+".png")

