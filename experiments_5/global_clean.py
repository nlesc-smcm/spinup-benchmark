import os
import numpy
import matplotlib
import math
import sys
from matplotlib import pyplot
from scipy import integrate

numpy.random.seed(123451)

from interface import POP
from amuse.units import units
from amuse.io import write_set_to_file

from visualisation import plot_hor, sst_plot, crossection, cross_vectorfield
from iemic_grid import depth_array

from iemic_grid import depth_array,depth_levels

if __name__=="__main__":
#prepare the plot stuff
  pyplot.ion()
  pyplot.show()

  dim = 12

  p=POP(number_of_workers=4, channel_type="sockets", mode='global',redirection='none')
  

  cwd=os.getcwd()

  d=depth_array(96,38,12)

  d=d[1:,:]

  print (d.shape)
  indices=numpy.indices(d.shape)

  dmask=d  
  depth=depth_levels(dim+1, 1.8)
  dz=depth[1:]-depth[:-1]
  print(depth, dz)
  p.parameters.topography_option='amuse'
  p.parameters.vert_grid_option='amuse'     
  p.parameters.vertical_layer_thicknesses=dz * (5000 | units.m)
  #p.parameters.surface_heat_flux_forcing = 'amuse'
  #p.parameters.surface_freshwater_flux_forcing = 'amuse'
  #p.parameters.windstress_monthly_file=os.path.join(cwd,'data/input/ws_monthly/ws.1958-2000.mon')
  print(p.parameters)
  print (p.elements)
  print (p.nodes)
  input()
  p.set_KMT(indices[0].flatten()+1,indices[1].flatten()+1, d.flat)  

  # wind forcing
  tau_x = numpy.zeros((96,40))
  tau_y = numpy.zeros((96,40))

  forcings = p.forcings.copy()
  channel=forcings.new_channel_to(p.forcings)

  f0 = 8.4e-5
  beta = 1.85e-11
  tau0 = 0.1
  
  for i in range(0,96):
    for j in range(0,40):
      cor_par = (f0 + beta*j*4e5)/f0
      tau_x[i][j] = tau0*cor_par*math.sin(2*math.pi*j/dim)

  forcings.tau_x = tau_x | units.Pa
  forcings.tau_y = tau_y | units.Pa
  channel.copy_attributes(["tau_x", "tau_y"])
  

  print (p.elements.lat.min().in_(units.deg),p.elements.lat.max().in_(units.deg))
  print (p.elements.lon.min().in_(units.deg),p.elements.lon.max().in_(units.deg))
  input()
  print (p.nodes.lat.min().in_(units.deg),p.nodes.lat.max().in_(units.deg))
  print (p.nodes.lon.min().in_(units.deg),p.nodes.lon.max().in_(units.deg))
  
  x=p.elements.lon.flat
  y=p.elements.lat.flat
  

  vel_lat_old = p.nodes3d.yvel[0][:][:]
  vel_z_old = p.nodes3d.zvel[0][:][:]
  vel_lat_old = vel_lat_old.value_in(units.cm / units.s)
  vel_z_old = vel_z_old.value_in(units.cm / units.s)
  density_old = p.elements3d.rho[:][:][:]
  density_old = density_old.value_in(units.g / units.cm**3)

  tnow=p.model_time
  dt = 10 | units.day
  tend=tnow+(365*5 | units.day)
  t = tnow.value_in(units.day)
  t = int(t/(365))
  
 
  fn = 0
  while tnow< tend-dt/2:
      print ("evolving to:", tnow+dt)
      p.evolve_model(tnow+dt)
      tnow=p.model_time
      forcings.tau_x = tau_x | units.Pa
      forcings.tau_y = tau_y | units.Pa
      channel.copy_attributes(["tau_x", "tau_y"])
      t = tnow.value_in(units.day)
      t = int(t/(365))
      fn += 1
     
      #Preparing for plotting different parameters in SN direction
      #zonal averaged 
      #coordinates and variables calculated in nodes
      lat = p.nodes3d.lat.value_in(units.deg)
      lon = p.nodes3d.lon.value_in(units.deg)
      z = p.nodes3d.z.value_in(units.m) 
      vel_lat_slice = p.nodes3d.yvel[0,:,:]
      vel_z_slice = p.nodes3d.zvel[0,:,:]

      for i in range(1,96):
        vel_lat_slice = vel_lat_slice + p.nodes3d.yvel[i,:,:]
        vel_z_slice = vel_z_slice + p.nodes3d.zvel[i,:,:]

      vel_lat = vel_lat_slice.value_in(units.m / units.s)/96
      vel_z = vel_z_slice.value_in(units.m / units.s)/96

      #Meridional overturning streamfunction
      cumsum=-vel_lat[:,11]*(5e3/12)
      z_strf=vel_lat.copy()
      z_strf[:,11]=cumsum

      for i in range(len(vel_lat[0,:])-1,0,-1):
        cumsum += vel_lat[:,i]*(5e3/12)   
        z_strf[:,i] = cumsum
      vel_lon_slice=p.nodes3d.xvel[:,:,0]

      for i  in range(1,12):
        vel_lon_slice = vel_lon_slice + p.nodes3d.xvel[:,:,i]
      vel_lon=vel_lon_slice.value_in(units.m / units.s)/12
      cumsum2=-vel_lon[:,0]*(5e3/39)
      bar_str=vel_lon.copy()
      bar_str[:,0]=cumsum2
      
      for i in range(0,len(vel_lon[0,:])): #-1,0,-1):
        cumsum2 += -vel_lon[:,i]*(5e3/39)
        bar_str[:,i] = cumsum2

      #coordinates and variables calcuated in elements
      lat_el = p.elements3d.lat.value_in(units.deg)
      lon_el = p.elements3d.lon.value_in(units.deg)
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
      salinity_cross=p.elements3d.salinity[5][:][:].value_in(units.g / units.kg)

      #Creating all plots and writing data to files, [] - range for colorbar  
      #Plot crossections
      ssh = p.nodes.ssh.value_in(units.m)
      sst = p.elements.temperature.value_in(units.C)
      sst[dmask==0]=sst.max()
      ssh[dmask==0]=ssh.max()
      sst_plot(p.elements,sst, "temp"+str(fn)+".png")
      sst_plot(p.nodes,ssh, "ssh"+str(fn)+".png")
      temp_top=p.elements3d.temperature[:,:,0].value_in(units.C)
      plot_hor(lat_el,lon_el,temp_top,[0.0,15.0],"temperature, C", \
            "temperature_top_"+str(fn)+".png", t)
      plot_hor(lat_el,lon_el,bar_str,[0.0,15.0],"barotropic_streamf", \
            "bar_str_"+str(fn)+".png", t)
      sal_sl = p.elements3d.salinity[2,:,:].value_in(units.g/units.kg)
      crossection(lat,z,sal_sl,[0.34,0.36],"salinity, g/kg","salinity_cross_"+  \
                  str(fn)+".png")
      crossection(lat,z,z_strf,[-6.0,6.0],"meridional overturning streamfunction, Sv", \
                  "merid_streamf"+str(fn)+".png")
      vel_x_hor = p.nodes3d.xvel[:,:,0].value_in(units.m/units.s)
      vel_y_hor = p.nodes3d.yvel[:,:,0].value_in(units.m/units.s)
      cross_vectorfield(lon,lat,vel_x_hor,vel_y_hor,"vel_total"+str(fn)+".png")

            
