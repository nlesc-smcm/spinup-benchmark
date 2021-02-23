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


#This script is for testing whether we are at steady state, 
#mean relative difference density, meridional streamfunction and velocities


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

  dim_h = 32
  dim_v = 16
  p=POP(number_of_workers=4, channel_type="sockets", mode='test')
  
  cwd=os.getcwd()
  
  depth=depth_levels(dim_v+1, 0.0)
  dz=depth[1:]-depth[:-1]
  print(depth, dz) 
  p.parameters.topography_option='internal'
  p.parameters.horiz_grid_option='amuse'
  p.parameters.vert_grid_option='amuse'     
  p.parameters.vertical_layer_thicknesses=dz * (5000 | units.m)
  print(p.parameters.vertical_layer_thicknesses)
  p.parameters.surface_heat_flux_forcing = 'amuse'
  p.parameters.surface_freshwater_flux_forcing = 'amuse'
  
# wind forcing 
  tau_x = numpy.zeros((dim_h,dim_h))
  tau_y = numpy.zeros((dim_h,dim_h))

  forcings = p.forcings.copy()
  channel=forcings.new_channel_to(p.forcings)

  f0 = 8.4e-5
  beta = 1.85e-11
  tau0 = 0.1
  
  for i in range(0,dim_h):
    for j in range(0,dim_h):
      cor_par = (f0 + beta*j*4e5)/f0
      tau_x[i][j] = tau0*cor_par*math.sin(2*math.pi*j/dim_h)

  forcings.tau_x = tau_x | units.Pa
  forcings.tau_y = tau_y | units.Pa
  channel.copy_attributes(["tau_x", "tau_y"])
  
  print (p.elements.lat.min().in_(units.deg),p.elements.lat.max().in_(units.deg))
  print (p.elements.lon.min().in_(units.deg),p.elements.lon.max().in_(units.deg))
  print ()
  print (p.nodes.lat.min().in_(units.deg),p.nodes.lat.max().in_(units.deg))
  print (p.nodes.lon.min().in_(units.deg),p.nodes.lon.max().in_(units.deg))
  
  x=p.elements.lon.flat
  y=p.elements.lat.flat
  
  #Block for calculating residuals
  vel_lat_old = p.nodes3d.yvel[0][:][:]
  vel_z_old = p.nodes3d.zvel[0][:][:]
  vel_lat_old = vel_lat_old.value_in(units.cm / units.s)
  vel_z_old = vel_z_old.value_in(units.cm / units.s)
  density_old = p.elements3d.rho[:][:][:]
  density_old = density_old.value_in(units.kg / units.m**3)
  temp_old = p.elements3d.temperature[:][:][:].value_in(units.C)
  vel_lat = vel_lat_old

  #Streamfunction
  y = numpy.linspace(0.0, 5e3, dim_v)
  cumsum=-vel_lat[:,dim_v-1]*(5e3/dim_v)
  z_strf=vel_lat.copy()
  z_strf[:,dim_v-1]=cumsum
  for i in range(len(vel_lat[0,:])-1,0,-1):
    cumsum += -vel_lat[:,i]*(5e3/dim_v)
    z_strf[:,i] = cumsum
  merid_old = z_strf
  vel_lon_slice=p.nodes3d.xvel[:][:][0]
  for i  in range(1,dim_v):
    vel_lon_slice = vel_lon_slice + p.nodes3d.xvel[:][:][i]
  vel_lon=vel_lon_slice.value_in(units.m / units.s)/dim_v
      
  cumsum2=-vel_lon[dim_h-1,:]*(5e3/dim_h)
  z_strf2=vel_lon.copy()
  z_strf2[dim_h-1,:]=cumsum2
  for i in range(len(vel_lon[:,0])-1,0,-1):
    cumsum2 += -vel_lon[i,:]*(5e3/dim_h)
    z_strf2[i,:] = cumsum2
  bar_old=z_strf2

  out_file1 = open("relative_temp.txt", "w+")
  out_file2 = open("relative_z.txt", "w+")
  out_file3 = open("relative_density.txt", "w+")
  out_file4 = open("relative_merid.txt","w+")  
  out_file5 = open("relative_bar.txt","w+")  
  rel_lat = []
  rel_z = []
  rel_dens = []
  rel_merid = []
  rel_bar = []
  rel_temp = []

  tnow=p.model_time
  dt = 365*10 | units.day
  dt_residual = 48  # hours
  tend=tnow+(365*1500 | units.day)

  j = 0
  k = 0
  start = True

  while tnow< tend-dt/2:
      j += 1
      if j%2 == 0:
         dt = 6 | units.hour
      else:
         dt = 365*10 | units.day
      if start == True:
         start = False
      else: 
         p.evolve_model(tnow+dt)
      tnow=p.model_time
      t = tnow.value_in(units.day)
      t = int(t/(365))
      print("evolving to" ,t)

#Preparing for plotting different parameters in SN direction
#all variables are zonally averaged 
#coordinates and variables calculated in nodes
      lat = p.nodes3d.lat.value_in(units.deg)
      lon = p.nodes3d.lon.value_in(units.deg)
      z = p.nodes3d.z.value_in(units.m) 
      vel_lat_slice = p.nodes3d.yvel[0][:][:]
      vel_z_slice = p.nodes3d.zvel[0][:][:]

      for i in range(1,dim_h):
        vel_lat_slice = vel_lat_slice + p.nodes3d.yvel[i][:][:]
        vel_z_slice = vel_z_slice + p.nodes3d.zvel[i][:][:]

  
      vel_lat = vel_lat_slice.value_in(units.m / units.s)/dim_h
      vel_z = vel_z_slice.value_in(units.m / units.s)/dim_h
#Meridional overturning streamfunction
      y = numpy.linspace(0.0, 5e3, dim_v)
      cumsum=-vel_lat[:,dim_v-1]*(5e3/dim_v)
      z_strf=vel_lat.copy()
      z_strf[:,dim_v-1]=cumsum
      for i in range(len(vel_lat[0,:])-1,0,-1):
        cumsum += -vel_lat[:,i]*(5e3/dim_v)
        z_strf[:,i] = cumsum

      vel_lon_slice=p.nodes3d.xvel[:][:][0]
      for i  in range(1,dim_h):
        vel_lon_slice = vel_lon_slice + p.nodes3d.xvel[:][:][i]
      vel_lon=vel_lon_slice.value_in(units.m / units.s)/dim_h
      
      cumsum2=-vel_lon[dim_h-1,:]*(5e3/dim_h)
      z_strf2=vel_lon.copy()
      z_strf2[dim_h-1,:]=cumsum2
      for i in range(len(vel_lon[:,0])-1,0,-1):
        cumsum2 += -vel_lon[i,:]*(5e3/dim_h)
        z_strf2[i,:] = cumsum2


#coordinates and variables calcuated in elements
      lat_el = p.elements3d.lat.value_in(units.deg)
      lon_el = p.elements3d.lon.value_in(units.deg)
      z_el = p.elements3d.z.value_in(units.m)
      
      density_slice = p.elements3d.rho[0][:][:]
      salinity_slice = p.elements3d.salinity[0][:][:]
      temp_slice = p.elements3d.temperature[0][:][:]

      for i in range(1,dim_h):
        temp_slice = temp_slice + p.elements3d.temperature[i][:][:]
        density_slice = density_slice + p.elements3d.rho[i][:][:]
        salinity_slice = salinity_slice + p.elements3d.salinity[i][:][:]

      temp_av = temp_slice.value_in(units.C)/dim_h
      density_av = density_slice.value_in(units.kg / units.m**3)/dim_h
      salinity = salinity_slice.value_in(units.g / units.kg)/dim_h
      salinity_cross=p.elements3d.salinity[5][:][:].value_in(units.g / units.kg)
#Creating all plots and writing data to files, [] - range for colorbar  
#Plot crossections
      if j % 2 != 0:
         temp_top=p.elements3d.temperature[:,:,0].value_in(units.C)
         temp_bot=p.elements3d.temperature[:,:,15].value_in(units.C)

         crossection(lat_el,z_el,vel_lat,[-6.0,6.0],"y_velocity, m/s","velocity_lat_"+ \
                     str(t)+".png")
         crossection(lat,z,vel_lat,[-6.0,6.0],"z_velocity, m/s","velocity_z_"+ \
                     str(t)+".png")
         crossection(lat,z,salinity,[0.34,0.36],"salinity, g/kg","salinity_av_"+  \
                     str(t)+".png")
         temp_cr=p.elements3d.temperature[10][:][:].value_in(units.C)
         crossection(lat_el,z_el,temp_av,[0.0,15.0],"temperature, C", \
                  "temperature_av_"+str(t)+".png")
         crossection(lat_el,z_el,density_av,[0.0,0.04],"density, kg/m3", \
                     "density_av_"+str(t)+".png")
         crossection(lat,z,z_strf,[-6.0,6.0],"meridional overturning streamfunction, Sv", \
                     "merid_streamf"+str(t)+".png")  

#Check residuals
      density = p.elements3d.rho[:][:][:].value_in(units.kg/units.m**3)
      temperature = p.elements3d.temperature[:][:][:].value_in(units.C)
      vel_lat_diff=vel_lat - vel_lat_old
      vel_z_diff = vel_z - vel_z_old
      density_diff = density - density_old
      vel_lat_old = vel_lat
      vel_z_old = vel_z
      temp_diff = temperature - temp_old
      temp_old = temperature 
      density_old = density
      merid_diff = z_strf - merid_old
      bar_diff = z_strf2 - bar_old
      merid_old = z_strf
      bar_old = z_strf2

      if j % 2 != 0:
         k += 1
         rel_lat.append(numpy.mean(abs(vel_lat_diff))/dt_residual)
         rel_z.append(numpy.mean(abs(vel_z_diff))/dt_residual)
         rel_dens.append(numpy.mean(abs(density_diff))/dt_residual)
         rel_merid.append(numpy.mean(abs(merid_diff))/dt_residual)
         rel_bar.append(numpy.mean(abs(bar_diff))/dt_residual)
         rel_temp.append(numpy.mean(abs(temp_diff))/dt_residual)

  for l in range(0,k-1):
      out_file1.write(str(rel_temp[l]) + '\n')
      out_file2.write(str(rel_z[l]) + '\n')
      out_file3.write(str(rel_dens[l]) + '\n')
      out_file4.write(str(rel_merid[l]) + '\n')
      out_file5.write(str(rel_bar[l])+ '\n')
      
  bar_x = p.nodes.vx_barotropic.value_in(units.cm / units.s ) 
  bar_y = p.nodes.vy_barotropic.value_in(units.cm / units.s ) 

#Write for restart
  ssh = p.nodes.ssh.value_in(units.cm)
  gradx = p.nodes.gradx.value_in(units.cm**2)
  grady = p.nodes.grady.value_in(units.cm**2)
  out_file1 = open("temperature.txt","w+")
  out_file2 = open("salinity.txt","w+")
  out_file3 = open("density.txt","w+")
  out_file4 = open("vel_x.txt","w+")
  out_file5 = open("vel_y.txt","w+")
  out_file6 = open("bar_x.txt","w+")
  out_file7 = open("bar_y.txt","w+")
  out_file8 = open("ssh.txt","w+")
  out_file9 = open("gradx.txt","w+")
  out_file10 = open("grady.txt","w+")
  for i in range(0,dim_h):
    for j in range(0,dim_h):
      for k in range(0,dim_v):
        out_file1.write(str(p.elements3d.temperature[i][j][k].value_in(units.C))  + '\n')
        out_file2.write(str(p.elements3d.salinity[i][j][k].value_in(units.g/units.kg))  + '\n')
        out_file3.write(str(p.elements3d.rho[i][j][k].value_in(units.kg/units.m**3))  + '\n')
        out_file4.write(str(p.nodes3d.xvel[i][j][k].value_in(units.cm/units.s))  + '\n')
        out_file5.write(str(p.nodes3d.yvel[i][j][k].value_in(units.cm/units.s))  + '\n')
  
  for i in range(0,dim_h):
    for j in range(0,dim_h):
      out_file6.write(str(bar_x[i][j]) + '\n')
      out_file7.write(str(bar_y[i][j]) + '\n')
      out_file8.write(str(ssh[i][j]) + '\n')
      out_file9.write(str(gradx[i][j]) + '\n')
      out_file10.write(str(grady[i][j]) + '\n')


