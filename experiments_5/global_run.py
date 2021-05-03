import os
import numpy
import matplotlib
import math
import sys
from matplotlib import pyplot
from scipy import integrate

numpy.random.seed(123451)
from interface import POP
from omuse.units import units
from amuse.io import write_set_to_file

from visualisation import plot_hor, sst_plot, cross_vectorfield
from iemic_grid import depth_array

from iemic_grid import depth_array,depth_levels

def z_from_cellcenterz(zc, firstlevel=None):
    if firstlevel is None:
      top=0*zc[0]
    z=numpy.zeros(len(zc)+1)*zc[0]
    direction=1
    if zc[0]<=zc[0]*0:
      direction=-1
    for i,_zc in enumerate(zc[::direction]):
      half=_zc-z[i]
      z[i+1]=z[i]+2*half
    return z[::direction]

def barotropic_streamfunction(u, dz, dy):
    """
    calculate barotropic stream function
    
    u: longitude velocity, 3 dim array(lon,lat, z)
    dz: z layer height (possibly array)
    dy: lattitude (physical) cell size (possibly array)
    """
    if len(u.shape)!=3:
        raise Exception("u dim !=3")

    # depth integration
    uint=(u*dz).sum(axis=-1)
    # lattitude integration (note the flip)
    uint=(uint*dy)[:,::-1].cumsum(axis=-1)[:,::-1]

    psib=numpy.zeros((u.shape[0],u.shape[1]))*uint[0,0]
    psib[:,:]=-uint   
    return psib  

def overturning_streamfunction(v, dz, dx):
    """ 
    calculate meriodional overturning streamfunction
    
    v: lattitudinal velocity, 3 dim array (lon, lat, z)
    dz: z layer height (possibly array)
    dx: longitudinal cell size (probably array for lattitude dependend)
    """
    if len(v.shape)!=3:
        raise Exception("v dim !=3")

    #integrate over longitude
    vint=(v.transpose((0,2,1))*dx).transpose((0,2,1))

    vint=vint.sum(axis=0)

    #depth integration
    vint=(vint*dz).cumsum(axis=-1)

    psim=numpy.zeros( (v.shape[1], v.shape[2]))*vint[0,0]
    psim[:,:]=-vint #1:
    return psim

if __name__=="__main__":
#prepare the plot stuff
  pyplot.ion()
  pyplot.show()

  dim = 12

  p=POP(number_of_workers=4, channel_type="sockets", mode='test') 
  

  cwd=os.getcwd()

  d=depth_array(96,118,12)

  d=d[1:,:]

  print (d.shape)
  indices=numpy.indices(d.shape)

  dmask=d  
  depth=depth_levels(dim+1, 1.8)
  dz=depth[1:]-depth[:-1]
  print(dz)
  p.parameters.topography_option='amuse'
  p.set_KMT(indices[0].flatten()+1,indices[1].flatten()+1, d.flat) 
  p.parameters.horiz_grid_option='internal'
  p.parameters.vert_grid_option='amuse'     
  p.parameters.vertical_layer_thicknesses=dz * (5000 | units.m)
  #p.parameters.surface_heat_flux_forcing = 'amuse'
  #p.parameters.surface_freshwater_flux_forcing = 'amuse'
  print(p.parameters)
  print (p.elements)
  print (p.nodes)
  input()
  

  print (p.elements.lat.min().in_(units.deg),p.elements.lat.max().in_(units.deg))
  print (p.elements.lon.min().in_(units.deg),p.elements.lon.max().in_(units.deg))
  input()
  print ()
  print (p.nodes.lat.min().in_(units.deg),p.nodes.lat.max().in_(units.deg))
  print (p.nodes.lon.min().in_(units.deg),p.nodes.lon.max().in_(units.deg))
  
  x=p.elements.lon.flat
  y=p.elements.lat.flat
  

  tnow=p.model_time
  dt = 50 | units.day
  tend = tnow+(365*100 | units.day)
  t = tnow.value_in(units.day)
  t = int(t/(365))
  

  sst = p.elements.temperature.value_in(units.C)
 
  print(sst.shape)
  pyplot.plot(x.value_in(units.deg),y.value_in(units.deg),'r+')
  pyplot.draw()
  sst_plot(p.elements,sst,"temp_init.png")
  fn = 0
  while tnow< tend-dt/2:
      print ("evolving to:", tnow+dt)
      p.evolve_model(tnow+dt)
      tnow=p.model_time
      t = tnow.value_in(units.day)
      t = int(t/(365))
      fn += 1
     
      print("KE", 0.5*(p.nodes3d.xvel**2+p.nodes3d.yvel**2).sum().in_(units.m**2/units.s**2))

#coordinates and variables calculated in nodes
      lat = p.nodes3d.lat.value_in(units.deg)
      lon = p.nodes3d.lon.value_in(units.deg)
      z = p.nodes3d.z.value_in(units.m) 
      
      dy=(p.nodes[0,1].lat-p.nodes[0,0].lat)*(1.|units.Rearth)

      zc=z_from_cellcenterz(p.nodes3d[0,0,:].z)
      dz=zc[1:]-zc[:-1]
      bar_str=barotropic_streamfunction(p.nodes3d.xvel,dz,dy)*1e-6

#coordinates and variables calcuated in elements
      lat_el = p.elements3d.lat.value_in(units.deg)
      lon_el = p.elements3d.lon.value_in(units.deg)
      z_el = p.elements3d.z.value_in(units.m)


#Creating all plots and writing data to files, [] - range for colorbar  
#Plotting
      
      ssh = p.nodes.ssh.value_in(units.m)
      ssh[dmask==0]=ssh.max()
      sst_plot(p.nodes,ssh, "ssh"+str(fn)+".png")

      temp_top=p.elements3d.temperature[:,:,0].value_in(units.C)
      plot_hor(lat_el,lon_el,temp_top,[0.0,15.0],"temperature, C", \
            "temperature_top_"+str(fn)+".png", t)

      plot_hor(lat_el,lon_el,bar_str.value_in(units.m**3/units.s),[-30.0,30.0],"barotropic_streamf", \
            "bar_str_"+str(fn)+".png", t)
#Plot total velocity field
      vel_x_hor = p.nodes3d.xvel[:,:,0].value_in(units.m/units.s)
      vel_y_hor = p.nodes3d.yvel[:,:,0].value_in(units.m/units.s)
      cross_vectorfield(lon,lat,vel_x_hor,vel_y_hor,"vel_total"+str(fn)+".png")
            
