import os
import numpy
import matplotlib
import cartopy.crs as ccrs
import sys

from matplotlib import pyplot
from scipy.interpolate import griddata

from amuse.units import units


def get_lats_lons(nodes):
    lats= nodes.lat.value_in(units.deg).flatten() 
    lons= nodes.lon.value_in(units.deg).flatten() 
    lats+=numpy.random.random(len(lats))*1.e-4
    lons+=numpy.random.random(len(lons))*1.e-4    

    a=numpy.where(lons>180)[0]
    lons[a]=lons[a]-360
    return lats,lons

def get_xy():
    N = 500j
    extent = (0,60.,0.,60.)
    xs,ys = numpy.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
    return xs,ys

def sst_plot(grid,sst,name):
    lats,lons=get_lats_lons(grid)
    xs,ys=get_xy()
    sst=sst.flatten()
    sst_=griddata((lons,lats),sst,(xs,ys), method='linear')
    
    pyplot.clf()
    pyplot.pcolormesh(xs,ys,sst_)
    pyplot.colorbar()
    pyplot.savefig(name,format='png')
   
    pyplot.show()
    pyplot.pause(0.5)
    pyplot.clf()

def crossection(lat, z, parameter, var_range, name):

# structure of arrays is [lon][lat][depth]
    xs = lat[0][:][:]
    ys = z[0][:][:]
#fast solution for plotting to rotate data(have depth in z
#next should invert z as surface now at the bottom(with [::-1] 

    par_=parameter.transpose()

    pyplot.clf()
    pyplot.pcolormesh(par_[::-1,:],vmin=var_range[0],vmax=var_range[1])
    pyplot.colorbar()
#plot contours
    ax=pyplot.contour(par_[::-1,:],colors='k')
#add labels
    pyplot.clabel(ax,inline=1, fontsize = 10)
    cbar.set_label(label,rotation=90)

    pyplot.savefig(name, format='png')
    
    pyplot.show()
    pyplot.pause(0.5)
    pyplot.clf()

def cross_vectorfield(lat, z, vel_lat, vel_z, name):
    print(lat.shape)
    xs = lat[0][:][:] 
    ys = z[0][:][:]

    pyplot.clf()
    pyplot.quiver(xs,ys,vel_lat,vel_z*200)
    pyplot.gca().invert_yaxis()
    pyplot.savefig(name,format='png')
    
    pyplot.show()
    pyplot.pause(0.5)
    pyplot.clf()

