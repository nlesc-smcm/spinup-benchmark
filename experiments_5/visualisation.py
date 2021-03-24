import os
import numpy
import matplotlib
import cartopy.crs as ccrs
import cartopy
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
    extent = (-180,180.,-90.,90.)
    xs,ys = numpy.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
    return xs,ys

def sst_plot(grid,sst,name,lowleft=[-80.,-180.] | units.deg,
               upperright=[80,180] | units.deg):
    lats,lons=get_lats_lons(grid)
    xs,ys=get_xy()
    sst=sst.flatten()
    sst_=griddata((lons,lats),sst,(xs,ys), method='linear')
    pyplot.clf()
    ax = pyplot.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    pyplot.pcolormesh(xs,ys,sst_)#,vmin=10,vmax=25)
    pyplot.ylim(ys.min(),ys.max())
    pyplot.colorbar()
    pyplot.savefig(name,format='png')
    pyplot.show()


def crossection(lat, z, parameter, var_range, label, name):

    xs = lat[0][:][:]
    ys = z[0][:][:]

    par_=parameter.transpose()

    pyplot.clf()
    pyplot.pcolormesh(par_[::-1,:]) #,vmin=var_range[0],vmax=var_range[1])
    cbar=pyplot.colorbar()
    #plot contours
    ax=pyplot.contour(par_[::-1,:],colors='k')
    #add labels
    pyplot.clabel(ax,inline=1, fontsize = 10)
    cbar.set_label(label,rotation=90)

    pyplot.savefig(name, format='png')
    
    pyplot.show()
    pyplot.pause(0.5)
    pyplot.clf()

def cross_vectorfield(lat, z, vel_lon, vel_lat, name):
    print(lat.shape)
    xs = lat[:,:,0] 
    ys = z[:,:,0]

    pyplot.clf()
    pyplot.quiver(xs,ys,vel_lon,vel_lat)
    pyplot.savefig(name,format='png')
    
    pyplot.show()
    pyplot.pause(0.5)
    pyplot.clf()

def plot_hor(lat, lon, parameter, var_range, label, name, time):

    xs,ys=get_xy()
    lats= lat[:,:,0].flatten()
    lons= lon[:,:,0].flatten()
    lats+=numpy.random.random(len(lats))*1.e-4
    lons+=numpy.random.random(len(lons))*1.e-4

    a=numpy.where(lons>180)[0]
    lons[a]=lons[a]-360

    par=parameter.flatten()
    par_=griddata((lons,lats),par,(xs,ys), method='linear')

    pyplot.clf()
    #plot coastlines and fields
    ax = pyplot.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    pyplot.pcolormesh(xs,ys,par_)    
    ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
    cbar=pyplot.colorbar()
    #add labels
    cbar.set_label(label,rotation=90)
    pyplot.text(-20.0, -110.0,'time = ' + str(time) + 'years')

    pyplot.savefig(name, format='png')
    
    pyplot.show()
    pyplot.pause(0.5)
    pyplot.clf()
