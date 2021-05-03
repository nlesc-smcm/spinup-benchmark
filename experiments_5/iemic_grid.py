from omuse.units import units

import numpy


# note that i-emic defines depth levels as negative numbers!

def depth_levels(N, stretch_factor=1.8): #1.8
    z=numpy.arange(N)/(1.*(N-1))
    if stretch_factor==0:
        return z
    else:
        return 1 - numpy.tanh(stretch_factor*(1-z))/numpy.tanh(stretch_factor)
    
#~ print h

def read_global_mask(Nx,Ny,Nz):
    filename="mask_global_{0}x{1}x{2}".format(Nx, Ny, Nz)
    
    mask=numpy.zeros((Nx+1,Ny+2,Nz+2), dtype='int')
    
    f=open(filename,'r')
    for k in range(Nz+2):
      line=f.readline() # ignored
      for j in range(Ny+2):
          line=f.readline()
          mask[:,j,k]=numpy.array([int(d) for d in line[:-2]]) # ignore last colum + newline
          
    return mask

def depth_array(Nx,Ny,Nz):
    mask=read_global_mask(Nx,Ny,Nz)
                
    depth=mask[:,:,:].argmin(axis=2)
    
    a=depth>0
    depth[a]=Nz-depth[a]
    
    depth=depth[::,::-1]
    
    return depth

if __name__=="__main__":


    d=depth_array(96, 38,12)

    import matplotlib
    matplotlib.use('TkAgg')
    from matplotlib import pyplot
 
 
    pyplot.imshow(d.T,origin="lower", cmap='jet')
    pyplot.show()

#~ print read_global_mask(96,38,12)

