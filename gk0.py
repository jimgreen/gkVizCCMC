import h5py

class File(h5py.File):
    def Time(self):
        return self['timeData'].attrs['vsTime']
    def DomainShape(self):
        """
        Returns:
            An array, say, (Nx,Ny,Nz)
        """
        return self['StructGrid'].attrs['vsNumCells']
    def Dimension(self):
        """
        Returns:
            0, 1, or 2 for 1D, 2D or 3D data
        """
        return len(self.DomainShape())
    def ComponentNumber(self):
        return self['StructGridField'].shape[-1]
    ########################
    # Retrieve coordinates #
    ########################
    def getCoordinates(self, crd, cc=True):
        """
        Args:
          crd: Direction of the coordinate to be retrieved. e.g., for 3d
               data, crd should be either one of the strings 'x', 'y', and
               'z', or one of the integers 0, 1, and 2.
           cc: Return cell-center coordintes. Default is True.
        Returns:
          A 1D arrray. If cc is True, a Nx-element array with cell center
          coordinates is returned; Otherwise, a (Nx+1)-element array with node
          center coordinates is returned, where first element is the lower
          bound, and last element is the upper bound.

        FIXME: The return values are different for uniform and nonuniform grids.
               For uniform grids, a 1D array is returned; for nonuniform grids,
               say, in the case of 3d, a 3d array of shape (Nz,Ny,Nx) is
               returned.
        """
        idx = ('x','y','z').index(crd)
        dim = self.Dimension()
        assert(idx < dim)
        lower = self['StructGrid'].attrs['vsLowerBounds'][idx]
        upper = self['StructGrid'].attrs['vsUpperBounds'][idx]
        num = self['StructGrid'].attrs['vsNumCells'][idx]
        import numpy as np
        crd = np.linspace(lower,upper,num+1,endpoint=True)
        if cc:
            crd = 0.5 * (crd[1:] + crd[:-1])
        return crd
    ########################
    # Retrieve actual data #
    ########################
    def __get1dField(self,fldID,ix=slice(None),dataPath='StructGridField'):
        return self[dataPath][ix,int(fldID)]
    def __get2dField(self,fldID,ix=slice(None),iy=slice(None),\
            dataPath='StructGridField',transpose=True):
        fldID = int(fldID)
        if transpose:
            return self[dataPath][ix,iy,int(fldID)].T
        else:
            return self[dataPath][ix,iy,int(fldID)]
    def __get3dField(self,fldID,ix=slice(None),iy=slice(None),iz=slice(None),\
            dataPath='StructGridField',transpose=True):
        fldID = int(fldID)
        if transpose:
            return self[dataPath][ix,iy,iz,int(fldID)].T
        else:
            return self[dataPath][ix,iy,iz,int(fldID)]
    def getField(self,fldID,**kwargs):
        """Lower level data reading interface.
        Args:
            fldID   : Integer ID number of the field component to be accessed;
                      e.g., 0 for electron mass density for 5-moment data
            **kwargs: **kwargs to be used with lower level functions
                - ix, iy, iz: Can be an integer of a slice object to slice data
                              in x, y, and y dimensions; e.g., ix=slice(None),
                              ix=1, ix=slice(1,-1,2); Defaults are slice(None),
                              i.e., do not slice.
                - dataPath  : Path to which data are stored; default is
                              'StructGridField'.
        Returns:
            A numpy array.
        """
        if self.Dimension() == 1:
            return self.__get1dField(fldID,**kwargs)
        elif self.Dimension() == 2:
            return self.__get2dField(fldID,**kwargs)
        elif self.Dimension() == 3:
            return self.__get3dField(fldID,**kwargs)

