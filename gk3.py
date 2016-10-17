import h5py

base_hydro_fields = [ # directly accessible physical quantities
    'rho','rhoux','rhouy','rhouz','e',]
extended_hydro_fields = [ # derivable physical quantities
    'ux','uy','uz','Pdynamic','P',
    ]

base_5m_fields = [ # directly accessible physical quantities
    'rho_e','rhoux_e','rhouy_e','rhouz_e','u_e',
    'rho_i','rhoux_i','rhouy_i','rhouz_i','u_i',
    'Ex','Ey','Ez','Bx','By','Bz']
extended_5m_fields = [ # derivable physical quantities
    'n_e','ux_e','uy_e','uz_e','Jx_e','Jy_e','Jz_e','Pdynamic_e','P_e','beta_e',
    'n_i','ux_i','uy_i','uz_i','Jx_i','Jy_i','Jz_i','Pdynamic_i','P_i','beta_i',
    'rho','ux','uy','uz','Jx','Jy','Jz','Pdynamic','P','beta',
    'Pmagnetic',
    ]

base_10m_fields = [ # directly accessible physical quantities
    'rho_e','rhoux_e','rhouy_e','rhouz_e',
    'pxx_e','pxy_e','pxz_e','pyy_e','pyz_e','pzz_e',
    'rho_i','rhoux_i','rhouy_i','rhouz_i',
    'pxx_i','pxy_i','pxz_i','pyy_i','pyz_i','pzz_i',
    'Ex','Ey','Ez','Bx','By','Bz']
extended_10m_fields = [ # derivable physical quantities
    'n_e','ux_e','uy_e','uz_e','Jx_e','Jy_e','Jz_e','Pdynamic_e',
    'Pxx_e','Pxy_e','Pxz_e','Pyy_e','Pyz_e','Pzz_e',
    'P_e','Ppar_e','Pperp_e','beta_e','betaPar_e','betaPerp_e','anisotropy_e',
    'Txx_e','Txy_e','Txz_e','Tyy_e','Tyz_e','Tzz_e',
    'T_e','Tpar_e','Tperp_e','u_e',
    'n_i','ux_i','uy_i','uz_i','Jx_i','Jy_i','Jz_i','Pdynamic_i',
    'Pxx_i','Pxy_i','Pxz_i','Pyy_i','Pyz_i','Pzz_i',
    'P_i','Ppar_i','Pperp_i','beta_i','betaPar_i','betaPerp_i','anisotropy_i',
    'Txx_i','Txy_i','Txz_i','Tyy_i','Tyz_i','Tzz_i',
    'T_i','Tpar_i','Tperp_i','u_i',
    'rho','ux','uy','uz','Jx','Jy','Jz','Pdynamic','P','beta',
    'Pmagnetic',
    ]

field_component_numbers = [5, 18, 28]
field_types = ['hydro', 'five-moment', 'ten-moment']
base_fields = [base_hydro_fields, base_5m_fields, base_10m_fields]
extended_fields = [extended_hydro_fields, extended_5m_fields, extended_10m_fields]

# LaTex names for different physical quantities
base_10m_field_names = [
    '$\\rho_e}$','$\\rho u_{x,e}$','$\\rho u_{y,e}$','$\\rho u_{z,e}$',
    '$\\mathcal{P}_{xx,e}$','$\\mathcal{P}_{xy,e}$','$\\mathcal{P}_{xz,e}$',
    '$\\mathcal{P}_{yy,e}$','$\\mathcal{P}_{yz,e}$','$\\mathcal{P}_{zz,e}$',
    '$\\rho_i}$','$\\rho u_{x,i}$','$\\rho u_{y,i}$','$\\rho u_{z,i}$',
    '$\\mathcal{P}_{xx,i}$','$\\mathcal{P}_{xy,i}$','$\\mathcal{P}_{xz,i}$',
    '$\\mathcal{P}_{yy,i}$','$\\mathcal{P}_{yz,i}$','$\\mathcal{P}_{zz,i}$',
    '$E_x$','$E_y$','$E_z$','$B_x$','$B_y$','$B_z$'
]

extended_10m_field_names = [
    '$n_e$','$u_{x,e}$','$u_{y,e}$','$u_{z,e}$',
    '$J_{x,e}$','$J_y,e}$','$J_{z,e}$','$\\rho_eu_e^2/2$',
    '$P_{xx,e}$','$P_{xy,e}$','$P_{xz,e}$',
    '$P_{yy,e}$','$P_{yz,e}$','$P_{zz,e}$',
    '$p_e$','$p_{\\parallel,e}$','$p_{\\perp,e}$',
    '$\\beta_e$','$\\beta_{\\parallel,e}$','$\\beta_{\\perp,e}$',
    '$p_{\\parallel,e}/p_{\\perp,e}$',
    '$T_{xx,e}$','$T_{xy,e}$','$T_{xz,e}$',
    '$T_{yy,e}$','$T_{yz,e}$','$T_{zz,e}$',
    '$T_e$','$T_{\\parallel,e}$','$T_{\\perp,e}$','$\epsilon_e$',
    '$n_i$','$u_{x,i}$','$u_{y,i}$','$u_{z,i}$',
    '$J_{x,i}$','$J_{y,i}$','$J_{z,i}$','$\\rho_iu_i^2/2$',
    '$P_{xx,i}$','$P_{xy,i}$','$P_{xz,i}$',
    '$P_{yy,i}$','$P_{yz,i}$','$P_{zz,i}$',
    '$p_i$','$p_{\\parallel,i}$','$p_{\\perp,i}$',
    '$\\beta_i$','$\\beta_{\\parallel,i}$','$\\beta_{\\perp,i}$',
    '$p_{\\parallel,i}/p_{\\perp,i}$',
    '$T_{xx,i}$','$T_{xy,i}$','$T_{xz,i}$',
    '$T_{yy,i}$','$T_{yz,i}$','$T_{zz,i}$',
    '$T_i$','$T_{\\parallel,i}$','$T_{\\perp,i}$','$\epsilon_i$',
    '$\\rho$','$u_x$','$u_y$','$u_z$','$J_x$','$J_y$','$J_z$',
    '$\\rho u^2/2$','$P$','$\\beta$','$B^2/2\\mu_0$',
]

def get_field_name(fld):
    """Gets the LaTex code of the name of the physical quantity
    """
    if fld in base_10m_fields:
        return base_10m_field_names[base_10m_fields.index(fld)]
    elif fld in extended_10m_fields:
        return extended_10m_field_names[extended_10m_fields.index(fld)]
    else:
        return str(fld)

class File(h5py.File):
    """Representation of a Gkeyll data file extended from h5py.File.
    """
    def __init__(self, filename, **kwargs):
        """Initialize a Gkeyll file object.
        Args:
            **kwargs: gamma, me, mi, qe, qi, etc., to be stored in params (a
                dict)
        """
        super(File, self).__init__(filename, 'r')
        self.params = {'gamma':5./3., 'transpose':True}
        self.params.update(**kwargs)
    #############################
    # Retrieve data information #
    #############################
    def Time(self):
        return self['timeData'].attrs['vsTime']
    def HasGhost(self):
        if (self.GridType() == 'uniform'):
            vsNumCells = self['StructGrid'].attrs['vsNumCells']
            fieldShape = self['StructGridField'].shape[:-1]
            # print('real cell: ', vsNumCells)
            # print('available cell: ', fieldShape)
            for i in range(self.Dimension()):
                if vsNumCells[i] != fieldShape[i]:
                    return True
            return False
        else:
            raise RuntimeError('Only support grid type uniform')
    def DomainShape(self):
        try:
            return self['StructGrid'].attrs['vsNumCells']
        except:
            return self['StructGridField'].shape[:-1]
    def Dimension(self):
        return len(self.DomainShape())
    def ComponentNumber(self):
        return self['StructGridField'].shape[-1]
    def FieldType(self):
        return field_types[field_component_numbers.index(self.ComponentNumber())]
    def Fields(self):
        return base_fields[field_types.index(self.FieldType())]
    def ExtendedFields(self):
        return extended_fields[field_types.index(self.FieldType())]
    def GridType(self):
        """
        Returns:
            A string, either 'uniform' or 'structured'
        """
        return self['StructGrid'].attrs['vsKind']
    def LowerBounds(self):
        """
        Returns:
            An array containing lower bounds in each direction. For example, for
            3D, (xlo, ylo, zlo) will be returned. Works for uniform grids only.
        """
        if self.GridType() == 'uniform':
            return self['StructGrid'].attrs['vsLowerBounds']
        elif self.GridType() == 'structured':
            if self.Dimension() == 1:
                return \
                  self['StructGrid'][0,0]
            if self.Dimension() == 2:
                return \
                  self['StructGrid'][0,0,0],\
                  self['StructGrid'][0,0,1]
            if self.Dimension() == 3:
                return \
                  self['StructGrid'][0,0,0,0],\
                  self['StructGrid'][0,0,0,1],\
                  self['StructGrid'][0,0,0,2]
    def UpperBounds(self):
        """
        Returns:
            An array containing upper bounds in each direction. For example, for
            3D, (xhi, yhi, zhi) will be returned. Works for uniform grids only.
        """
        if self.GridType() == 'uniform':
            return self['StructGrid'].attrs['vsUpperBounds']
        elif self.GridType() == 'structured':
            if self.Dimension() == 1:
                return \
                  self['StructGrid'][-1,0]
            if self.Dimension() == 2:
                return \
                  self['StructGrid'][-1,0,0],\
                  self['StructGrid'][0,-1,1]
            if self.Dimension() == 3:
                return \
                  self['StructGrid'][-1,0,0,0],\
                  self['StructGrid'][0,-1,0,1],\
                  self['StructGrid'][0,0,-1,2]
    def addParam(self, name, value):
        self.params[name] = value
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
        gridType = self['StructGrid'].attrs['vsKind']
        dim = self.Dimension()
        assert(idx < dim)
        if gridType in ['uniform']:
            lower = self['StructGrid'].attrs['vsLowerBounds'][idx]
            upper = self['StructGrid'].attrs['vsUpperBounds'][idx]
            num = self['StructGrid'].attrs['vsNumCells'][idx]
            import numpy as np
            crd = np.linspace(lower,upper,num+1,endpoint=True)
        elif gridType in ['structured']:
            crd = self.getField(idx, dataPath='StructGrid', transpose=True)
            if dim == 1:
                if idx == 0:
                    crd = crd[:]
            elif dim == 2:
                if idx == 0:
                    crd = crd[0,:]
                elif idx == 1:
                    crd = crd[:,0]
            elif dim == 3:
                if idx == 0:
                    crd = crd[0,0,:]
                elif idx == 1:
                    crd = crd[0,:,0]
                elif idx == 2:
                    crd = crd[:,0,0]
        if cc:
            crd = 0.5 * (crd[1:] + crd[:-1])
        return crd
    def getCoordinates_Dense(self, crd, cc=True, **kwargs):
        assert(self.GridType() == 'structured')
        idx = ('x','y','z').index(crd)
        dim = self.Dimension()
        assert(idx < dim)
        crd = self.getField(idx, dataPath='StructGrid', **kwargs)
        if cc:
            if dim == 1:
                if idx == 0:
                    crd = 0.5 * (crd[1:] + crd[:-1])
            elif dim == 2:
                if idx == 0:
                    crd = 0.5 * (crd[1:,1:] + crd[1:,:-1])
                elif idx == 1:
                    crd = 0.5 * (crd[1:,1:] + crd[:-1,1:])
            elif dim == 3:
                if idx == 0:
                    crd = 0.5 * (crd[1:,1:,1:] + crd[1:,1:,:-1])
                elif idx == 1:
                    crd = 0.5 * (crd[1:,1:,1:] + crd[1:,:-1,1:])
                elif idx == 2:
                    crd = 0.5 * (crd[1:,1:,1:] + crd[:-1,1:,1:])
        return crd
    ########################
    # Retrieve actual data #
    ########################
    def nearestIndex(self, val, dir):
        """
        Args:
            val: Value of coordinate.
            dir: Direction of calculation.
        Returns:
            The index where coordinate is closest to val.
        """
        assert(self.GridType() == 'uniform')
        val = float(val)
        crdIdx = ('x','y','z').index(dir)
        lo = float(self.LowerBounds()[crdIdx])
        hi = float(self.UpperBounds()[crdIdx])
        if val<lo or val>hi:
           raise ValueError('%g not in range[%g, %g]'%(val,lo,hi))
        N = self.DomainShape()[crdIdx]
        return int(N*(val-lo)/(hi-lo))
    def calcIndex(self, idx, dir):
        """
        Args:
            idx: Can be one of the following:
                 - An integer or a slice object
                 - None
                 - A string whose first two chars give an operator ('=:', '>:',
                   or '<:') and the remaining chars give the coordinate as a
                   number.
            dir: Direction of calculation.
        Returns:
            An integer of a slice object.
            - If idx is an integer or a slice object, the idx itself is
              returned.
            - If idx is None, slice(None) is returned.
            - If idx is a string, a slice object that approximate the
              corresponding relation is returned. For example, if idx='<:val',
              then the function returns slice(None, thisIndex, None), where
              thisIndex is the index where coordinate is closest to val.
        """
        if type(idx) in [slice, int]:
            return idx
        elif idx in [None, slice(None)]:
            return slice(None)
        elif type(idx) in [str]:
            op = idx[0:2]
            thisIdx = self.nearestIndex(float(idx[2:]), dir)
            if op == '=:':
                return thisIdx
            elif op == '<:':
                return slice(None, thisIdx, None)
            elif op == '>:':
                return slice(thisIdx, None, None)
            else:
                raise ValueError('op %s not recognized'%op)
    def __get1dField(self,fldID,ix=slice(None),dataPath='StructGridField'):
        ix = self.calcIndex(ix, 'x')
        fldID = int(fldID)
        return self[dataPath][ix,fldID]
    def __get2dField(self,fldID,ix=slice(None),iy=slice(None),\
            dataPath='StructGridField',transpose=None):
        ix = self.calcIndex(ix, 'x')
        iy = self.calcIndex(iy, 'y')
        fldID = int(fldID)
        if transpose == None:
            transpose = self.params['transpose']
        if transpose:
            return self[dataPath][ix,iy,fldID].T
        else:
            return self[dataPath][ix,iy,fldID]
    def __get3dField(self,fldID,ix=slice(None),iy=slice(None),iz=slice(None),\
            dataPath='StructGridField',transpose=None):
        ix = self.calcIndex(ix, 'x')
        iy = self.calcIndex(iy, 'y')
        iz = self.calcIndex(iz, 'z')
        fldID = int(fldID)
        if transpose == None:
            transpose = self.params['transpose']
        if transpose:
            return self[dataPath][ix,iy,iz,fldID].T
        else:
            return self[dataPath][ix,iy,iz,fldID]
    def __getField(self,fldID,**kwargs):
        """Lower level data reading interface.
        Args:
            fldID   : Integer ID number of the field component to be accessed;
                      e.g., 0 for electron mass density for 5-moment data
            **kwargs: **kwargs to be used with lower level functions
                - ix, iy, iz: Can be an integer of a slice object to slice data
                              in x, y, and y dimensions; e.g., ix=slice(None),
                              ix=1, ix=slice(1,-1,2); Defaults are slice(None),
                              i.e., do not slice.

                              Can also be a string whose first two chars give an
                              operator ('=:', '>:', or '<:') and the remaining
                              chars give the coordinate as a number. The slice
                              object will approximate the corresponding relation
                              is returned. For example, if idx='<:val', then the
                              function returns slice(None, thisIndex, None),
                              where thisIndex is the index where coordinate is
                              closest to val.

                              For 1d data, ix is optional.
                              For 2d data, ix and iy are optional.
                              For 3d data, ix, iy, and iz are optional.
                - dataPath  : Path to which data are stored; default is
                              'StructGridField'.
                - transpose : For 2d/3d data, transpose the data array or not;
                              default is True, i.e., transpose the data.
        Returns:
            A numpy array.
        """
        fldID = int(fldID)
        if self.Dimension() == 1:
            return self.__get1dField(fldID,**kwargs)
        elif self.Dimension() == 2:
            return self.__get2dField(fldID,**kwargs)
        elif self.Dimension() == 3:
            return self.__get3dField(fldID,**kwargs)
    def getField(self,fld,**kwargs):
        """High level data reading interface.
        
        Reading data of the physical quantity named fld (if fld is a string) or
        with the sequential id number fld (if fld is an integer).

        Call Fields() to get list of directly accessible physical quantities;
        Call ExtendedFields() to get list of physical quantities that can be
        derived from directly accessible physical quantities.
  
        For five-moment data, the string name fld can be either
        1) One of the following directly stored quantities:
          # electron moments:
            'rho_e','rhoux_e','rhouy_e','rhouz_e','u_e'
          # ion moments:
            'rho_i','rhoux_i','rhouy_i','rhouz_i','u_i'
          # electromagnetic field components:
            'Ex','Ey','Ez','Bx','By','Bz'
        Notes:
        e = p/(gamma-1) + rho * v**2 / 2
      
        or
        
        2)One of the following derived quantities:
          # electron moments:
            'n_e','ux_e','uy_e','uz_e','Jx_e','Jy_e','Jz_e',
            'Pdynamic_i','P_e','beta_e'
          # ion moments:
            'n_i','ux_i','uy_i','uz_i','Jx_i','Jy_i','Jz_i',
            'Pdynamic_i','P_i','beta_i'
          # total moments:
            'rho','ux','uy','uz','Jx','Jy','Jz','Pdynamic','P','beta'
          # electromagnetic quantities
            'Pmagnetic','E2','B2'

        Args:
            fld     : Name or sequential id number of the physical quantity. If
                      fld is a string, it should be one of the supported field
                      names.
            **kwargs: **kwargs to be used with lower level function __getField.
                      See the documentation for__getField for its options.
        Returns:
            A numpy array.
        """
        #print('getField(%s)'%fld)
        try:
            fldID = int(fld)
            return self.__getField(fldID,**kwargs)
        except:
            if fld in self.Fields():
                return self.__getField(self.Fields().index(fld),**kwargs)
            elif fld in self.ExtendedFields():
                return self.__calcField(fld,**kwargs)
            else:
                return self.eval(fld,**kwargs)
    def __calcField(self,fld,**kwargs):
        def eval__(fld):
            return self.eval(fld, **kwargs)
        import numpy as np
        if fld in ['n_e']:
            return eval__('rho_e/me')
        elif fld in ['ux_e']:
            return eval__('rhoux_e/rho_e')
        elif fld in ['uy_e']:
            return eval__('rhouy_e/rho_e')
        elif fld in ['uz_e']:
            return eval__('rhouz_e/rho_e')
        elif fld in ['Jx_e']:
            return eval__('rhoux_e*qe/me')
        elif fld in ['Jy_e']:
            return eval__('rhouy_e*qe/me')
        elif fld in ['Jz_e']:
            return eval__('rhouz_e*qe/me')
        elif fld in ['Pdynamic_e']:
            return eval__('0.5*(rhoux_e**2+rhouy_e**2+rhouz_e**2)/rho_e')
        elif fld in ['P_e']: # e = P/(gamma-1) + rho * v**2 / 2
            if self.FieldType() in ['five-moment']:
                return eval__('(gamma-1)*(u_e - Pdynamic_e)')
        elif fld in ['beta_e']:
            return eval__('P_e/Pmagnetic')
        if fld in ['n_i']:
            return eval__('rho_i/mi')
        elif fld in ['ux_i']:
            return eval__('rhoux_i/rho_i')
        elif fld in ['uy_i']:
            return eval__('rhouy_i/rho_i')
        elif fld in ['uz_i']:
            return eval__('rhouz_i/rho_i')
        elif fld in ['Jx_i']:
            return eval__('rhoux_i*qi/mi')
        elif fld in ['Jy_i']:           
            return eval__('rhouy_i*qi/mi')
        elif fld in ['Jz_i']:           
            return eval__('rhouz_i*qi/mi')
        elif fld in ['Pdynamic_i']:
            return eval__('0.5*(rhoux_i**2+rhouy_i**2+rhouz_i**2)/rho_i')
        elif fld in ['P_i']: # e = P/(gamma-1) + rho * v**2 / 2
            if self.FieldType() in ['five-moment']:
                return eval__('(gamma-1)*(u_i - Pdynamic_i)')
        elif fld in ['beta_i']:
            return eval__('P_i/Pmagnetic')
        elif fld in ['rho']:
            return eval__('rho_e+rho_i')
        elif fld in ['ux']:
            if self.FieldType() in ['five-moment','ten-moment']:
                return eval__('(rhoux_e+rhoux_i)/(rho_e+rho_i)')
            elif self.FieldType() in ['hydro']:
                return eval__('rhoux/rho')
        elif fld in ['uy']:
            if self.FieldType() in ['five-moment','ten-moment']:
                return eval__('(rhouy_e+rhouy_i)/(rho_e+rho_i)')
            elif self.FieldType() in ['hydro']:
                return eval__('rhouy/rho')
        elif fld in ['uz']:
            if self.FieldType() in ['five-moment','ten-moment']:
                return eval__('(rhouz_e+rhouz_i)/(rho_e+rho_i)')
            elif self.FieldType() in ['hydro']:
                return eval__('rhouz/rho')
        elif fld in ['Jx']:
            return eval__('Jx_e+Jx_i')
        elif fld in ['Jy']:
            return eval__('Jy_e+Jy_i')
        elif fld in ['Jz']:
            return eval__('Jz_e+Jz_i')
        elif fld in ['Pdynamic']:
            if self.FieldType() in ['five-moment','ten-moment']:
                return eval__('Pdynamic_e+Pdynamic_i')
            elif self.FieldType() in ['hydro']:
                return eval__('0.5*(rhoux**2+rhouy**2+rhouz**2)/rho')
        elif fld in ['P']:
            if self.FieldType() in ['five-moment','ten-moment']:
                return eval__('P_e+P_i')
            elif self.FieldType() in ['hydro']:
                return eval__('(gamma-1)*(e - Pdynamic)')
        elif fld in ['beta']:
            return eval__('P/Pmagnetic')
        elif fld in ['Pmagnetic']:
            return eval__('(0.5/mu0)*(Bx**2+By**2+Bz**2)')
        else:
            raise NameError('%s is not implemented in __calcField' % fld)
    def eval(self, eqn, **kwargs):
        """
        Args:
            eqn     : A base field term or any combination of base terms.
            **kwargs: Keyword arguments to be passed to getField.
        Returns:
            A numpy array.

        FIXME: Possibly redundant getter, e.g., when evaulating Jz_e+Jz, Jz_e is
               accessed twice!
        """
        import numexpr as ne
        import re
        import random
        import string
        salt = ''.join(random.SystemRandom().choice(string.ascii_uppercase
            + string.ascii_lowercase) for _ in range(16))
        _symbol_re = r'\b([_A-Za-z][_a-zA-Z0-9]*)\b'
        var_re = _symbol_re + r'(?!\s*\()'
        flds = []
        # for security
        eqn = eqn.replace('__', '')
        local_dict = dict()
        
        def getter(fld):
            if fld in self.params.keys():
                return self.params[fld]
            else:
                return self.getField(fld, **kwargs)
        def var_salter(symbols):
            symbol = symbols.groups()[0]
            if (symbol in ['gamma','me','mi','qe','qi'] and \
                not symbol in self.params.keys()):
                   raise RuntimeError('%s unspecified'%symbol)
            else:
                if not symbol in self.Fields() + self.ExtendedFields() +\
                        self.params.keys():
                    raise RuntimeError('%s is not an eligible quantity'%symbol)
            salted_symbol = salt+'_'+symbol
            # print(salted_symbol)
            if not salted_symbol in local_dict:
                this_fld = getter(symbol)
                local_dict[salted_symbol] = this_fld
                if len(flds) == 0:
                    if symbol in self.Fields() + self.ExtendedFields() +\
                            self.params.keys():
                        flds.append(this_fld)
                    else:
                        raise RuntimeError('Reduced to scalar.')
            return salted_symbol
        
        salted_eqn = re.sub(var_re, var_salter, eqn)
        
        return ne.evaluate(salted_eqn, local_dict=local_dict,
                        global_dict={'__builtins__': {}})
    def getScalar(self, fld, r=None, **kwargs):
        """
        Args:
            fld     : name of scalar
            r       : if r is a rectilinear grid, then set r.point_data.scalars
                      and r.dimensions
            **kwargs: slicing info passed to getField
        Returns:
            vtkDoubleArray
        """
        import vtk
        from vtk.util import numpy_support
        arr_numpy = self.getField(fld, transpose=False, **kwargs)
        arr_vtk = nump_support.numpy_to_vtk(
                arr_numpy.ravel(order='F'),
                deep=True, array_type=vtk.VTK_DOUBLE)
        arr_vtk.SetName(fld)
        if type(r) == tvtk.tvtk_classes.rectilinear_grid.RectilinearGrid:
            r.point_data.scalars = arr_vtk
            r.dimensions = arr_numpy.shape
        del(arr_numpy)
        return arr_vtk
    def getVector3(self, flds, name=None, r=None, **kwargs):
        """
        Args:
            flds    : list containing names of all three components
            name    : name of vector
            r       : if r is a rectilinear grid, then set r.point_data.vectors
                      and r.dimensions
            **kwargs: slicing info passed to getField
        Returns:
            vtkDoubleArray
        """
        import vtk
        from vtk.util import numpy_support
        import numpy as np
        u = self.getField(flds[0], transpose=False, **kwargs)
        v = self.getField(flds[1], transpose=False, **kwargs)
        w = self.getField(flds[2], transpose=False, **kwargs)
        arr_vtk = numpy_support.numpy_to_vtk(
                np.ascontiguousarray(np.array([
                    u.ravel(order='F'),
                    v.ravel(order='F'),
                    w.ravel(order='F')
                    ]).T),
                deep=True, array_type=vtk.VTK_DOUBLE,
                )
        if name != None:
            arr_vtk.SetName(name)
        import tvtk
        if type(r) == tvtk.tvtk_classes.rectilinear_grid.RectilinearGrid:
            r.point_data.vectors = arr_vtk
            r.dimensions = u.shape
        del(u,v,w)
        return arr_vtk
    def getRectilinearGrid(self, scalar=None, vector=None, vector_name=None, **kwargs):
        from tvtk.api import tvtk
        r = tvtk.RectilinearGrid()
        if scalar!=None:
            self.getScalar(scalar, r, **kwargs)
        if vector!=None:
            self.getVector3(vector, vector_name, r, **kwargs)
        r.x_coordinates = self.getCoordinates('x') # TODO: slicing
        r.y_coordinates = self.getCoordinates('y')
        r.z_coordinates = self.getCoordinates('z')
        return r
