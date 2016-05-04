#!/usr/bin/env python 
''' VIC Image Driver testing '''

import xarray as xr
import numpy as np 
import numpy.testing as npt 

def assert_nan_equal(ds_domain, ds_output): 
    """ 
    test to see if nans in two data arrays are the same in order to check domain 
    first dataarray is domain file, second dataarray is a variable in the output
    dataset for image driver 
    """

    # get lat/lng dim names from domain file 
    lng_name, lat_name = ds_domain.coords

    match_vars = [lng_name, lat_name] 
    
    # check to be sure that mask, lats and lons match between domain file and output file
    for var in match_vars: 
        # raise AssertionError if they don't match
        npt.assert_allclose(ds_output[var], ds_domain[var], equal_nan=True) 

    # check that nans are occurring in the same place in the arrays 
    
    # check all variables in the dataset 
    for da in ds_output.data_vars: 
        
        # check all timesteps 
        for ts in range(0, len(ds_output.time)):
     
            # slice time step of output DataArray
            da_sel = ds_output[da].isel(time=ts) 

            # check all layers 
            if (len(ds_output[da].dims) > 3):  
                for layer in range(0, len(ds_output[da].nlayer)):
                    # npt.assert_allclose(da_sel.isel(nlayer=layer), ds_domain['mask'], equal_nan=True)
                    npt.assert_array_equal(np.isnan(da_sel.isel(nlayer=layer).values), np.isnan(ds_domain['mask'].values))                    
                    
            else: 
                # npt.assert_allclose(da_sel, ds_domain['mask'], equal_nan=True)
                npt.assert_array_equal(np.isnan(da_sel.values), np.isnan(ds_domain['mask'].values))
