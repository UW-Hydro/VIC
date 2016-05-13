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
    match_vars = ds_domain.coords
    
    # check to be sure that mask, lats and lons match between domain file and output file
    
    for var in match_vars: 
        # raise AssertionError if they don't match
        npt.assert_allclose(ds_output[var], ds_domain[var], equal_nan=True) 

    # check that nans are occurring in the same place in the arrays 
    
    # check all variables in the dataset 
    for da in ds_output.data_vars: 
        
        if (len(ds_output[da].dims) > 3):
            # check all layers and timesteps 
            da_isnull_reduced = ds_output[da].isnull().all(dim=('time', 'nlayer')) 
            npt.assert_array_equal(da_isnull_reduced.values, 
                                    np.isnan(ds_domain['mask']))

        else: 
            # check all timesteps 
            da_isnull_reduced = ds_output[da].isnull().all(dim=('time')) 
            npt.assert_array_equal(da_isnull_reduced.values, 
                                    np.isnan(ds_domain['mask'])) 
