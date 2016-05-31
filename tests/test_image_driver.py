#!/usr/bin/env python 
''' VIC Image Driver testing '''

import xarray as xr
import numpy as np 
import numpy.testing as npt 

def test_image_driver_no_output_file_nans(fnames, domain_file):
    '''
    Test that all VIC image driver output files have the same nan structure as
    the domain file
    '''
    for fname in fnames:
        ds_domain = xr.open_dataset(domain_file)
        ds_output = xr.open_dataset(fname)
        assert_nan_equal(ds_domain, ds_output)

def assert_nan_equal(ds_domain, ds_output): 
    """ 
    test to see if nans in two data arrays are the same in order to check domain 
    first dataarray is domain file, second dataarray is a variable in the output
    dataset for image driver 
    """

    # get lat/lng coordinate names from domain file mask 
    match_vars = ds_domain['mask'].coords
    
    # check to be sure that mask, lats and lons match between domain file and output file
    
    for var in match_vars: 
        # raise AssertionError if they don't match
        npt.assert_allclose(ds_output[var], ds_domain[var], equal_nan=True) 

    # check that nans are occurring in the same place in the arrays 
    
    # check all variables in the dataset 
    for da in ds_output.data_vars: 
        
        # get dimensions to reduce DataArray on
        dim_diff = set(ds_domain['mask'].dims).symmetric_difference(set(ds_output[da].dims))
        
        # reduce DataArray
        da_null_reduced = ds_output[da].isnull().all(dim=dim_diff)
        
        # raise AssertionError if NaNs do not match 
        npt.assert_array_equal(da_null_reduced.values, np.isnan(ds_domain['mask'])) 
