c ====================================================================
c In editing this file : please make sure there are spaces between
c each word and character in the lines defining the variables.
c This is in order to allow the program calc_memory.c to run.
c Thus keep the format:
c parameter ( MAX_VAR = x )
c Also make sure parameter is spelled in lower case.
c Do not change the line with SNW_FLG.
c $Id$
c ====================================================================

      integer MAX_PIX,MAX_SPP,MAX_STA,MAX_SOI,MAX_VEG
      integer MAX_CAT,MAX_ROW,MAX_COL,MAX_FIL,MAX_SER
      integer MAX_ATB,MAX_TST,MAX_VST,MAX_TOP,MAX_LAN
      integer MAX_PP1,MAX_PP2,MOS_FLG,UST_FLG,SNW_FLG
      integer LAK_FLG

      parameter ( MAX_PIX = 3960 )
      parameter ( MAX_PP1 = 3960 )
      parameter ( MAX_PP2 = 1 )
      parameter ( MAX_SPP = 4 )
      parameter ( MAX_STA = 253 )
      parameter ( MAX_SOI = 376 )
      parameter ( MAX_VEG = 20 )
      parameter ( MAX_VST = 3 )
      parameter ( MAX_LAN = 1 )
      parameter ( MAX_CAT = 101 )
      parameter ( MAX_ROW = 61 )
      parameter ( MAX_COL = 67 )
      parameter ( MAX_FIL = 400 )
      parameter ( MAX_SER = 2 )
      parameter ( MAX_ATB = 10 )
      parameter ( MAX_TOP = 1 )
      parameter ( MAX_TST = 1 )
      parameter ( MOS_FLG = 1 )
      parameter ( UST_FLG = 1 )
      parameter ( LAK_FLG = 1 )

      parameter ( SNW_FLG = SNOW_RUN*(MOS_FLG+UST_FLG-MOS_FLG*UST_FLG) )

c ====================================================================
c Parameters for array dimensions for weighting and station files.
c
c MAX_PIX - Maximum number of pixels.
c MAX_PP1 - Maximum number of pixels in statistical mode.  If running
c	    in statistical mode this number should equal MAX_PIX.
c	    In distributed mode this number should equal 1.
c MAX_PP2 - Maximum number of pixels in statistical mode for fractional
c	    coverage.  If running in distributed mode or in statistical
c	    mode without representation of rainfall fractional coverage
c	    this number should equal 1.  If rainfall fractional coverage
c	    is represented this number should equal MAX_PIX.
c MAX_SPP - Maximum number of stations used per pixel.
c MAX_STA - Maximum number of stations for each forcing variable.
c MAX_SOI - Maximum number of soil classes.
c MAX_VEG - Maximum number of land cover classes.
c MAX_VST - Maximum number of vegetation classes per pixel in
c	    statistical mode.
c	    If running the model in statistical mode, this number
c	    should equal MAX_VEG or less.
c	    If running the model in distributed mode, this number
c	    should equal 1.
c MAX_LAN - Maximum number of vegetation classes per pixel.
c	    If running the model in distributed mode, this number
c	    should equal 1.
C	    If considering fractional coverage of rainfall than this
c	    number should equal MAX_VST.
c	    If assuming rainfall to be uniformly distributed over
c	    each grid than this paramters should equal 1.
c MAX_CAT - Maximum number of catchments.
c MAX_ROW - Maximum number of rows in the image.
c MAX_COL - Maximum number of columns in the image.
c MAX_FIL - Highest input or output file number.
c MAX_SER - Maximum number of series of printing the output
c           images (e.g. one serie is from time step 2 through 7, a
c           second serie from time step 25 through 33, ...).
c MAX_ATB - Maximum number of intervals in the atb distributions.
c MAX_TOP - Maximum number of topindex intervals per pixel.
C	    If considering fractional coverage of rainfall than this
c	    number should equal MAX_ATB.
c	    If assuming rainfall to be uniformly distributed over
c	    each grid than this paramters should equal 1.
c MAX_TST - Maximum number of time steps to be solved for.
c MOS_FLG - 1 if there is a moss layer in at least one land cover
c	    class.  zero if all land cover classes have no moss layer.
c UST_FLG - 1 if there is an understory layer in at least one land cover
c	    class.  zero if all land cover classes have no understory layer.
c LAK_FLG - 1 if there is at least one land cover class that is a lake.
c           0 if none of the land cover classes are lakes.
c ====================================================================
