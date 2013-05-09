#!/usr/bin/perl
#
# fix_resid.pl - Does quality check on residual moisture in VIC soil file.
#                Also supplies default values for bulk_density, soil_density,
#                Wcr_FRACT, and Wpwp_FRACT if they are missing.
#
# author: Ted Bohn, some time in 2005 (around the release of VIC 4.0.5)
#         tbohn@hydro.washington.edu
#
# usage:
#	fix_resid.pl soil_file nlayers > new_soil_file
# where
#	soil_file	= name of input VIC soil file
#	nlayers		= number of soil layers
#	new_soil_file	= name of output VIC soil file
#
# The quality check has 3 parts:
# 1. Check for zeros in bulk density, soil density, Wcr_FRACT, and Wpwp_FRACT.
#    Replace zeros with default values.  Also, if Wcr_FRACT or Wpwp_FRACT are > 1.0,
#    replace them with default values (you can change these default values if you want).
# 2. Make sure Wpwp_FRACT <= 0.99*Wcr_FRACT.
# 3. Make sure residual moisture <= 0.99*Wpwp_fract
#
# If you would like to use different default values for bulk_density,
# soil_density, Wcr_FRACT, and Wpwp_FRACT, change them below.
#----------------------------------------------------------------------------

# Get cmdline parameters
$soilfile = shift;
$nlayer = shift;

# Default values, in case any are zero
for ($idx=0; $idx<$nlayer; $idx++) {
  if ($idx == 0) {
    # Defaults for top soil layer
    $def_bulkdens[$idx] = 1300;
    $def_soildens[$idx] = 2500;
    $def_Wcr_FRACT[$idx] = 0.42;
    $def_Wpwp_FRACT[$idx] = 0.28;
  }
  else {
    # Defaults for all other layers
    $def_bulkdens[$idx] = 1500;
    $def_soildens[$idx] = 2700;
    $def_Wcr_FRACT[$idx] = 0.43;
    $def_Wpwp_FRACT[$idx] = 0.28;
  }
}

# Read soil param file and compute new values
open (FILE, $soilfile) or die "Error: cannot open $soilfile\n";
foreach (<FILE>) {
  chomp;
  @fields = split /\s+/;
  $cellid = $fields[1];
  for ($idx=0; $idx<$nlayer; $idx++) {

    # Read values
    $bulkdens[$idx] = $fields[7*$nlayer+12+$idx];
    $soildens[$idx] = $fields[8*$nlayer+12+$idx];
    $Wcr_FRACT[$idx] = $fields[9*$nlayer+13+$idx];
    $Wpwp_FRACT[$idx] = $fields[10*$nlayer+13+$idx];
    $resid[$idx] = $fields[11*$nlayer+16+$idx];

    # Handle zeros
    if ($bulkdens[$idx] <= 0) {
      printf STDERR "%s: WARNING: cell %d invalid bulk density for layer %d (%f); setting to %f\n", $0, $cellid, ($idx+1), $bulkdens[$idx], $def_bulkdens[$idx];
      $bulkdens[$idx] = $def_bulkdens[$idx];
    }
    if ($soildens[$idx] <= 0) {
      printf STDERR "%s: WARNING: cell %d invalid soil density for layer %d (%f); setting to %f\n", $0, $cellid, ($idx+1), $soildens[$idx], $def_soildens[$idx];
      $soildens[$idx] = $def_soildens[$idx];
    }
    if ($Wcr_FRACT[$idx] <= 0 || $Wcr_FRACT[$idx] > 1.0) {
      printf STDERR "%s: WARNING: cell %d invalid Wcr_FRACT for layer %d (%f); setting to %f\n", $0, $cellid, ($idx+1), $Wcr_FRACT[$idx], $def_Wcr_FRACT[$idx];
      $Wcr_FRACT[$idx] = $def_Wcr_FRACT[$idx];
    }
    if ($Wpwp_FRACT[$idx] <= 0 || $Wpwp_FRACT[$idx] > 1.0) {
      printf STDERR "%s: WARNING: cell %d invalid Wpwp_FRACT for layer %d (%f); setting to %f\n", $0, $cellid, ($idx+1), $Wpwp_FRACT[$idx], $def_Wpwp_FRACT[$idx];
      $Wpwp_FRACT[$idx] = $def_Wpwp_FRACT[$idx];
    }

    # Make sure Wpwp_FRACT doesn't exceed 0.99*Wcr_FRACT
    if ($Wpwp_FRACT[$idx] > 0.99*$Wcr_FRACT[$idx]) {
      printf STDERR "%s: WARNING: cell %d layer %d Wpwp_FRACT (%f) is greater than 0.99*Wcr_FRACT (Wcr_FRACT == %f); setting Wpwp_FRACT to 0.99*Wcr_FRACT\n", $0, $cellid, ($idx+1), $Wpwp_FRACT[$idx], $Wcr_FRACT[$idx];
      $Wpwp_FRACT[$idx] = 0.99*$Wcr_FRACT[$idx];
      # Truncate to desired precision
      $Wpwp_FRACT[$idx] *= 10000;
      $Wpwp_FRACT[$idx] = int $Wpwp_FRACT[$idx];
      $Wpwp_FRACT[$idx] /= 10000;
    }

    # Compute wilting point
    $Wpwp[$idx] = $Wpwp_FRACT[$idx] * (1 - $bulkdens[$idx]/$soildens[$idx]);

    # Make sure residual moisture doesn't exceed 0.99*Wpwp
    if ($resid[$idx] > 0.99*$Wpwp[$idx]) {
      printf STDERR "%s: WARNING: cell %d layer %d residual moisture (%f) is greater than 0.99*Wpwp (Wpwp == %f); setting residual moisture to 0.99*Wpwp\n", $0, $cellid, ($idx+1), $resid[$idx], $Wpwp[$idx];
      $resid[$idx] = 0.99*$Wpwp[$idx];
      # Truncate to desired precision
      $resid[$idx] *= 10000;
      $resid[$idx] = int $resid[$idx];
      $resid[$idx] /= 10000;
    }

  }

  # Print out new values
  for ($idx=0; $idx<(7*$nlayer+12); $idx++) {
    print "$fields[$idx] ";
  }
  for ($idx=0; $idx<$nlayer; $idx++) {
    print "$bulkdens[$idx] ";
  }
  for ($idx=0; $idx<$nlayer; $idx++) {
    print "$soildens[$idx] ";
  }
  for ($idx=(9*$nlayer+12); $idx<(9*$nlayer+13); $idx++) {
    print "$fields[$idx] ";
  }
  for ($idx=0; $idx<$nlayer; $idx++) {
    print "$Wcr_FRACT[$idx] ";
  }
  for ($idx=0; $idx<$nlayer; $idx++) {
    print "$Wpwp_FRACT[$idx] ";
  }
  for ($idx=(11*$nlayer+13); $idx<(11*$nlayer+16); $idx++) {
    print "$fields[$idx] ";
  }
  for ($idx=0; $idx<$nlayer; $idx++) {
    print "$resid[$idx] ";
  }
  print "$fields[$#fields]\n";
}
close (FILE);
