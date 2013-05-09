#!/usr/bin/perl -w
# set_moist.pl - Sets initial soil moisture in VIC soil file
# Author: T. Bohn, edited A. Wood may 2006
#
# usage:
#	set_moist.pl soil_file nlayers init_moist_level > new_soil_file
# where
#	soil_file	= name of input VIC soil file
#	nlayers		= number of soil layers
#	init_moist_level   = desired initial soil moisture level
#                               options: "max","fc","wp", fraction [0.0-1.0] in range wp to max
#	new_soil_file	= name of output VIC soil file
#
# This assumes that thickness, bulk_density, and soil_density are correct
# in your input soil file.  In particular, soil_density must be non-zero.
#
# you can run fix_resid.pl to fix those first.
#----------------------------------------------------------------------------

# Get cmdline parameters
$soilfile = shift;
$nlayer = shift;
$init_moist_level = shift;

# Read soil param file and compute new values
open (FILE, $soilfile) or die "Error: cannot open $soilfile\n";
foreach (<FILE>) {
  chomp;
  @fields = split /\s+/;
  for ($idx=0; $idx<$nlayer; $idx++) {

    # Read values
    $thickness[$idx] = $fields[4*$nlayer+10+$idx];
    $bulkdens[$idx] = $fields[7*$nlayer+12+$idx];
    $soildens[$idx] = $fields[8*$nlayer+12+$idx];
    $Wcr_FRACT[$idx] = $fields[9*$nlayer+13+$idx];
    $Wpwp_FRACT[$idx] = $fields[10*$nlayer+13+$idx];

    # Compute soil moisture
    # Compute max moisture: thickness(m) * porosity * 1000(mm/m)
    $max_moist = $thickness[$idx] * (1 - $bulkdens[$idx]/$soildens[$idx]) * 1000;
    if ($init_moist_level eq "max") {
      # set soil moisture to max possible
      $init_moist[$idx] = $max_moist;
    }
    elsif ($init_moist_level eq "fc") {
      # set soil moisture to field capacity
      $init_moist[$idx] = $max_moist * $Wcr_FRACT[$idx] / 0.7;  # Wcr_FRACT is 70% of field capacity
    }
    elsif ($init_moist_level eq "wp") {
      # set soil moisture to wilting point
      $init_moist[$idx] = $max_moist * $Wpwp_FRACT[$idx];
    }
    elsif ($init_moist_level =~ /\d/ && $init_moist_level >= 0.0 && $init_moist_level <= 1.0) {
      # set soil moisture to fraction between max moisture and wilting point
      $init_moist[$idx] = ($max_moist * $Wpwp_FRACT[$idx]) +
        ($init_moist_level * ($max_moist-($max_moist * $Wpwp_FRACT[$idx])));
    }
    else {
      print STDERR "$0: ERROR: unrecognized init_moist_level $init_moist_level\n";
      print STDERR "   must be one of max, fc, wp, or a fraction (0.0-1.0) of the range between wp and max\n";
      die;
    }

    # Set precision
    $init_moist[$idx] *= 1000;
    $init_moist[$idx] = int $init_moist[$idx];
    $init_moist[$idx] /= 1000;

  }

  # Print out new values
  for ($idx=0; $idx<(3*$nlayer+9); $idx++) {
    print "$fields[$idx] ";
  }
  for ($idx=0; $idx<$nlayer; $idx++) {
    print "$init_moist[$idx] ";
  }
  for ($idx=(4*$nlayer+9); $idx<(7*$nlayer+12); $idx++) {
    print "$fields[$idx] ";
  }
  for ($idx=0; $idx<$nlayer; $idx++) {
    print "$bulkdens[$idx] ";
  }
  for ($idx=0; $idx<$nlayer; $idx++) {
    print "$soildens[$idx] ";
  }
  for ($idx=(9*$nlayer+12); $idx<$#fields; $idx++) {
    print "$fields[$idx] ";
  }
  print "$fields[$#fields]\n";
}
close (FILE);
