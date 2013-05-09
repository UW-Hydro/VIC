#!/usr/bin/perl
#
# convert_state_file.pl - script to convert VIC state files from the format
#                         used by old versions of VIC to the format used by
#                         the current version of VIC
#
# usage: see usage() function below
#
# Author: Ted Bohn
# $Id: $
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminary Stuff
#-------------------------------------------------------------------------------

# This allows us to use sophisticated command-line argument parsing
use Getopt::Long;

# Default values
$nfrost = 1;
$nband_global = 1;

# This is for validation
@supported_input_versions = ("4.0.3", "4.0.4", "4.0.5", "4.0.6", "4.1.0_r1", "4.1.0_r2", "4.1.0_r3", "4.1.0_r4", "4.1.0_r5", "4.1.1", "4.1.2");

#-------------------------------------------------------------------------------
# Parse the command line
#-------------------------------------------------------------------------------

# Hash used in GetOptions function
# format: option => \$variable_to_set
%options_hash = (
  h   => \$help,
  i   => \$infile,
  o   => \$outfile,
  vi  => \$version_in,
  nfrost  => \$nfrost,
  nband  => \$nband_global,
  snowband  => \$snowband_file,
);

# This parses the command-line arguments and sets values for the variables in %option_hash
$status = &GetOptions(\%options_hash,"h","i=s","o=s","vi=s","nfrost=s","nband=s", "snowband=s");

#-------------------------------------------------------------------------------
# Validate the command-line arguments
#-------------------------------------------------------------------------------

# Help option
if ($help) {
  &usage("full");
  exit(0);
}

# Other options
if ($version_in !~ /^4\.\d\.\d/) {
  printf STDERR "$0: ERROR: the version specified for the input state file ($version_in) is not a valid version descriptor; versions must be in the form 4.x.y, where x and y are integers, and y can be followed by an optional _rz, where z is another integer\n";
  &usage("short");
  exit(1);
}
else {
  $found = 0;
  foreach $supp_ver (@supported_input_versions) {
    if ($version_in eq $supp_ver) {
      $found = 1;
      last;
    }
  }
  if (!$found) {
    printf STDERR "$0: ERROR: the version specified for the input state file ($version_in) is not supported by this script\n";
    &usage("short");
    exit(1);
  }
}
if ($nfrost !~ /^\d+$/) {
  printf STDERR "$0: ERROR: the supplied number of spatial frost fractions ($nfrost) is not valid; you must supply a positive integer; if the SPATIAL_FROST feature is turned off, you must supply a \"1\" for this.\n";
  &usage("short");
  exit(1);
}
if ($nfrost == 1) {
  $nfrost_in = 1;
}
elsif ($nfrost > 1) {
  if ($version_in =~ /^4\.0/) {
    $nfrost_in = 1;
  }
  else {
    $nfrost_in = $nfrost;
  }
}

if ($nband_global !~ /^\d+$/ || $nband_global <= 0) {
  printf STDERR "$0: ERROR: the supplied global number of snow bands ($nband_global) is not valid; you must supply a positive integer; if the SNOWBAND feature is turned off, you must supply a \"1\" for this.\n";
  &usage("short");
  exit(1);
}

#-------------------------------------------------------------------------------
# Get any information necessary to parse the state file
#-------------------------------------------------------------------------------

if ($nband_global > 1 && $version_in eq "4.0.6") {
  # Read snowband param file and get relevant information for reading state file
  open (SNOW, $snowband_file) or die "$0: ERROR: cannot open snowband parameter file $snowband_file for reading\n";
  foreach (<SNOW>) {
    chomp;
    @fields = split /\s+/;
    $cellid = $fields[0];
    $count_nonzero = 0;
    for ($i=0; $i<$nband_global; $i++) {
      $area = $fields[$i+1];
      if ($area > 0) {
        $count_nonzero++;
      }
    }
    $nband_per_cell{$cellid} = $count_nonzero;
  }
  close(SNOW);
}

#-------------------------------------------------------------------------------
# Read the state file and print the re-formatted information
#-------------------------------------------------------------------------------

# Open input state file
open (STATE, $infile) or die "$0: ERROR: cannot open state file $infile for reading\n";

# Open output state file
open (OUT, ">$outfile") or die "$0: ERROR: cannot open output state file $outfile for writing\n";

# Read state file "header"
$line = <STATE>; # Date
print OUT $line;
$line = <STATE>; # Nlayer, Nnodes
chomp $line;
($nlayer,$nnode) = split /\s+/, $line;
print OUT "$line\n";

# Read the rest of the state file
$cell_header = 1;
foreach (<STATE>) {
  chomp;
  @fields = split /\s+/;
  if ($cell_header) {
    $cellid = shift @fields;
    $nveg = shift @fields;
    if ($version_in =~ /4\.1\.0_r(3|4|5)/) {
      $extra_veg = shift @fields;
    }
    $nband = shift @fields;
    for ($i=0; $i<$nnode; $i++) {
      $dz_node[$i] = $fields[$i];
    }
    if ($version_in =~ /4\.1/ && ($version_in !~ /4\.1\.0/ || $version_in =~ /r(4|5)/)) {
      for ($i=0; $i<$nnode; $i++) {
        $dz_sum[$i] = $fields[$i+$nnode];
      }
    }
    else {
      $dz_sum[0] = 0;
      $sum = 0;
      for ($i=1; $i<$nnode; $i++) {
        $sum += ($dz_node[$i]+$dz_node[$i-1])*0.5;
        $dz_sum[$i] = $sum;
      }
    }
    print OUT "$cellid $nveg $nband";
    for ($i=0; $i<$nnode; $i++) {
      print OUT " $dz_node[$i]";
    }
    for ($i=0; $i<$nnode; $i++) {
      print OUT " $dz_sum[$i]";
    }
    print OUT "\n";
    $cell_header = 0;
    $veg = 0;
    $band = 0;
    if ($nband_global > 1 && $version_in eq "4.0.6") {
      $nband_this_cell = $nband_per_cell{$cellid};
    }
    else {
      $nband_this_cell = $nband_global;
    }
    if ($version_in =~ /4\.0/) {
      $distprec_line = 1;
    }
    else {
      $veg_header = 1;
    }
  }
  elsif ($distprec_line) {
    ($still_storm,$dry_time) = @fields;
    $distprec_line = 0;
    $veg_header = 1;
  }
  elsif ($veg_header) {
    $mu = shift @fields;
    if ($version_in =~ /4\.1/) {
      $still_storm = shift @fields;
      $dry_time = shift @fields;
    }
    print OUT "$mu $still_storm $dry_time\n";
    $veg_header = 0;
    $vegband_block = 1;
  }
  elsif ($vegband_block) {

    $veg_tmp = shift @fields;
    $band_tmp = shift @fields;
    print OUT "$veg_tmp $band_tmp";
    for ($i=0; $i<$nlayer; $i++) {
      $tmp = shift @fields;
      print OUT " $tmp";
    }
    for ($i=0; $i<$nlayer; $i++) {
      for ($j=0; $j<$nfrost_in; $j++) {
        $tmp = shift @fields;
        print OUT " $tmp";
      }
      for ($j=$nfrost_in; $j<$nfrost; $j++) {
        print OUT " $tmp";
      }
    }
    if ($veg < $nveg) {
      $tmp = shift @fields;
      print OUT " $tmp";
    }
    for ($i=0; $i<11+$nnode; $i++) {
      print OUT " $fields[$i]";
    }
    print OUT "\n";

    # Increment band and veg indexes
    $band++;
    if ($band == $nband_this_cell) {

      # Fill in 0-area bands that were omitted from 4.0.6 state file
      if ($version_in eq "4.0.6") {
        while ($band < $nband_global) {
          print OUT "$veg_tmp $band";
          for ($i=0; $i<$nlayer; $i++) {
            print OUT " 0.0";
          }
          for ($i=0; $i<$nlayer; $i++) {
            for ($j=0; $j<$nfrost; $j++) {
              print OUT " 0.0";
            }
          }
          if ($veg < $nveg) {
             print OUT " 0.0";
          }
          for ($i=0; $i<11+$nnode; $i++) {
            print OUT " 0.0";
          }
          print OUT "\n";
          $band++;
        }
      }

      $band = 0;
      $veg++;
      $vegband_block = 0;

      if ($veg > $nveg) {
        $cell_header = 1;
      }
      else {
        $veg_header = 1;
      }

    }

  }
#how to handle lakes
}
close(STATE);

#-------------------------------------------------------------------------------
# Usage
#-------------------------------------------------------------------------------

sub usage() {

  print "\n";
  print "$0: script to convert VIC state files from the format used by old versions of VIC to the format used by the current version of VIC\n";
  print "\n";
  print "usage:\n";
  print "  $0 [-h] -i <infile> -vi <version_in> [-nfrost <nfrost>] [-nband <nband> [-snowband <snowband_file>]] -o <outfile> \n";
  print "\n";
  if ($_[0] eq "full") {
    print "  -h\n";
    print "    prints this usage message\n";
    print "\n";
    print "  -i <infile>\n";
    print "    <infile>  = filename of input VIC state file to be converted.\n";
    print "\n";
    print "  -vi <version_in>\n";
    print "    <version_in>  = the VIC release version corresponding to the input state file; should be of the form \"4.x.y[_rz]\", where x, y, and z are integers and the \"_rz\" is optional, only occurring in some VIC releases.\n";
    print "    The VIC release string for your copy of VIC can be found by typing \"vicNl -v\".\n";
    print "\n";
    print "  -o <outfile>\n";
    print "    <outfile>  = filename of output VIC state file.\n";
    print "\n";
    print "  -nfrost <nfrost>\n";
    print "    (optional) this only needs to be specified if the input state file was created by VIC 4.1.0 or later with SPATIAL_FROST = TRUE.\n";
    print "    <nfrost>  = Number of bins in the spatial frost distribution\n";
    print "\n";
    print "  -nband <nband>\n";
    print "    (optional) this only needs to be specified if the input state file was created with SNOWBANDS = TRUE.\n";
    print "    <nband>  = *Global* number of snow elevation bands, i.e. the maximum number of snow bands allowed for any grid cell, also the number of area fields in the snowband parameter file.\n";
    print "\n";
    print "  -snowband <snowband_file>\n";
    print "    (optional) this only needs to be specified if nbands was specified and is > 1\n";
    print "    <snowband_file>  = filename of the snowband parameter file.\n";
  }
}
