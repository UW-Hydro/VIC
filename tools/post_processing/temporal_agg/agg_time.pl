#!/usr/bin/perl
#
# agg_time.pl - script to temporally aggregate ASCII data files
#
# usage: see usage() function below
#
# Author: Ted Bohn
# $Id: agg_time.pl,v 1.5 2005/08/27 06:32:15 vicadmin Exp $
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminary Stuff
#-------------------------------------------------------------------------------

# This allows us to use sophisticated command-line argument parsing
use Getopt::Long;

# This is for validation
@supported_input_intervals = ("minutely", "hourly", "daily", "monthly");
@supported_output_intervals = ("hourly", "daily", "monthly", "yearly", "annual");

#-------------------------------------------------------------------------------
# Parse the command line
#-------------------------------------------------------------------------------

# Hash used in GetOptions function
# format: option => \$variable_to_set
%options_hash = (
  h   => \$help,
  i   => \$infile,
  o   => \$outfile,
  in  => \$in_interval,
  out => \$out_interval,
  # There should be an entry for each aggregation method here
  avg => \@fields_avg,    # array of fields to be averaged
  sum => \@fields_sum,    # array of fields to be summed
  beg => \@fields_beg,    # array of fields for which we store value at beginning of interval
  end => \@fields_end,    # array of fields for which we store value at end of interval
  min => \@fields_min,    # array of fields for which we store min value within interval
  max => \@fields_max,    # array of fields for which we store max value within interval
);

# This parses the command-line arguments and sets values for the variables in %option_hash
$status = &GetOptions(\%options_hash,"h","i=s","o=s","in=s","out=s",
  "avg=s@","sum=s@","beg=s@","end=s@","min=s@","max=s@");

#-------------------------------------------------------------------------------
# Validate the command-line arguments
#-------------------------------------------------------------------------------

# Help option
if ($help) {
  usage("full");
  exit(0);
}

# Validate input/output files
if (!$infile || !$outfile) {
  printf STDERR "$0: ERROR: input and output filenames must be specified.\n";
  usage("short");
  exit(1);
}

# Validate the input interval
if (!$in_interval) {
  $in_interval = "daily";
}
else {
  $supported = 0;
  foreach $interval (@supported_input_intervals) {
    if ($in_interval =~ /^$interval$/i) {
      $in_interval = $interval;
      $supported = 1;
    }
  }
  if (!$supported) {
    printf STDERR "$0: ERROR: specified input interval \"$in_interval\" not supported.\n";
    usage("short");
    exit(1);
  }
}

# Validate the aggregation interval
if (!$out_interval) {
  $out_interval = "monthly";
}
else {
  $supported = 0;
  foreach $interval (@supported_output_intervals) {
    if ($out_interval =~ /^$interval$/i) {
      $out_interval = $interval;
      $supported = 1;
    }
  }
  if (!$supported) {
    printf STDERR "$0: ERROR: specified aggregation interval \"$out_interval\" not supported.\n";
    usage("short");
    exit(1);
  }
}
if ($out_interval eq "annual") {
  $out_interval = "yearly";
}
if (($in_interval eq "daily" || $in_interval_interval eq "monthly") && $out_interval eq "daily") {
  printf STDERR "$0: ERROR: we cannot aggregate to daily interval if input data is not sub-daily.\n";
  usage("short");
  exit(1);
}

# Determine first data field (indexing starting at 0)
if ($in_interval eq "monthly") {
  $first_data_field = 2;
}
elsif ($in_interval eq "daily") {
  $first_data_field = 3;
}
elsif ($in_interval eq "hourly") {
  $first_data_field = 4;
}
else {
  $first_data_field = 5;
}

# Validate fields and assign aggregation methods to them
assign_agg_methods("avg", \@fields_avg, $first_data_field, \%agg_method);
assign_agg_methods("sum", \@fields_sum, $first_data_field, \%agg_method);
assign_agg_methods("beg", \@fields_beg, $first_data_field, \%agg_method);
assign_agg_methods("end", \@fields_end, $first_data_field, \%agg_method);
assign_agg_methods("min", \@fields_min, $first_data_field, \%agg_method);
assign_agg_methods("max", \@fields_max, $first_data_field, \%agg_method);

#-------------------------------------------------------------------------------
# Aggregation
#-------------------------------------------------------------------------------

# Open output file for writing
open (OUTFILE, ">$outfile") or die "$0: ERROR: cannot open $outfile\n";

# Read input file and perform aggregation
open (INFILE, $infile) or die "$0: ERROR: cannot open $infile\n";
$first_line = 1;
$num_recs = 0;
foreach (<INFILE>) {

  # Remove leading spaces
  s/^\s+//;

  # Parse line
  @fields = split /\s+/;
  ($year, $month) = ($fields[0], $fields[1]);
  if ($in_interval eq "daily"
    || $in_interval eq "hourly"
    || $in_interval eq "minutely") {
    $day = $fields[2];
  }
  if ($in_interval eq "hourly"
    || $in_interval eq "minutely") {
    $hour = $fields[3];
  }
  if ($in_interval eq "minutely") {
    $minute = $fields[4];
  }

  # If this is the first line of the file, do some initialization.
  if ($first_line) {

    # Initialize the aggregation array, now that we know how many
    # fields are in a line
    for ($i=$first_data_field; $i<=$#fields; $i++) {
      push (@agg_data, 0);
    }

    # Initialize start-of-interval date/time
    $start_year = $year;
    $start_month = $month;
    if ($in_interval eq "daily"
      || $in_interval eq "hourly") {
      $start_day = $day;
    }
    if ($in_interval eq "hourly") {
      $start_hour = $hour;
    }

    # Initialize $prev_*
    $prev_year = $year;
    $prev_month = $month;
    if ($in_interval eq "daily"
      || $in_interval eq "hourly"
      || $in_interval eq "minutely") {
      $prev_day = $day;
    }
    if ($in_interval eq "hourly"
      || $in_interval eq "minutely") {
      $prev_hour = $hour;
    }
    if ($in_interval eq "minutely") {
      $prev_minute = $minute;
    }

    # Turn off the first_line flag
    $first_line = 0;

  }

  # Compare current record's time to previous record's time.
  # If this is the beginning of the next aggregation interval,
  # do end-of-interval stuff
  if ( ($out_interval eq "hourly" && $hour != $prev_hour)
    || ($out_interval eq "daily" && $day != $prev_day)
    || ($out_interval eq "monthly" && $month != $prev_month)
    || ($out_interval eq "yearly" && $year != $prev_year)
    ) {

    # Print interval date
    printf OUTFILE "%04d", $start_year;
    if ($out_interval eq "monthly"
      || $out_interval eq "daily"
      || $out_interval eq "hourly") {
      printf OUTFILE " %02d", $start_month;
    }
    if ($out_interval eq "daily"
      || $out_interval eq "hourly") {
      printf OUTFILE " %02d", $start_day;
    }
    if ($out_interval eq "hourly") {
      printf OUTFILE " %02d", $start_hour;
    }

    # Finish aggregation and print data
    for ($i=0; $i<=$#fields-$first_data_field; $i++) {
      if ($agg_method{$i} eq "sum" || $agg_method{$i} eq "beg" || $agg_method{$i} eq "end"
        || $agg_method{$i} eq "min" || $agg_method{$i} eq "max") {
        # We already have what we need, so do nothing
        ;
      }
      else {
        # Default is average
        $agg_data[$i] /= $num_recs;
      }
      printf OUTFILE " %12.4f", $agg_data[$i];
    }
    printf OUTFILE "\n";

    # Re-initialize
    $start_year = $year;
    $start_month = $month;
    if ($in_interval eq "daily"
      || $in_interval eq "hourly"
      || $in_interval eq "minutely") {
      $start_day = $day;
    }
    if ($in_interval eq "hourly"
      || $in_interval eq "minutely") {
      $start_hour = $hour;
    }
    if ($in_interval eq "minutely") {
      $start_minute = $minute;
    }
    for ($i=0; $i<=$#fields-$first_data_field; $i++) {
      $agg_data[$i] = 0;
    }
    $num_recs = 0;

  }

  # Start aggregation for each field
  for ($i=0; $i<=$#fields-$first_data_field; $i++) {
    if ($agg_method{$i} eq "beg") {
      # Instantaneous value at beginning of interval
      if ($num_recs == 0) {
        $agg_data[$i] = $fields[$i+$first_data_field];
      }
    }
    elsif ($agg_method{$i} eq "end") {
      # Instantaneous value at end of interval; keep updating until final value
      $agg_data[$i] = $fields[$i+$first_data_field];
    }
    elsif ($agg_method{$i} eq "min") {
      # Min value within interval
      if ($num_recs == 0) {
        $agg_data[$i] = $fields[$i+$first_data_field];
      }
      elsif ($fields[$i+$first_data_field] < $agg_data[$i]) {
        $agg_data[$i] = $fields[$i+$first_data_field];
      }
    }
    elsif ($agg_method{$i} eq "max") {
      # Max value within interval
      if ($num_recs == 0) {
        $agg_data[$i] = $fields[$i+$first_data_field];
      }
      elsif ($fields[$i+$first_data_field] > $agg_data[$i]) {
        $agg_data[$i] = $fields[$i+$first_data_field];
      }
    }
    else {
      # Avg and sum require keeping a running total
      $agg_data[$i] += $fields[$i+$first_data_field];
    }
  }

  # Update $prev_*
  $prev_year = $year;
  $prev_month = $month;
  if ($in_interval eq "daily"
    || $in_interval eq "hourly"
    || $in_interval eq "minutely") {
    $prev_day = $day;
  }
  if ($in_interval eq "hourly"
    || $in_interval eq "minutely") {
    $prev_hour = $hour;
  }
  if ($in_interval eq "minutely") {
    $prev_minute = $minute;
  }

  # Update $num_recs
  $num_recs++;

}
close (INFILE);

# Handle the end of the final interval

# Print interval date
printf OUTFILE "%04d", $start_year;
if ($out_interval eq "monthly"
  || $out_interval eq "daily"
  || $out_interval eq "hourly") {
  printf OUTFILE " %02d", $start_month;
}
if ($out_interval eq "daily"
  || $out_interval eq "hourly") {
  printf OUTFILE " %02d", $start_day;
}
if ($out_interval eq "hourly") {
  printf OUTFILE " %02d", $start_hour;
}

# Finish aggregation and print data
for ($i=0; $i<=$#fields-$first_data_field; $i++) {
  if ($agg_method{$i} eq "sum" || $agg_method{$i} eq "beg" || $agg_method{$i} eq "end"
    || $agg_method{$i} eq "min" || $agg_method{$i} eq "max") {
    # We already have what we need, so do nothing
    ;
  }
  else {
    # Default is average
    $agg_data[$i] /= $num_recs;
  }
  printf OUTFILE " %12.4f", $agg_data[$i];
}
printf OUTFILE "\n";

# Close output file
close (OUTFILE);

#-------------------------------------------------------------------------------
# Subroutines
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Validate fields and assign aggregation methods to them
#-------------------------------------------------------------------------------

sub assign_agg_methods() {

  my ($agg_method, $field_list_ref, $first_data_field, $agg_method_ref) = @_;
  my $string, @fields, $field;

  foreach $string (@$field_list_ref) {

    @fields = split /,/, $string;

    foreach $field (@fields) {

      # Validate the field

      # Make sure it's a positive integer
      if ($field !~ /^\d+$/ || $field < 1) {
        printf STDERR "$0: ERROR: -$agg_method $field: fields must be positive integers.\n";
        usage("short");
        exit(1);
      }

      # Convert field from user's indexing (first field of line has index=1)
      # to our internal indexing (first data field has index=0)
      $field -= $first_data_field + 1;

      # Must be a data field, not a date field
      if ($field < 0) {
        printf STDERR "$0: ERROR: fields containing the date are not valid for aggregation.\n";
        usage("short");
        exit(1);
      }

      # Record this field's aggregation method
      $$agg_method_ref{$field} = $agg_method;

    }

  }

}

#-------------------------------------------------------------------------------
# Usage
#-------------------------------------------------------------------------------

sub usage() {

  print "\n";
  print "$0: script to temporally aggregate ASCII data files\n";
  print "\n";
  print "usage:\n";
  print "  $0 [-h] -i <infile> -o <outfile> [-in <in_interval>] [-out <out_interval>]\n";
  print "    [ [ -<agg_method> <field>[,<field>...] ] [ -<agg_method> <field>[,<field>...] ] ... ]\n";
  print "\n";
  if ($_[0] eq "full") {
    print "  -h\n";
    print "    prints this usage message\n";
    print "\n";
    print "  -i <infile>\n";
    print "    <infile>  = Input filename.  Files should be of the format:\n";
    print "      YEAR MONTH [DAY] [HOUR] DATA DATA DATA\n";
    print "    where\n";
    print "      YEAR  = 4-digit year\n";
    print "      MONTH = 2-digit month\n";
    print "      DAY   = 2-digit day (optional - see below)\n";
    print "      HOUR  = 2-digit hour (optional - see below)\n";
    print "      DATA  = data fields\n";
    print "    Fields may be separated by multiple spaces or tabs.\n";
    print "\n";
    print "  -o <outfile>\n";
    print "    <outfile> = Output filename.\n";
    print "\n";
    print "  -in <in_interval>\n";
    print "    <in_interval> = Input file record interval.  Must be one of the\n";
    print "    following strings (case-insensitive):\n";
    print "      minutely= input time step is any integer number of\n";
    print "                minutes.  Input file has format:\n";
    print "                YEAR MONTH DAY HOUR MIN DATA DATA DATA ...\n";
    print "      hourly  = input time step is sub-daily (any integer number of\n";
    print "                hours).  Input file has format:\n";
    print "                YEAR MONTH DAY HOUR DATA DATA DATA ...\n";
    print "      daily   = input time step is daily.  Input file has format:\n";
    print "                YEAR MONTH DAY DATA DATA DATA ...\n";
    print "      monthly = input time step is monthly.  Input file has format:\n";
    print "                YEAR MONTH DATA DATA DATA ...\n";
    print "    Default: daily.\n";
    print "    NOTE: it is important to specify this option correctly to ensure that\n";
    print "    this script can read the input file correctly.\n";
    print "\n";
    print "  -out <out_interval>\n";
    print "    <out_interval> = Interval to aggregate to.  Must be one of the following\n";
    print "    strings (case-insensitive):\n";
    print "      hourly  = output file has format:\n";
    print "                YEAR MONTH DAY HOUR DATA DATA DATA ...\n";
    print "      daily   = output file has format:\n";
    print "                YEAR MONTH DAY DATA DATA DATA ...\n";
    print "      monthly = output file has format:\n";
    print "                YEAR MONTH DATA DATA DATA ...\n";
    print "      yearly (or annual) = output file has format:\n";
    print "                YEAR DATA DATA DATA ...\n";
    print "    Default: monthly.\n";
    print "\n";
    print "  [ [ -<agg_method> <field>[,<field>,...] ] [ -<agg_method> <field>[,<field>,...] ] ... ]\n";
    print "    List of aggregation methods and fields to be processed with those methods,\n";
    print "    where\n";
    print "      <agg_method> = one of the following strings (case-insensitive):\n";
    print "        avg = average the field over the aggregation interval.\n";
    print "        sum = sum the field over the aggregation interval.\n";
    print "        beg = use field\'s value at beginning of aggregation interval.\n";
    print "        end = use field\'s value at end of aggregation interval.\n";
    print "        min = use field\'s minimum value within aggregation interval.\n";
    print "        max = use field\'s maximum value within aggregation interval.\n";
    print "      <field>[,<field>,...] = comma-separated list of one or more fields\n";
    print "        to be processed via the specified aggregation method.  Field indexing\n";
    print "        starts with 1 and all fields are included (i.e. the YEAR field is field 1).\n";
    print "    Default: All fields not specified will be aggregated using \"avg\".\n";
    print "    NOTE: If you list a field more than once, the FINAL aggregation method you\n";
    print "          list it with will be the one that is used.\n";
    print "\n";
    print "  Examples:\n";
    print "    To aggregate a file of daily records to monthly values, where the 5th\n";
    print "    field is to be summed instead of averaged:\n";
    print "      $0 -i in_file -o out_file -sum 5\n";
    print "\n";
    print "    To aggregate a file of hourly records to monthly values, where fields\n";
    print "    5 and 6 are to be summed instead of averaged:\n";
    print "      $0 -i in_file -o out_file -in hourly -sum 5,6\n";
    print "\n";
    print "    To aggregate a file of 3-hourly records to yearly values, where fields\n";
    print "    5 and 6 are to be summed instead of averaged:\n";
    print "      $0 -i in_file -o out_file -in hourly -out yearly -sum 5,6\n";
    print "\n";
  }

}
