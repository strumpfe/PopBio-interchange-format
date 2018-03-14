#!/usr/bin/perl
#########################################################################################################
##  convert to interchange format                                                                      ##
##  South Walton County Mosquito Control District (SWCMCD)                                             ##
##                                                                                                     ##
##  2018 Daniel Lawson (daniel.lawson@imperial.ac.uk)                                                  ##
##                                                                                                     ##
#########################################################################################################

use Text::CSV;
use Getopt::Long;
use IO::Handle;
use DateTime;
use strict;

my $verbose;
my $moreverbose;
my $file;
my $prefix;
my $summary;
my $test;
my %species;
my %collection_history;
my %ID_count;
my @columns;		# Column definitions
my $collection = 0;
my $sample     = 0;
my $total_count;
my $line_count;

##
## Hardcoded
##
my $start_column = 14;
my $total_column = 44;


#---------------------------------------------------------#
GetOptions (
    # Misc
    "verbose"       => \$verbose,
    "verboser"      => \$moreverbose,
    # Data
    "file=s"        => \$file,
    "colstart|cs=s"	=> \$start_column,
    "coltotal|ct=s"	=> \$total_column,
    # Output
    "prefix=s"	=> \$prefix,
    # Mode
    "summary" 	=> \$summary,
    "test"	=> \$test,
    );
#---------------------------------------------------------#

unless ( $file ) {
	print " No input file supplied\n";
	exit(0);
}

unless ( $prefix ) {
	print " No identifier prefix supplied\n";
	exit(0);
}

##
## define the nomenclature and output file name
##

my $project           = "${prefix}_surveillance";
my $prefix_collection = "${prefix}_collection";
my $prefix_sample     = "${prefix}_sample";
my $outfile           = "${prefix}.pop";

if ( ($verbose) || ($summary) ) {
	print "// Project:           $project\n";
	print "// File:              $file\n";
	print "// Species columns:   Start=$start_column, collection total=$total_column\n";
}

open (OUTPUT, "> $outfile");
print OUTPUT "collection,sample,collection_date_start,collection_date_end,trap_id,GPS_lat,GPS_lon,location_description,trap_type,attractant,trap_number,trap_duration,species,identification_method,sex,sample_count,assay_type,assay_method,assay_result,assay_units\n";

my $csv = Text::CSV->new({ binary => 1 } );              # create a new object

##
## Main loop through input csv file
##

open (my $data , "< $file") || "Failed to open file: $file\n";
while (my $line = <$data>) {
	chomp $line;

	if ($csv->parse($line)) {
		my @f = $csv->fields();

		# record the number of rows (collections)
		$line_count++;

		# Get column definitions from the header line
		if  ( $line_count == 1 ) {
			@columns = @f;
			next;
		}

		###################################################################################################################################
		##
		## All the main variable in one place......
		## Need to bind columns to data fields and/or set default values
		##
		###################################################################################################################################

		my $collection_date_start            = $f[12];		# Collection end data will be calculated using $collection_timespan
		my $collection_timespan              = 1;		# Number of days trap deployed
		my $collection_GPS_lat               = $f[4];		# Latitude in decimal degrees
		my $collection_GPS_lon               = $f[5];		# Longitude in decimal degrees
		my $collection_location_description  = $f[6];		# Description of collection location
		my $collection_trap_ID               = $f[0];		# Trap ID
		my $collection_trap_type             = $f[1];		# Trap type/protocol
		my $collection_trap_attractant       = $f[2];		# Attractant used in trap
		my $collection_trap_number           = 1;		# Number of traps
		my $collection_trap_duration         = 1;		# Number of days trap operating
		my $collection_species_id_protocol   = "SPECIES_MORPHO";	# Species identification method
		my $collection_devstage              = "adult";		# Developmental stage
		my $collection_sex                   = "female";		# Sex

		# Track collections by portmanteau of time/location/trap_type
		my $id = "${collection_trap_ID}_${collection_date_start}_${collection_GPS_lat}_${collection_GPS_lon}";

		####################################################################################################################################

		print "ID :: $id\n" if ( $moreverbose );
		print  "// WARNING - seen this ID before, assign to existing collection \n" if ( ($collection_history{$id} ne "" ) && ($verbose) );

		if ( $test ) {
			print "$id\n" if ($verbose);
			$ID_count{$id}++;
		}
		else {
			## Species column counts
			for ( my $i = $start_column; $i <= $total_column; $i++ ) {

				if ( $i == $total_column )  {
					if ( $f[$i] == 0 ) {
						print "// No mosquitoes collected for $id $collection_date_start $i count=$f[$i]\n" if ($moreverbose);
						# Housekeeping for ID ordinals
						$collection++;
						$collection_history{$id} = $collection;
						$sample++;
						# Determine 1 day period
						my ($year,$month,$day) = $collection_date_start =~ (/(\d{4})-(\d{2})-(\d{2})/);
						my $dt1 = DateTime->new( year => $year, month => $month, day => $day );
						my $dt2 = $dt1->clone->add( days => ¢collection_timespan);
						# Print collection details
						printf OUTPUT ("${prefix_collection}_%05d,${prefix_sample}_%05d,", $collection, $sample);
						print OUTPUT $dt1->ymd . "," . $dt2->ymd . ",";
						print OUTPUT "$collection_trap_ID,";
						print OUTPUT "$collection_GPS_lat,$collection_GPS_lon,\"$collection_location_description\",";
						print OUTPUT "$collection_trap_type,\"$collection_trap_attractant\",";
						print OUTPUT "$collection_trap_number,$collection_trap_duration,";
						print OUTPUT "BLANK,";
						print OUTPUT "$collection_species_id_protocol,$collection_devstage,$collection_sex,$f[$i],\n";
						print OUTPUT ",,,\n";
						}
					next;
				}

				# discard zero value cells
				if ( $f[$i] eq "" ) {
					print "// Confirmed absence, zero sample size, for $id $f[12] $i count=$f[$i]\n" if ($moreverbose);
				}
				else {
					# Stable ID for collection
					if ( $collection_history{$id} ne "") {
						$sample++;
						printf OUTPUT ("${prefix_collection}_%05d,${prefix_sample}_%05d,", $collection_history{$id}, $sample);
					}
					else {
						$collection++;
						$collection_history{$id} = $collection;
						$sample++;
						printf OUTPUT ("${prefix_collection}_%05d,${prefix_sample}_%05d,", $collection, $sample);
					}

					# [WHEN] Collection times
					my ($year,$month,$day) = $collection_date_start =~ (/(\d{4})-(\d{2})-(\d{2})/);
					my $dt1 = DateTime->new( year => $year, month => $month, day => $day );
					my $dt2 = $dt1->clone->add( days => $collection_timespan);

					print OUTPUT $dt1->ymd . "," . $dt2->ymd . ",";
					print OUTPUT "$collection_trap_ID,";
					print OUTPUT "$collection_GPS_lat,$collection_GPS_lon,\"$collection_location_description\",";
					print OUTPUT "$collection_trap_type,\"$collection_trap_attractant\",";
					print OUTPUT "$collection_trap_number,$collection_trap_duration,";
					print OUTPUT "$columns[$i],";
					print OUTPUT "$collection_species_id_protocol,$collection_devstage,$collection_sex,$f[$i],";
					print OUTPUT ",,,\n";

					# add sample_count to total_count
					$total_count = $total_count + $f[$i];
					$species{$columns[$i]} = $species{$columns[$i]} + $f[$i];
				}
			}
		}
	}
	else {
	 	print "// parsing line failed [$line]\n" if ($verbose);
	}
}
close FILE;

close OUTPUT;

if ( ($verbose) || ($summary) ) {
	print "\n// Parsed $line_count rows and created " . $sample  . " samples from " . $collection . " collections\n";
	print "// Total number of mosquitoes: $total_count from " . keys(%species) ." species\n\n";
	foreach my $i (sort keys %species) {
		printf ("// %-32s %4d\n", $i, $species{$i});
	}
}

if ( $test ) {
	my $errors;
	foreach my $i (sort keys %ID_count) {
		if ( $ID_count{$i} > 1 ) {
			printf ("// ID: %-48s is not unique (n=$ID_count{$i})\n", $i);
			$errors++;
		}
	}
	if ( ! defined $errors ) {
		print "// All IDs are unique, good to process this file\n";
	}
}
