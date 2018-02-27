#!/usr/bin/perl
#########################################################################################################
##  convert to interchange format                                                                      ##
##  South Walton County Mosquito District (SWCMCD)                                                     ##
##                                                                                                     ##
##  2018 Daniel Lawson (daniel.lawson@imperial.ac.uk)                                                  ##
##                                                                                                     ##
#########################################################################################################

use Text::CSV;
use Getopt::Long;
use IO::Handle;
use strict;

my $verbose;
my $moreverbose;
my $file;
my $summary;
my $prefix;

#---------------------------------------------------------#
GetOptions (
    # Misc
    "verbose"       => \$verbose,
    "verboser"      => \$moreverbose,
#    "help"          => \$help,
    # Data
    "file=s"        => \$file,
    # Output
    "prefix=s"	=> \$prefix,
    # Mode
    "summary" 	=> \$summary,
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
my $output            = "${prefix}.pop";

my $writesamples     = 0;
my $writecollections = 0;

# Define initial values for collection and samples
my $collection = 1;
my $sample     = 1;

my $collection_id;
my $sample_id;
my %collection_history;

my $rows_discarded = 0;
my $rows_processed = 0;
my $rows_in_total  = 0;
my $total_count    = 0;

my %species_seen;
my @species;

# We are parsing csv files
my $csv = Text::CSV->new({ binary => 1 });              # create a new object

## open output file
open (OUTPUT, "> $output");
print OUTPUT "collection,sample,collection_date_start,collection_date_end,trap_id,GPS_lat,GPS_lon,location_description,trap_type,attractant,trap_number,trap_duration,species,identification_method,developmental_stage,sex,sample_count\n";


## read from input file
## "",
## $f[1]	"variable",
## $f[2]	"Collection site county",
## $f[3]	"Protocol REF",
## $f[4]	"Collection site latitude",
## $f[5]	"Attractant(s) Used",
## $f[6]	"Collection site longitude",
## $f[7]	"Collection site city",
## $f[8]	"Collection site province",
## $f[9]	"Date",
## $f[10]	"Collection site locality",
## $f[11]	"collection",
## $f[12]	"Sample size",
## $f[13]	"Characteristics",
## $f[14]	"sample",
## $f[15]	"VBsp",
## $f[16]	"MosquitoNetList"

open (my $data , "< $file");
while (my $line = <$data>) {
	chomp $line;
	$rows_in_total++;

	# ignore header lines (different for each file?)
	next if  ( $line =~ /Date of collection/ );

	if ($csv->parse($line)) {
		my @f = $csv->fields();

		# Track collections by portmanteau of date/location/trap_type/lat/lon
		my $id = "$f[9]_$f[7]_$f[3]_$f[4]_$f[6]";

		if ( $verbose ) {
			print "ID :: $id\t";
			if (  $collection_history{$id} ne "" ) {
				print  "WARNING - seen this ID before\n"
			}
			else {
				print "\n";
			}
		}

		# discard zero value cells
		if ( $f[12] == 0 ) {
			print "// Confirmed absence, zero sample size, ignore\n" if ($moreverbose);
			$rows_discarded++;
			next;
		}

		# Stable ID for collection
		if ( $collection_history{$id} ne "") {
			# Assert collection and sample identifiers
			$collection_id = sprintf ("${prefix_collection}_%04d", $collection_history{$id});
			$sample_id     = sprintf ("${prefix_sample}_%04d", $sample);

			print OUTPUT "$collection_id,$sample_id,";
			$sample++;
		}
		else {
			# Assert collection and sample identifiers
			$collection_id = sprintf ("${prefix_collection}_%04d", $collection);
			$sample_id     = sprintf ("${prefix_sample}_%04d", $sample);

			print OUTPUT "$collection_id,$sample_id,";
			$collection_history{$id} = $collection;
			$collection++;
			$sample++;
		}

		# [WHEN] Collection times
		# [collection_date_start] [collection_date_end]
		print OUTPUT "$f[9],$f[9],";

		# Trap ID		[No trap IDs for SWCMCD]
		# [trap ID]
		print OUTPUT ",";

		# [WHERE] Collection location
		# [latitude] [longitude] [location]
		print OUTPUT "$f[4],$f[6],\"$f[10]\",";

		# [HOW] Trap details	[Default single trap for single night collection]
		# [Trap type] [Attractant][Trap number] [Trap duration]
		print OUTPUT "$f[3],$f[5],1,1,";

		# [WHAT] Species and sample_count	Assert default values for [SPECIES_MORPHO] [adult]
		# [Species] [Species identification protocol] [Developmental stage] [Sex] [Sample size]
		print OUTPUT "$f[1],SPECIES_MORPHO,adult,";
		if    ( $f[13] eq "female count" ) { print OUTPUT "female,"; }
		elsif ( $f[13] eq "male count" )   { print OUTPUT "male,";   }
		print OUTPUT "$f[12]\n";

		# add sample_count to total_count
		$rows_processed++;
		$total_count = $total_count + $f[12];
		push @{ $species_seen{$collection_id} }, $f[1];
		push @species, $f[1];

	}
	else {
		print "// parsing line failed [$line]\n" if ($verbose);
	}
}
close FILE;

##
## Summary of parsing the raw file
##

print "\n// Convert to interchange format : $project\n\n";
print "Rows in total:  $rows_in_total\n";
print "Rows processed: $rows_processed\n";
print "Rows discarded: $rows_discarded\n";
if ( ($rows_processed + $rows_discarded) != $rows_in_total ) { print "// WARNING - disparity in processing of input file rows\n"; }
print "\n";
print "No. of collections: " . ($collection -1) . "\n";
print "No. of samples:     " . ($sample-1)      . "\n";
print "Total count:        $total_count\n\n";


if ( defined $summary ) {
	foreach my $i (sort keys %species_seen) {

		my @full_species = uniq( @species );

		# print array contents (redundant list of species seen)
#		printf ("%-36s :: [", $i);
#		print scalar @{ $species_seen{$i} };
#		print "] :: ";
#		print join ',', @{ $species_seen{$i} };
#		print "\n";

		# print distinct array contents (non-redundant list of species_seen)
		my @distinct = uniq( @{ $species_seen{$i} } );
		printf ("%-36s :: [", $i);
		print scalar @distinct;
		print "] :: ";
		print join ',', @distinct;
		print "\n";

		# print list of AWOL species based on the full enumerated list in %species
		printf ("%-36s :: [", $i);
		my @missing = grep { ! ({ $_, 0 } ~~ @distinct) } @full_species;
		print scalar @missing;
		print "] :: ";
		print join ", ", @missing;
		print "\n";

		print "\n";
	}
}
close OUTPUT;

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
