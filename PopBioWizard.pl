#!/usr/bin/perl
##
## PopBioWizard.pl
## daniel.lawson@imperial.ac.uk
## 2018-01-17 v0.1

use strict;
use Getopt::Long;
use IO::Handle;
use Text::ParseWords;
use Text::CSV;
use Data::Dumper;

# Hash name
our %ISA;
our %meta;
our %species;
our @species;
our %ontology_lookup;
our %species_seen;		# List of reported species for each collection
our %collections2add;	# Hash of collection IDs for which confirmed absence needs to be added
our %collection_meta;
our %geoname;

# Options
my $verbose;             	# Verbose messages
my $moreverbose;         	# More verbosity
my $debug;               	# Write out data parsing debug information (very verbose)
my $help;                	# Help documentation via POD
my $local;               	# Read from local hash (.dat format)
my $file;                	# Read from local data file (.pop format)
my $configfile;          	# Local configuration file with project metadata (akin to the i_investigation sheet)
my $add_zeros;	     	# Programmatically add zero sized sample for each collection (impute species from full list)
my $bg_counter;         # Use BG-Counter term for Protocol REF
my $max_species_num;	# Total number of species reported/tested for
my $max_sample_id;  	# Maximum ordinal in sample ID, needed for when adding verified absences (zero sample size)
my $s_sample;
my $a_collection;
my $a_species;
my $a_virus;
#---------------------------------------------------------#

#---------------------------------------------------------#
GetOptions (
    # Misc
    "verbose"       => \$verbose,
    "verboser"      => \$moreverbose,
    "help"          => \$help,
    # Data
    "local"         => \$local,
    "file=s"        => \$file,
    "config=s"      => \$configfile,
    # Process
    "zeros"         => \$add_zeros,
    "bg-counter"    => \$bg_counter,
    # Output
    "sample"        => \$s_sample,
    "collection"    => \$a_collection,
    "species"       => \$a_species,
    "virus"	=> \$a_virus,
    );
#---------------------------------------------------------#

my $dir         = $ENV{"HOME"} . "/Dropbox/VectorBase/PopBio/DATA_FILES";    	# Define local directory for data files
my ($data_file) = $file =~ (/(\S+)\.pop/); $data_file = $data_file . ".dat";	# Define output dat file for local copy of hash

##
## Get project metadata from configuration file
##

&get_project_metadata;

print "// PopBioWizard run ". gmtime( time()) ."\n";
print "//  configuration from $configfile\n";
print "//  data from file $file\n" if ( $file );
print "//  data from hash $dir/$data_file\n" if ( $local );
print "\n";
print "// Study identifier = $meta{Study_identifier}\n";
print "// No. of species tested : $max_species_num\n";
#print "// Species list : [" . ( join ',', @species) . "]\n";
foreach my $i (sort keys %species) {
	printf ("//  %-30s $species{$i}\n", $i);
	}
my $ontoterms = keys %ontology_lookup;
print "\n// No. of ontology terms : $ontoterms\n";
foreach my $i (sort keys %ontology_lookup) {
	printf ("//  %-30s $ontology_lookup{$i}\n", $i);
	}
print "\n";



##
## Read collection information
## either from .dat hash table or local file
##

# retrieve data from local file (.pop format)
if ($file) {
  print "// Using data from file :: $file\n\n" if ( $moreverbose );
  &get_data_from_file;
}
# retrieve data from local hash structure
elsif ($local) {
  # retrieve the hash from the file.
  print "// Using data from hash :: $dir/$data_file\n\n" if ( $moreverbose );
  &get_data_from_local_hash;
}
else {
  print "# $0 " . gmtime( time()) ."\n\n";
  print "No data source declared. Aborting.\n";
  print "Use one of the following options:\n";
  print " '-local' option to read from local hash structure\n";
  print " '-file <filename>' to read from local file (.pop format)\n\n";
  print " '-help' for full documentation\n\n";
  exit(0);
}


##
## impute confirmed absence of mosquitoes
##

if ( $add_zeros ) {

	# Loop through hash to collate list of seen species
	foreach my $i ( keys %ISA ) {

		foreach my $j ( @{ $ISA{$i}{Species} } ) {
			print "// $j collected on $ISA{$i}{ColStart} $ISA{$i}{ColID} $ISA{$i}{SamID}\n" if ( $moreverbose );
			push @{ $species_seen{$ISA{$i}{ColID}} }, $j;
			$collections2add{$ISA{$i}{ColID}} = 1;		# Mark collection for addition
		}
	}

	# Loop through collections that need augmenting
	foreach my $i ( sort keys %collections2add ) {

		my @distinct = uniq( @{ $species_seen{$i} } );
		my @missing  = grep { ! ({ $_, 0 } ~~ @distinct) } @species;

		my $number_spp = scalar @distinct;
		print "// Collection $i ($number_spp) :: [" . (join ', ', @distinct) . "]\n" if ( $moreverbose );

		# Check whether this collection is complete (has all species been found?)
		# If true then we don't need to make a new sample for the collection to store confirmed absences
		next if ($number_spp == $max_species_num);

		# Discard if this is already a blank entry (no mosquitoes collected)
		if ( $distinct[0] eq "BLANK") {
			print "// Not making a zero entry as this is a BLANK collection $i\n\n" if ( $verbose );
			next;
		}

		if ( ($missing[0] =~ /genus/ ) or ($missing[0] eq "Culicinae") or ($missing[0] eq "Culicidae") ) {
			print "// Not making a zero entry as missing species is a generic term $i\n\n" if ( $verbose );
			next;
		}

		# Make a new sample for confirmed absences
		print "// Need to add confirmed absence for collection $i, " . ($max_species_num - $number_spp ) . "\n" if ( $verbose );

		print "// Collection $i ($number_spp) :: [" . (join ', ', @missing) . "]\n" if ( $moreverbose );
		$max_sample_id++;

		my $new_sample_id = sprintf ("$meta{Sample_nomenclature}_%.5d", $max_sample_id);
		print "// Create new sample \"$new_sample_id\" for collection $i\n\n" if ( $verbose );

		# Collection IDs
		$ISA{$new_sample_id}{ColID}           = $i;
		$ISA{$new_sample_id}{SamID}           = $new_sample_id;
		$ISA{$new_sample_id}{SamDesc}         = "Record of absence of some species of mosquito";   # sample description
		# Collection date
		$ISA{$new_sample_id}{ColStart}        = $collection_meta{$i}{ColStart};  # collection start date
		$ISA{$new_sample_id}{ColEnd}          = $collection_meta{$i}{ColEnd};    # collection end date
		# Collection site location
		$ISA{$new_sample_id}{ColLat}          = $collection_meta{$i}{ColLat};   # collection site latitude
		$ISA{$new_sample_id}{ColLon}          = $collection_meta{$i}{ColLon};   # collection site longitude
		$ISA{$new_sample_id}{ColDesc}         = $collection_meta{$i}{ColDesc};  # collection site description
		# Collection trap metadata
		$ISA{$new_sample_id}{TrapType}        = $collection_meta{$i}{TrapType};        # trap type
		$ISA{$new_sample_id}{TrapAttractant}  = $collection_meta{$i}{TrapAttractant};  # trap attractant
		$ISA{$new_sample_id}{TrapDuration}    = $collection_meta{$i}{TrapDuration};    # Duration of trap deployment
		$ISA{$new_sample_id}{TrapQuantity}    = $collection_meta{$i}{TrapQuantity};    # No. of traps deployed
		# Collected material metadate
		$ISA{$new_sample_id}{SamQuantity}     = 0;   # No. of animals collected
		foreach my $j (@missing) {
			next if ( $j eq "Culicidae" );	# Ignore generic species, don't make confirmed zero sample size for generic terms
			next if ( $j eq "Culicinae" );
			next if ( $j =~ /genus/ );

			push @{ $ISA{$new_sample_id}{Species} },  $j;   #
		}
		$ISA{$new_sample_id}{SpeciesProc}          = "SPECIES_MORPHO";   	# species identification method
		$ISA{$new_sample_id}{SamSex}               = "female";   		# sex
		$ISA{$new_sample_id}{SamStage}             = "adult";    		# developmental stage
	}
}

##
## s_sample sheet output
##

if ( $s_sample ) {

	open (OUTPUT, "> ./s_sample.csv");

	# Header
	print OUTPUT "Source Name,Sample Name,Description,Comment [comment],Material Type,Term Source Ref,Term Accession Number,Comment [age],Characteristics [sex (EFO:0000695)],Term Source Ref,Term Accession Number,Characteristics [developmental stage (EFO:0000399)],Term Source Ref,Term Accession Number,Characteristics [sample size (VBcv:0000983)]\n";

	foreach my $i (keys %ISA) {
		printf OUTPUT ("$meta{Study_identifier},$ISA{$i}{SamID},$ISA{$i}{SamDesc},,", $i );
		# Material type, default "pool"
		my $look = "pool";
		my ($ontology,$accession) = $ontology_lookup{$look} =~ (/(\S+?)\:(\S+)/);
		print OUTPUT "$look,$ontology,$accession,,";
		# Sex
		my $look = $ISA{$i}{SamSex};
		my ($ontology,$accession) = $ontology_lookup{$look} =~ (/(\S+?)\:(\S+)/);
		print OUTPUT "$look,$ontology,$accession,";
		# Developmental stage
		my $look = $ISA{$i}{SamStage};
		my ($ontology,$accession) = $ontology_lookup{$look} =~ (/(\S+?)\:(\S+)/);
		print OUTPUT "$look,$ontology,$accession,";
		# Sample size
		print OUTPUT "$ISA{$i}{SamQuantity}\n";

	}
	close OUTPUT;
}


##
## a_collection sheet output
##

# [DL] Currently has no provision for adding collection location, GAZ term etc.

if ( $a_collection ) {

	open (OUTPUT, "> ./a_collection.csv");

	# Header
	print OUTPUT "Sample Name,Assay Name,Description,Protocol REF,Performer,Date,Characteristics [sampling time (EFO:0000689)],Characteristics [Temperature at time of collection (EFO:0001702)],Unit,Term Source Ref,Term Accession Number,Comment [household ID],Characteristics [Collection site (VBcv:0000831)],Term Source Ref,Term Accession Number,Characteristics [Collection site latitude (VBcv:0000817)],Characteristics [Collection site longitude (VBcv:0000816)],Characteristics [Collection site altitude (VBcv:0000832)],Comment [collection site coordinates],Characteristics [Collection site location (VBcv:0000698)],Characteristics [Collection site village (VBcv:0000829)],Characteristics [Collection site locality (VBcv:0000697)],Characteristics [Collection site suburb (VBcv:0000845)],Characteristics [Collection site city (VBcv:0000844)],Characteristics [Collection site county (VBcv:0000828)],Characteristics [Collection site district (VBcv:0000699)],Characteristics [Collection site province (VBcv:0000700)],Characteristics [Collection site country (VBcv:0000701)]\n";

	foreach my $i (keys %ISA) {
		printf OUTPUT ("$ISA{$i}{SamID},$ISA{$i}{ColID},\"$ISA{$i}{ColDesc}\",COLLECT_$ISA{$i}{TrapType},,$ISA{$i}{ColStart},,,,,,,,GAZ,,$ISA{$i}{ColLat},$ISA{$i}{ColLon},,IA,,,,,,,,,\n", $i );
	}
	close OUTPUT;
}


##
## a_species sheet output
##


if ( $a_species ) {

	open (OUTPUT, "> ./a_species.csv");

	print OUTPUT "Sample Name,Assay Name,Description,Protocol REF,Performer,Date,Characteristics [species assay result (VBcv:0000961)],Term Source Ref,Term Accession Number\n";

	my $species_proc;  # species identification method

	if ( $bg_counter ) {
		$species_proc = "SPECIES";
	}
	else {
		$species_proc = "SPECIES_MORPHO";
	}

	foreach my $i (keys %ISA) {
		printf OUTPUT ("$ISA{$i}{SamID},$ISA{$i}{SamID}.spp,,$species_proc,,$ISA{$i}{ColStart},", $i );
		my ($sp_species,$sp_onto,$sp_acc);
		foreach my $j ( @{ $ISA{$i}{Species} } ) {
			my ($onto,$acc) = $species{$j} =~ ( /^(\S+?)\:(\S+)/ );
			if ( $acc eq "" ) { print "// WARNING: No ontology term for $j $i [$ISA{$i}{ColID} : $ISA{$i}{SamID}]\n"; }
			$sp_species = $sp_species . $j . ";";
			$sp_onto    = $sp_onto . $onto . ";";
			$sp_acc     = $sp_acc . $acc   . ";";
		}
		chop ($sp_species);
		chop ($sp_onto);
		chop ($sp_acc);
		print OUTPUT "$sp_species,$sp_onto,$sp_acc\n";
		undef ($sp_species);
		undef ($sp_onto);
		undef ($sp_acc);
	}
	close OUTPUT;
}


##
## a_ASSAY sheet output
##

if ( $a_virus ) {
	open (OUTPUT_A, "> ./a_virus.csv");
	open (OUTPUT_P, "> ./p_virus.csv");

	# Header
	print OUTPUT_A "Sample Name,Assay Name,Description,Protocol REF,Raw Data File\n";
	print OUTPUT_P "Assay Name,Phenotype Name,Observable,Term Source Ref,Term Accession Number,Attribute,Term Source Ref,Term Accession Number,Value,Term Source Ref,Term Accession Number\n";

	foreach my $i (keys %ISA) {

		foreach my $j (@{ $ISA{$i}{Assay} }) {
			my ($assay_type,$assay_proc,$assay_result,$assay_units) = split /;/, $j;
			print "// $ISA{$i}{SamID} [$assay_type,$assay_proc,$assay_result,$assay_units] $j\n" if ( $moreverbose );


			printf OUTPUT_A ("$ISA{$i}{SamID},$ISA{$i}{SamID}.${assay_type},,${assay_type},p_virus.txt\n", $i );

			if ( $assay_result eq "Positive") {
				printf OUTPUT_P ("$ISA{$i}{SamID}.${assay_type},\"${assay_type} infected\",\"arthropod infection status\",VSMO,0000009,${assay_type},VSMO,0000535,present,PATO,0000467\n");

			}
			elsif ( $assay_result eq "Negative") {
 				printf OUTPUT_P ("$ISA{$i}{SamID}.${assay_type},\"${assay_type} infection not detected\",\"arthropod infection status\",VSMO,0000009,${assay_type},VSMO,0000882,absent,PATO,0000462\n");
			}
		}
	}
	close OUTPUT_A;
	close OUTPUT_P;
}


##
## Write out perl hash of data to local disk
##

&write_data_to_local_file;

exit(0);

#-----------------------------------------------------------------------------------------------#


##
## get_data_from_hash
##

sub get_data_from_local_hash {
	open (FH, "< $dir/$data_file") or die "$dir/$data_file : $!";
	undef $/;
	my $data = <FH>;
	eval $data;
	die if $@;
	$/ = "\n";
	close FH;
}

##
## get_data_from_file
##

# f[ 0] = collection ID
# f[ 1] = sample ID
# f[ 2] = collection start date (ISO 8601)
# f[ 3] = collection end date (ISO 8601)
# f[ 4] = trap ID
# f[ 5] = latitude, GPS lat (decimal)
# f[ 6] = longitude, GPS lon (decimal)
# f[ 7] = location description (text string, may contain commas)
# f[ 8] = trap type
# f[ 9] = attractant
# f[10] = duration of deployment (number of nights, interger
# f[11] = number of traps deployed (interger)
# f[12] = species collected (usually binomial name)
# f[13] = species identification method (text string)
# f[14] = collected animals developmental stage (larva,pupa,adult)
# f[15] = collected animals sex (female|male|mixed)
# f[16] = number of animals (interger)
# f[17] = assay type/name
# f[18] = assay protocol/method
# f[19] = assay result
# f[20] = units of measurement (where needed)

##
## primary key is the sample ID
##

sub get_data_from_file {
	my $csv = Text::CSV->new( { binary => 1 } );              	# create a new object
	my $row_count;					# count number of rows processed for reporting
	my $row_parsed;					# count number of rows parsed for reporting

	open (my $data, "< $file") or die "Can't open file '$file' : $!\n";
	while (my $line = <$data>) {

		$row_count++;
		chomp $line;
		next if ( $line =~ /^\/\// );  # discard header line
		next if ( $line =~ /^collection\,sample/ );

		if ($csv->parse($line)) {
			my @f = $csv->fields();

			$ISA{$f[1]}{ColID}                = $f[0];   # collection ID
			$ISA{$f[1]}{SamID}                = $f[1];   # sample ID
			# Collection date
			$ISA{$f[1]}{ColStart}             = $f[2];   # collection start date
			$ISA{$f[1]}{ColEnd}               = $f[3];   # collection end date
			# Collection site location
			$ISA{$f[1]}{ColLat}               = $f[5];   # collection site latitude
			$ISA{$f[1]}{ColLon}               = $f[6];   # collection site longitude
			$ISA{$f[1]}{ColDesc}              = $f[7];   # collection site description
			# Collection trap metadata
			$ISA{$f[1]}{TrapType}             = $f[8];   # trap type

			# If trap_type is AWOL set to MIRO:30000045 "catch of live specimens"
			if ( $f[8] eq "" ) {
				$ISA{$f[1]}{TrapType} = "LIVE";
			}

			$ISA{$f[1]}{Attractant}           = $f[9];   # trap type

			$ISA{$f[1]}{TrapDuration}         = $f[10];   # Duration of trap deployment
			$ISA{$f[1]}{TrapQuantity}         = $f[11];   # No. of traps deployed
			# Collected material metadate
			$ISA{$f[1]}{SamQuantity}          = $f[16];   # No. of animals collected

			if  ( $f[12] eq "BLANK" ) {
				print "// Assert all species for the Blank collection $f[0] $f[1]\n" if ($verbose);
				$ISA{$f[1]}{SamDesc} = "Record of absence of some species of mosquito";
				foreach my $j (@species) {
					# Ignore generic species assertion, don't make confirmed zero sample size for generic terms
					next if ( $j eq "Culicidae" );
					next if ( $j eq "Culicinae" );
					next if ( $j =~ /genus/ );
					push @{ $ISA{$f[1]}{Species} },  $j;
				}
			}
			else {
				push @{ $ISA{$f[1]}{Species} },     $f[12];   # species
			}

			$ISA{$f[1]}{SpeciesProc}          = $f[13];   # species identification method
			$ISA{$f[1]}{SamStage}             = $f[14];   # developmental stage
			$ISA{$f[1]}{SamSex}               = $f[15];   # sex

			# Assays metadata and result
			push @{ $ISA{$f[1]}{Assay} },      "$f[17];$f[18];$f[19];$f[20];";   # assay type


			# Parallel tracking of collection metadata ()
			#
			# Lazy overwriting of (presumed) consistent collection metadata

			$collection_meta{$f[0]}{ColStart}      = $f[2];    # collection start date
			$collection_meta{$f[0]}{ColEnd}        = $f[3];    # collection end date

			$collection_meta{$f[0]}{ColLat}        = $f[5];    # collection site latitude
			$collection_meta{$f[0]}{ColLon}        = $f[6];    # collection site longitude
			$collection_meta{$f[0]}{ColDesc}       = $f[7];    # collection site description

			$collection_meta{$f[0]}{TrapType}      = $f[8];    # trap type
			# If trap_type is AWOL set to MIRO:30000045 "catch of live specimens"
			if ( $f[8] eq "" ) {
				$collection_meta{$f[0]}{TrapType} = "LIVE";
			}
			$collection_meta{$f[0]}{TrapAttractant} = $f[9];    # trap type
			$collection_meta{$f[0]}{TrapDuration}   = $f[10];   # Duration of trap deployment
			$collection_meta{$f[0]}{TrapQuantity}   = $f[11];   # No. of traps deployed

			# Calculating the maximum sample ID ordinal
			my ( $ord ) = $f[1] =~ ( /sample_(\d{5})/ );
			if ( $ord > $max_sample_id ) {
				$max_sample_id = $ord;
			}

			# Tracking for debug/reporting
			$row_parsed++;
		}
	}

	# Report number of rows parsed from the input file
	if ( $verbose ) {
		print "// Parsed file      : $file\n";
		print "// No. rows in file : $row_count\n";
		print "// No. rows parsed  : $row_parsed\n";
		print "// Max. sample ID   : $max_sample_id\n";
	}
}

##
## get_project_metadata
##

sub get_project_metadata {
	open (my $data, "< $configfile") or die "Can't open file '$configfile' : $!\n";
	while (my $line = <$data>) {
		chomp $line;
		# Study identifier
		if ( ( $line =~ /^Study_identifier\s+\:\s+(\S+.+)$/) ) {
			$meta{Study_identifier} = $1;
		}
		if ( ( $line =~ /^Sample_nomenclature\s+\:\s+(\S+.+)$/) ) {
			$meta{Sample_nomenclature} = $1;
		}
		if ( ( $line =~ /^Collection_nomenclature\s+\:\s+(\S+.+)$/) ) {
			$meta{Collection_nomenclature} = $1;
		}
		# Enumerated list of species
		if ( ( $line =~ /^Study_species\s+\:\s+\'(\S+.+?)\',\s+\'(\S+)\'$/ ) ) {
			$species{$1} = $2;
		}
		# Ontology terms assertions
		if ( ( $line =~ /^Study_ontology\s+\:\s+\'(.+?)\',\s+\'(\S+)\'$/ ) ) {
			$ontology_lookup{$1} = $2;
		}

	}

	# Make an array of the species names as well
	foreach my $j (sort keys %species) {
		 push @species, $j;
          }

	# Store the total number of species reported
	$max_species_num = scalar @species;
}


##
## calculate diff of elements in 2 arrays
##

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}


##
## write_data_to_local_file
##

sub write_data_to_local_file {
	my $str = Data::Dumper->Dump( [\%ISA],[qw/*ISA/] );
	open (FH, "> $dir/$data_file") or die "Can't open $dir/$data_file as output .dat\n";
	print FH $str;
	close FH;
	print "\n// Wrote hash %ISA to '$dir/$data_file'\n\n" if ($verbose);
}

#-----------------------------------------------------------------------------------------------#
