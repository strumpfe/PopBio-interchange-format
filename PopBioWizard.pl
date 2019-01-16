#!/bin/env perl
##
## PopBioWizard.pl
## daniel.lawson@imperial.ac.uk
## 2018-01-17 v0.1

use strict;
use warnings;
use Getopt::Long;
use IO::Handle;
use Text::ParseWords;
use Text::CSV;
use Text::CSV::Hashify 0.11;
use Data::Dumper;
use Tie::Hash::Indexed;
use YAML::XS qw/LoadFile/;
use Bio::Parser::ISATab;

my %ISA;  # main hash to store data read in from $file - main key is sample ID
tie %ISA, 'Tie::Hash::Indexed'; # ordered hash - will return samples in order
my %collection_meta;

# Options
my $verbose;                    # Verbose messages
my $moreverbose;                # More verbosity
my $help;                       # Help documentation via POD
my $file;                       # Read from local data file (.csv or tsv interchange format)
my $configfile;                 # Local configuration file with project metadata (akin to the i_investigation sheet)
my $add_zeros;	     	        # Programmatically add zero sized sample for each collection (impute species from full list)
my $bg_counter;                 # Use BG-Counter term for Protocol REF
my $output_delimiter = ",";     # should be "tab", "TAB", "comma", "COMMA", or ","  (providing an actual TAB on commandline is a hassle)
my $output_directory = "./temp-isa-tab";  # where the output ISA-Tab files should go
my $output_suffix = "txt";      # ISA-Tab output files file suffix - not a commandline option - will be set automatically to 'csv' if needed
my $i_investigation;
my $s_sample;
my $a_collection;
my $a_species;
my $a_virus;
my $all_regular_sheets;

#---------------------------------------------------------#

#---------------------------------------------------------#
GetOptions (
    # Misc
    "verbose"       => \$verbose,
    "verboser"      => \$moreverbose,
    "help"          => \$help,
    # Data
    "file=s"        => \$file,
    "config=s"      => \$configfile,
    # Process
    "zeros"         => \$add_zeros,
    "bg-counter"    => \$bg_counter,
    # Output
    "investigation" => \$i_investigation,  # output the i_investigation.txt sheet (to $output_directory)
    "samples"        => \$s_sample,        # output the s_samples.txt sheet
    "collection"    => \$a_collection,     # and
    "species"       => \$a_species,        # so
    "virus"	    => \$a_virus,          # on
    "isatab"        => \$all_regular_sheets, # shortcut for -investigation -samples -collection --species --output_delimiter TAB
    "output-delimiter|delimiter=s" => \$output_delimiter,
    "output-directory|directory=s" => \$output_directory,
    );
#---------------------------------------------------------#

die "must provide --file <input_filename> and --config <config_filename> options\n" unless ($file && $configfile);

# handle some implied commandline args
$verbose = 1 if ($moreverbose);
if ($all_regular_sheets) {
    ($i_investigation, $s_sample, $a_collection, $a_species, $output_delimiter) = (1,1,1,1,"TAB");
}

# convert $output_delimiter option (any other values will fall through)
$output_delimiter = "\t" if ($output_delimiter =~ /tab/i);
$output_delimiter = "," if ($output_delimiter =~ /comma/i);
my $csv_formatter = Text::CSV_XS->new ({ binary => 1, eol => $/, sep_char => $output_delimiter });
$output_suffix = "csv" if ($output_delimiter eq ",");

##
## Get project metadata from configuration file
##

my $config = LoadFile($configfile); # read from YAML format
die "problem reading config file '$configfile' - is it YAML formatted?\n" unless (keys %$config);

# print Dumper($config); exit;

print "// PopBioWizard run ". gmtime( time()) ."\n";
print "//  configuration from $configfile\n";
print "//  data from file $file\n" if ( $file );
print "\n";
print "// Study identifier = '$config->{study_identifier}'\n";
my @expected_species = keys %{$config->{study_species}};
my $max_species_num = scalar @expected_species;

print "// No. of species tested : $max_species_num\n";

foreach my $i (@expected_species) {
    printf ("//  %-30s $config->{study_species}{$i}\n", $i);
}
my $ontoterms = keys %{$config->{study_terms}};
print "\n// No. of ontology terms : $ontoterms\n";
foreach my $i (sort keys %{$config->{study_terms}}) {
    printf ("//  %-30s $config->{study_terms}{$i}\n", $i);
}
print "\n";


##
## handle i_investigation sheet (only needs data from config file)
##

if ($i_investigation) {
    print "// Writing $output_directory/i_investigation.txt sheet\n" if ($verbose);
    my $isa_parser = Bio::Parser::ISATab->new(directory => $output_directory);
    $config->{study_file_name} = "s_samples.$output_suffix";

    # fill in STUDY DESIGN DESCRIPTORS section if not provided already in config file
    unless ($config->{study_designs}) {
	$config->{study_designs} = [
	    { study_design_type => 'observational design',
	      study_design_type_term_source_ref => 'EFO',
	      study_design_type_term_accession_number => '0000629' },
	    { study_design_type => 'strain or line design',
	      study_design_type_term_source_ref => 'EFO',
	      study_design_type_term_accession_number => '0001754' },
	    ];
    }

    # fill in STUDY ASSAYS section if not provided already in config file
    unless ($config->{study_assays}) {
	push @{$config->{study_assays}}, { study_assay_measurement_type => 'field collection',
					   study_assay_file_name => "a_collection.$output_suffix" } if ($a_collection);
	push @{$config->{study_assays}}, { study_assay_measurement_type => 'species identification assay',
					   study_assay_file_name => "a_species.$output_suffix" } if ($a_species);
	push @{$config->{study_assays}}, { study_assay_measurement_type => 'phenotype assay',
					   study_assay_file_name => "a_virus.$output_suffix" } if ($a_virus);
    }


    # now write the i_investigation sheet using data in $config (any additional non-ISA-Tab data will be ignored)
    $isa_parser->write( { ontologies => [], studies => [ $config ] } );
    print "// Done investigation sheet\n" if ($moreverbose);
}


##
## Read collection information
##

# retrieve data from local file (.pop format)
if ($file) {
    print "// Using data from file :: $file\n\n" if ( $moreverbose );
    &get_data_from_file;
} else {
    print "# $0 " . gmtime( time()) ."\n\n";
    print "No data source declared. Aborting.\n";
    print "Use one of the following options:\n";
    print " '-file <filename>' to read from local file (.csv or .tsv interchange format)\n\n";
    print " '-help' for full documentation\n\n";
    exit(0);
}


##
## impute confirmed absence of mosquitoes
##

if ( $add_zeros ) {

    my %species_seen;		# List of reported species for each collection
    my %collections2add;	# Hash of collection IDs for which confirmed absence needs to be added
    my %zero_sample_count;      # counter for the ordinal added to zero sample IDs: <collection_ID>_zero_sample_NNN

    # Loop through hash to collate list of seen species
    foreach my $i ( keys %ISA ) {

	foreach my $j ( @{ $ISA{$i}{Species} } ) {
	    print "// $j collected on $ISA{$i}{collection_start_date} $ISA{$i}{collection_ID} $ISA{$i}{sample_ID}\n" if ( $moreverbose );
	    push @{ $species_seen{$ISA{$i}{collection_ID}} }, $j;
	    $collections2add{$ISA{$i}{collection_ID}} = 1;		# Mark collection for addition
	}
    }

    # Loop through collections that need augmenting
    foreach my $i ( sort keys %collections2add ) {

	my @distinct = uniq( @{ $species_seen{$i} } );
	my @missing  = grep { ! ({ $_, 0 } ~~ @distinct) } @expected_species;

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

	my $new_sample_id = sprintf ("${i}_zero_sample_%03d", ++$zero_sample_count{$i});
	print "// Create new sample \"$new_sample_id\" for collection $i\n\n" if ( $verbose );

	# Collection IDs
	$ISA{$new_sample_id}{collection_ID}           = $i;
	$ISA{$new_sample_id}{sample_ID}           = $new_sample_id;
	$ISA{$new_sample_id}{sample_description}  = "Record of absence of some species of mosquito";
	# Collection date
	$ISA{$new_sample_id}{collection_start_date}  = $collection_meta{$i}{collection_start_date};  # collection start date
	$ISA{$new_sample_id}{collection_end_date}    = $collection_meta{$i}{collection_end_date};    # collection end date
	# Collection site location
	$ISA{$new_sample_id}{GPS_latitude}           = $collection_meta{$i}{GPS_latitude};   # collection site latitude
	$ISA{$new_sample_id}{GPS_longitude}          = $collection_meta{$i}{GPS_longitude};   # collection site longitude
	$ISA{$new_sample_id}{collection_description} = $collection_meta{$i}{collection_description};  # collection event description
	# Collection trap metadata
	$ISA{$new_sample_id}{trap_type}              = $collection_meta{$i}{trap_type};        # trap type
	$ISA{$new_sample_id}{attractant}             = $collection_meta{$i}{attractant};  # trap attractant
	$ISA{$new_sample_id}{trap_duration}          = $collection_meta{$i}{trap_duration};    # Duration of trap deployment
	$ISA{$new_sample_id}{trap_number}            = $collection_meta{$i}{trap_number};    # No. of traps deployed
	# Collected material metadate
	$ISA{$new_sample_id}{sample_count}           = 0;   # No. of animals collected
	foreach my $j (@missing) {
	    next if ( $j eq "Culicidae" );	# Ignore generic species, don't make confirmed zero sample size for generic terms
	    next if ( $j eq "Culicinae" );      # I wonder if we could generalise this to all single-word taxonomic terms? <<<<<<
	    next if ( $j =~ /genus/ );

	    push @{ $ISA{$new_sample_id}{Species} },  $j;   #
	}

	# TO SORT OUT STILL:
	$ISA{$new_sample_id}{species_identification_method} = "SPECIES_MORPHO?????";   	# species identification method
	$ISA{$new_sample_id}{sex}                           = "female";   		# sex
	$ISA{$new_sample_id}{developmental_stage}           = "adult";    		# developmental stage
    }
}

##
## s_sample sheet output
##

if ( $s_sample ) {

    my $s_tab = []; # output data structure reference to array of arrays

    open(my $s_fh, ">$output_directory/s_samples.$output_suffix") || die;

    push @{$s_tab}, [ 'Source Name', 'Sample Name', 'Description',
    'Material Type', 'Term Source Ref', 'Term Accession Number',
    'Characteristics [sex (EFO:0000695)]', 'Term Source Ref', 'Term Accession Number',
    'Characteristics [developmental stage (EFO:0000399)]', 'Term Source Ref', 'Term Accession Number',
    'Characteristics [sample size (VBcv:0000983)]' ];

    foreach my $row (values %ISA) {
	push @{$s_tab}, [
	    $config->{study_identifier},
	    $row->{sample_ID},
	    $row->{sample_description} // '',
	    ontology_triplet_lookup("pool", $config->{study_terms}, "strict"),
	    ontology_triplet_lookup($row->{sex}, $config->{study_terms}, "strict"),
	    ontology_triplet_lookup($row->{developmental_stage}, $config->{study_terms}, "strict"),
	    $row->{sample_count}
	];
    }
    print_table($s_fh, $s_tab);
    close($s_fh);
}


##
## a_collection sheet output
##

# [DL] Currently has no provision for adding collection location, GAZ term etc.

if ( $a_collection ) {

    my $c_tab = []; # output data structure reference to array of arrays
    open(my $c_fh, ">$output_directory/a_collection.$output_suffix") || die;

    push @{$c_tab}, [ 'Sample Name', 'Assay Name', 'Description',
		      'Protocol REF', 'Performer', 'Date',
		      'Characteristics [Collection site (VBcv:0000831)]', 'Term Source Ref', 'Term Accession Number',
		      'Characteristics [Collection site latitude (VBcv:0000817)]',
		      'Characteristics [Collection site longitude (VBcv:0000816)]'
    ];
    foreach my $row (values %ISA) {
	push @{$c_tab}, [ $row->{sample_ID}, $row->{collection_ID}, $row->{collection_description} // '',
			  $row->{trap_type},
			  '', # blank Performer
			  $row->{collection_start_date},
			  # pending GAZ retirement...
			  # temporary handling of "Collection site" column using the first available of
			  # location_ADM2>location_ADM1>location_country>location_description
                          # and looking up GAZ terms in the study_terms lookup if available
			  ontology_triplet_lookup($row->{location_ADM2} || 
						  $row->{location_ADM1} ||
						  $row->{location_country} ||
						  $row->{location_description}, $config->{study_terms}, 'relaxed'),
			  $row->{GPS_latitude}, $row->{GPS_longitude},
	];
	# TO DO: collection site coordinates qualifier code or ontology term
	
    }
    print_table($c_fh, $c_tab);
    close($c_fh);
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
	printf OUTPUT ("$ISA{$i}{sample_ID},$ISA{$i}{sample_ID}.spp,,$species_proc,,$ISA{$i}{collection_start_date},", $i );
	my ($sp_species,$sp_onto,$sp_acc);
	foreach my $j ( @{ $ISA{$i}{Species} } ) {
	    my ($onto,$acc) = $config->{study_species}{$j} =~ ( /^(\S+?)\:(\S+)/ );
	    if ( $acc eq "" ) { print "// WARNING: No ontology term for $j $i [$ISA{$i}{collection_ID} : $ISA{$i}{sample_ID}]\n"; }
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
	    print "// $ISA{$i}{sample_ID} [$assay_type,$assay_proc,$assay_result,$assay_units] $j\n" if ( $moreverbose );


	    printf OUTPUT_A ("$ISA{$i}{sample_ID},$ISA{$i}{sample_ID}.${assay_type},,${assay_type},p_virus.txt\n", $i );

	    if ( $assay_result eq "Positive") {
		printf OUTPUT_P ("$ISA{$i}{sample_ID}.${assay_type},\"${assay_type} infected\",\"arthropod infection status\",VSMO,0000009,${assay_type},VSMO,0000535,present,PATO,0000467\n");

	    }
	    elsif ( $assay_result eq "Negative") {
		printf OUTPUT_P ("$ISA{$i}{sample_ID}.${assay_type},\"${assay_type} infection not detected\",\"arthropod infection status\",VSMO,0000009,${assay_type},VSMO,0000882,absent,PATO,0000462\n");
	    }
	}
    }
    close OUTPUT_A;
    close OUTPUT_P;
}


exit(0);

#-----------------------------------------------------------------------------------------------#


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
    # use the object oriented interface to get the sample_IDs in order
    my $hashify = Text::CSV::Hashify->new( { file => $file,
					     key => 'sample_ID',
					     sep_char => $file =~ /\.csv$/ ? ',' : "\t",
					   } );

    my $input_ref = $hashify->all;
    my $sample_IDs = $hashify->keys;

    die "no data read from $file - is sample_ID column present?" unless (defined $input_ref && keys %{$input_ref});

    # copy over the data into global %ISA - using ordered sample_IDs
    foreach my $sample_ID (@$sample_IDs) {
	$ISA{$sample_ID} = $input_ref->{$sample_ID};
    }

    foreach my $key (keys %ISA) {
	my $sample_data = $ISA{$key};

	# TO DO: validate trap_type is in $config hash

	# If trap_type is AWOL set to MIRO:30000045 "catch of live specimens"
	if ($sample_data->{trap_type} eq "") {
	    $sample_data->{trap_type} = "LIVE";  # Note: this will just be the value for the Protocol REF column in the a_collection ISA-sheet
	}

	# handle special case "BLANK" species
	if ($sample_data->{species} eq "BLANK") {
	    print "// Assert all species for the Blank collection $sample_data->{collection_ID} $sample_data->{sample_ID}\n" if ($verbose);
	    $sample_data->{sample_description} = "Record of absence of some species of mosquito";
	    $sample_data->{species} = [];
	    foreach my $j (@expected_species) {
		# Ignore generic species assertion, don't make confirmed zero sample size for generic terms
		next if ( $j eq "Culicidae" );
		next if ( $j eq "Culicinae" );
		next if ( $j =~ /genus/ );
		push @{ $sample_data->{species} },  $j;
	    }
	}

	# Parallel tracking of collection metadata ()
	#
	# Lazy overwriting of (presumed) consistent collection metadata
	for my $column ('collection_start_date', 'collection_end_date', 'GPS_latitude', 'GPS_longitude', 'location_description', 'trap_type', 'attractant', 'trap_duration', 'trap_number') {
	    $collection_meta{$key}{$column} = $sample_data->{$column};
	}
    }

    # Assays metadata and result
    # push @{ $ISA{$f[1]}{Assay} },      "$f[17];$f[18];$f[19];$f[20];";   # assay type

    #
    # # Calculating the maximum sample ID ordinal
    # my ( $ord ) = $f[1] =~ ( /sample_(\d+)/ );
    # if ( $ord > $max_sample_id ) {
    # 	$max_sample_id = $ord;
    # }

    # Tracking for debug/reporting
    # $row_parsed++;
    # 	}
    # }

    # Report number of rows parsed from the input file
    # if ( $verbose ) {
    # 	print "// Parsed file      : $file\n";
    # 	print "// No. rows in file : $row_count\n";
    # 	print "// No. rows parsed  : $row_parsed\n";
    # 	print "// Max. sample ID   : $max_sample_id\n";
    # }
}


##
## calculate diff of elements in 2 arrays
##

sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}


=head2 print_table

args: filehandle, arrayref

prints tab-delimited data to the file.

If you need extra newlines, print them explicitly.

=cut

sub print_table {
    my ($filehandle, $arrayref) = @_;
    foreach my $row (@{$arrayref}) {
	$csv_formatter->print($filehandle, $row);
    }
}


=head2 ontology_triplet_lookup


=cut


sub ontology_triplet_lookup {
    my ($value, $lookup, $mode) = @_;
    my $strict = $mode =~ /strict/i;

    my ($onto, $acc) = ('', '');

    if (defined $value && length($value)) {
	
	my $term_acc = $lookup->{$value};

	if ($term_acc) {
	    ($onto, $acc) = $term_acc =~ (/(\S+?)\:(\S+)/);
	    if (defined $onto && defined $acc) {
		return ($value, $onto, $acc);
	    } elsif ($strict) {
		die "malformed ontology term accession '$term_acc' in config file\n";
	    } elsif ($verbose) {
		warn "malformed ontology term accession '$term_acc' in config file\n";
	    }
	} elsif ($strict) {
	    die "ontology term '$value' not defined in config file\n";
	} elsif ($verbose) {
	    warn "ontology term '$value' not defined in config file\n";
	}
    } elsif ($strict) {
	die "empty value passed to ontology_triplet_lookup()\n";
    } elsif ($verbose) {
	warn "empty value passed to ontology_triplet_lookup()\n";
	$value = '';
    } else {
	$value = '';
    }
    return ($value, $onto, $acc);
}


#-----------------------------------------------------------------------------------------------#
