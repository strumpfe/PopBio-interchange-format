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
GetOptions(
    # Misc
    "verbose"       => \$verbose,
    "verboser"      => \$moreverbose,
    "help"          => \$help,
    # Data
    "file=s"        => \$file,
    "config=s"      => \$configfile,
    # Process
    "zeros|zeroes"         => \$add_zeros,
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
print "\n// Sexes expected:\n";
foreach my $i (keys %{$config->{study_sexes}}) {
    printf ("//  %-30s $config->{study_sexes}{$i}\n", $i);
}
print "\n// Developmental stages expected:\n";
foreach my $i (keys %{$config->{study_developmental_stages}}) {
    printf ("//  %-30s $config->{study_developmental_stages}{$i}\n", $i);
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
    # collection_ID => species => sex => developmental_stage => count
    my %collection_species_sex_stage;
    tie %collection_species_sex_stage, 'Tie::Hash::Indexed'; # preserve collection order

    # collection_ID => a 'row' record in %ISA 
    my %collection_row;

    # counter for the ordinal added to zero sample IDs: <collection_ID>_zero_sample_NNN
    my %zero_sample_count;

    # Loop through hash to collate seen species-sex-stage combos
    # also look out for species eq 'BLANK' rows
    foreach my $sample_ID ( keys %ISA ) {
        my $row = $ISA{$sample_ID};
        my $collection_ID = $row->{collection_ID};
        my $species = $row->{species};

        if ($species eq 'BLANK') {
            # make an empty hash for this $collection_ID so we know it needs
            # zero samples creating for it
            $collection_species_sex_stage{$collection_ID} = {};
        } else {
            $collection_species_sex_stage{$collection_ID}{$species}{$row->{sex}}{$row->{developmental_stage}}++;
        }

        $collection_row{$collection_ID} //= $row;
    }

    # Loop through collections that may need augmenting
    foreach my $collection_ID ( keys %collection_species_sex_stage ) {
        # loop through all species-sex-stage combinations and figure out which species are missing for each sex-stage combo
        my %sex_stage_missing_species; # sex => stage => species => 1
        foreach my $species (@expected_species) {
            foreach my $sex (keys %{$config->{study_sexes}}) {
                foreach my $stage (keys %{$config->{study_developmental_stages}}) {
                    unless (exists $collection_species_sex_stage{$collection_ID}{$species}{$sex}{$stage}) {
                        $sex_stage_missing_species{$sex}{$stage}{$species} = 1;
                    }
                }
            }
        }

        # create new samples in %ISA for the missing species for each sex-stage combo
        foreach my $sex (keys %sex_stage_missing_species) {
            foreach my $stage (keys %{$sex_stage_missing_species{$sex}}) {
                my @missing_species = keys %{$sex_stage_missing_species{$sex}{$stage}};

                # remove higher taxonomic levels - typically single word names (so look for whitespace)
                # and also those containing 'genus'
                @missing_species = grep !/\bgenus\b/, grep /\s/, @missing_species;

                if (@missing_species) {
                    my $n_species = @missing_species;
                    my $new_sample_id = sprintf ("${collection_ID}_zero_sample_%03d", ++$zero_sample_count{$collection_ID});
                    print "// Created new sample '$new_sample_id' for collection $collection_ID $sex $stage $n_species species\n" if ( $verbose );

                    # do a full copy of the representative row for this collection
                    my $new_row = $ISA{$new_sample_id} = { %{$collection_row{$collection_ID}} };

                    # set the sample ID
                    $new_row->{sample_ID} = $new_sample_id;
                    # set the stage and sex
                    $new_row->{sex} = $sex;
                    $new_row->{developmental_stage} = $stage;
                    # and the species as an arrayref
                    $new_row->{species} = \@missing_species;
                    # and of course the sample size
                    $new_row->{sample_count} = 0;
                    # and a description
                    $new_row->{sample_description} = "Record of absence of some species ($sex, $stage)";
                }
            }
        }
    }
}

##
## s_sample sheet output
##

if ( $s_sample ) {
    print "// Writing $output_directory/s_samples.$output_suffix sheet\n" if ($verbose);

    my $s_tab = []; # output data structure reference to array of arrays

    open(my $s_fh, ">$output_directory/s_samples.$output_suffix") || die;

    push @{$s_tab}, [
        'Source Name', 'Sample Name', 'Description',
        'Material Type', 'Term Source Ref', 'Term Accession Number',
        'Characteristics [sex (EFO:0000695)]', 'Term Source Ref', 'Term Accession Number',
        'Characteristics [developmental stage (EFO:0000399)]', 'Term Source Ref', 'Term Accession Number',
        'Characteristics [sample size (VBcv:0000983)]',
        'Comment [sample_comment]'
    ];

    foreach my $row (values %ISA) {
        next if ($row->{species} eq 'BLANK');
        push @{$s_tab}, [
            $config->{study_identifier},
            $row->{sample_ID},
            $row->{sample_description} // '',
            ontology_triplet_lookup("pool", $config->{study_terms}, "strict"),
            ontology_triplet_lookup($row->{sex}, $config->{study_sexes}, "strict"),
            ontology_triplet_lookup($row->{developmental_stage}, $config->{study_developmental_stages}, "strict"),
            $row->{sample_count},
            $row->{sample_comment} // '',
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
    print "// Writing $output_directory/a_collection.$output_suffix sheet\n" if ($verbose);

    my $c_tab = []; # output data structure reference to array of arrays
    open(my $c_fh, ">$output_directory/a_collection.$output_suffix") || die;

    push @{$c_tab}, [
        'Sample Name', 'Assay Name', 'Description',
        'Protocol REF', 'Performer', 'Date',
        'Characteristics [Collection duration in days (VBcv:0001009)]',
        'Characteristics [Number of traps (VBcv:0001122)]',
        'Comment [collection_comment]',
        'Comment [Trap ID]',
        'Characteristics [Attractant (IRO:0000034)]', 'Term Source Ref', 'Term Accession Number',
        'Characteristics [Collection site (VBcv:0000831)]', 'Term Source Ref', 'Term Accession Number',
        'Characteristics [Collection site latitude (VBcv:0000817)]',
        'Characteristics [Collection site longitude (VBcv:0000816)]',
        'Comment [collection site coordinates]',
    ];

    foreach my $row (values %ISA) {
        next if ($row->{species} eq 'BLANK');

        ### process the placename info
        # Do the following, pending GAZ retirement...
        # try these in the following order
        # location_description>location_ADM2>location_ADM1>location_country
        # first GAZ term found in study_terms is put into the "Collection site" column
        # if placenames are provided but no GAZ term found, then a warning will be emitted
        my ($location_name, $location_onto, $location_acc) = ('', '', '');
        my $tries = 0;
        foreach my $place ($row->{location_description}, $row->{location_ADM2},
                           $row->{location_ADM1}, $row->{location_country}) {
            next unless defined $place && length($place);

            my ($name, $onto, $acc) = ontology_triplet_lookup($place, $config->{study_terms}, 'relaxed');
            if ($onto && length($acc)) {
                ($location_name, $location_onto, $location_acc) = ($place, $onto, $acc);
                last;
            }
            $tries++;
        }
        warn sprintf "WARNING: couldn't find placename ontology term for descrip: '%s' ADM2: '%s' ADM1: '%s' country: '%s'\n",
                     $row->{location_description} // '', $row->{location_ADM2} // '', $row->{location_ADM1} // '', $row->{location_ADM1} // ''
              if ($tries && !length($location_acc)); 

        push @{$c_tab}, [
            $row->{sample_ID}, $row->{collection_ID}, $row->{collection_description} // '',
            $row->{trap_type},
            '', # blank Performer
            $row->{collection_start_date} eq $row->{collection_end_date} ? $row->{collection_end_date} : "$row->{collection_start_date}/$row->{collection_end_date}",
            $row->{trap_duration} // '',
            $row->{trap_number} // '',
            $row->{collection_comment} // '',
            $row->{trap_ID} // '',
            ontology_triplet_lookup($row->{attractant}, $config->{study_terms}, 'relaxed'),
            $location_name, $location_onto, $location_acc,
            # Lat and long easy
            $row->{GPS_latitude}, $row->{GPS_longitude},
            $row->{GPS_qualifier} // '',
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
    print "// Writing $output_directory/a_species.$output_suffix sheet\n" if ($verbose);

    my $sp_tab = []; # output data structure reference to array of arrays
    open(my $sp_fh, ">$output_directory/a_species.$output_suffix") || die;

    push @{$sp_tab}, [
        'Sample Name', 'Assay Name', 'Description',
        'Protocol REF', 'Performer', 'Date',
        'Characteristics [species assay result (VBcv:0000961)]', 'Term Source Ref', 'Term Accession Number',
        'Comment [species_comment]'
    ];

    foreach my $row (values %ISA) {
        next if ($row->{species} eq 'BLANK');

        if (ref($row->{species}) eq 'ARRAY') {
            # do all the ontology lookups for the multiple species
            my (@ontos, @accs);
            foreach my $species (@{$row->{species}}) {
                my ($temp_spp, $onto, $acc) = ontology_triplet_lookup($species, $config->{study_species}, "strict");
                push @ontos, $onto;
                push @accs, $acc;
            }

            push @{$sp_tab}, [
                $row->{sample_ID}, "$row->{sample_ID}.spp", '',
                $row->{species_identification_method}, '', '',
                join(';', @{$row->{species}}), join(';', @ontos), join(';', @accs),
                $row->{species_comment} // '',
            ];
        } else {
            push @{$sp_tab}, [
                $row->{sample_ID}, "$row->{sample_ID}.spp", '',
                $row->{species_identification_method}, '', '',
                ontology_triplet_lookup($row->{species}, $config->{study_species}, "strict"),
                $row->{species_comment} // '',
            ];
        }
    }
    print_table($sp_fh, $sp_tab);
    close($sp_fh);
}


##
## a_virus sheet output
##

if ( $a_virus ) {
    print "// Writing $output_directory/a_virus.$output_suffix sheet\n" if ($verbose);

    die "virus sheets handling TO BE REIMPLIMENTED...";
    # TO DO:
    # put data for multiple assays (e.g. tests for different viruses) in one row of SAF
    # possibly using vertical bar as the separator

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
            } elsif ( $assay_result eq "Negative") {
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

##
## primary key is the sample ID
##

sub get_data_from_file {
    # use the object oriented interface to get the sample_IDs in order
    my $hashify = Text::CSV::Hashify->new({
        file => $file,
        key => 'sample_ID',
        sep_char => $file =~ /\.csv$/ ? ',' : "\t",
    });

    my $input_ref = $hashify->all;
    my $sample_IDs = $hashify->keys;

    die "no data read from $file - is sample_ID column present?" unless (defined $input_ref && keys %{$input_ref});

    # copy over the data into global %ISA - using ordered sample_IDs
    foreach my $sample_ID (@$sample_IDs) {
        $ISA{$sample_ID} = $input_ref->{$sample_ID};
    }

    # VALIDATION
    foreach my $sample_ID (keys %ISA) {
        my $row = $ISA{$sample_ID};

        if ($row->{trap_number} && $row->{trap_number} > 1) {
            die "trap_number values greater than 1 not yet handled in ISA-Tab/Chado\n";
        }

        # TO DO: validate trap_type is in $config hash
        #        validate GPS coords
        #        check collection fields are unique per collection_ID
        #        check species are in $config->{study_species}
        #        can check sex, attractant, etc fields too
        #        ...?
    }


    # Report number of rows parsed from the input file
    if ( $verbose ) {
        print "// Parsed file      : $file\n";
        printf "// No. rows parsed  : %d\n", scalar keys %ISA;
    }
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
