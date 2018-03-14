#!/usr/bin/perl

##
## InterchangeValidator.pl
##
## 2018 Dan Lawson (daniel.lawson@imperial.ac.uk)
## v1.0 [DL] pragmatic perl

use Text::CSV;
use Getopt::Long;
use IO::Handle;
use strict;

my $verbose;
my $file;
my $help;
my $summary;
my %species;		# Hash of species present in interchange file

my %collection_date;	# Hash of collection dates
my %collection_desc;	# Hash of collection descriptions
my %collection_protocol;	# Hash of collection protocols (traap type)
my %VBsp;			# Hash of valid VBsp names and accessionss
my $rows;			# Count of rows in file
my $errors;  		# Number of errors found in file
my $warnings;		# Count of warnings
my $data_rows;		# Count of data rows
my $comments;		# Count of comment rows
my $blanks;		# Count of blank rows

##
## Parameters for validation checks
##

my $expected_columns = 21;		# Number of columns expected in the interchange format.


#---------------------------------------------------------#
GetOptions (
	# Misc
	"verbose"		=> \$verbose,
	"help"		=> \$help,
	# Data
	"file=s"		=> \$file,
	# Mode
	"summary" 	=> \$summary,
);
#---------------------------------------------------------#

# pod documentation for help
if ( $help ) {
  exec ('perldoc',$0);
}

# catch simple errors
unless ( $file ) {
	print "No input file supplied\n";
	exit(0);
}


&parse_spreadsheet;


# Assign error output filename
my $output = $file . "_output.txt";


open (OUTPUT, "> $output");


# We are parsing csv files
my $csv = Text::CSV->new({ binary => 1 });              # create a new object

open (my $data , "< $file");
while (my $line = <$data>) {
	chomp $line;
	$rows++;

	# Ignore column headings (first line of a valid file)
	next if ( $rows == 1);

	if ($csv->parse($line)) {
		my @f = $csv->fields();

		##
		## Ignore column constraint for comment lines
		##
		if ( $line =~ /^\/\// ) {
			$comments++;
			next;
		}

		##
		## Ignore column constraint for blank lines
		##
		if ( $line eq "" ) {
			$blanks++;
			next;
		}

		##
		## Check for correct number of columns (21)
		##
		if ( scalar @f != $expected_columns ) {
			print OUTPUT "// ERROR: Incorrect number of columns, " . (scalar @f) . " in line $rows\n";
			$errors++;
		}

		##
		## Missing GPS coordinates
		##
		if ( ( $f[5] eq "" ) or ( $f[6] eq "" ) ) {
			print OUTPUT "// ERROR: GPS lat-lon coordinates are missing for sample $f[1] at $rows\n";
			$errors++;
		}

		##
		## Min/Max values for latitude and longitude
		##
		if ( ($f[5] < -90) or ($f[5] > 90) ) {
			print OUTPUT "// WARNING: Latitude coordinate is out of range (-90 < Lat < 90), $f[5] at $rows \n";
			$warnings++;
		}
		if ( ($f[6] < -180) or ($f[6] > 180) ) {
			print OUTPUT "// WARNING: Longitude coordinate is out of range (-180 < Lat < 180), $f[6] at $rows\n";
			$warnings++;
		}

		##
		## Consistency of collection information
		##
		if ( ! defined $collection_date{$f[0]} ) {
			$collection_date{$f[0]} = $f[2] . "-" . $f[3];
		}
		elsif ( $collection_date{$f[0]} ne "$f[2]-$f[3]" ) {
			print OUTPUT "// ERROR: Collection $f[0] has multiple dates [$collection_date{$f[0]} v $f[2]-$f[3]]\n";
			$errors++;
		}

		if ( ! defined $collection_desc{$f[0]} ) {
			$collection_desc{$f[0]} = $f[7];
		}
		elsif ( $collection_desc{$f[0]} ne $f[7] ) {
			print OUTPUT "// ERROR: Collection $f[0] has multiple descriptions [\"$collection_desc{$f[0]}\" v \"$f[7]]\"\n";
			$errors++;
		}

		if ( ! defined $collection_protocol{$f[0]} ) {
			$collection_protocol{$f[0]} = $f[8];
		}
		elsif ( $collection_protocol{$f[0]} ne $f[8] ) {
			print OUTPUT "// ERROR: Collection $f[0] has multiple protocols (trap type) [\"$collection_protocol{$f[0]}\" v \"$f[8]\"]\n";
			$errors++;
		}


		##
		## Check validity of the species assingment
		##
		if ( (! defined $VBsp{$f[12]}) && ( $f[12] ne "BLANK") ) {
			print OUTPUT "// WARNING: Asserted species $f[12] does not have a VBsp accession\n";
			$warnings++;
		}


		##
		## Sample size must be an integer (warn if zero?)
		##
		if ( $f[16] eq "" ) {
			print OUTPUT "// ERROR: Collected sample size is absent ($f[1] => Species $f[12]) in line $rows \n";
			$errors++;
		}
		elsif ( $f[16] == 0 ) {
			print OUTPUT "// WARNING: Collected sample size is zero ($f[1] => Species: $f[12]) in line $rows \n";
			$warnings++;
		}

		##
		## Summary information
		##
		if ( $summary ) {
			$species{$f[12]}++;
		}


		# Increment valid data row count
		$data_rows++;
	}

}
close FILE;
close OUTPUT;

##
## Write summary information to terminal
##

print "// File: $file\n";
print "//\n";
printf ("// No. lines        : %5d\n", $rows);
printf ("// No. comment rows : %5d\n", $comments);
printf ("// No. blank rows   : %5d\n", $blanks);
printf ("// No. data rows    : %5d\n", $data_rows);
print "//\n";
printf ("// No. warnings     : %5d\n", $warnings);
printf ("// No. errors       : %5d\n", $errors);
print "//\n";

##
## Write errors and warnings to terminal if required
##

if ( $verbose ) {
	open (ERROR, "< $output");
	while (<ERROR>) {
		print;
	}
}


##
## Summary of all species
##

if ( $summary ) {

	print "// No. species present in $file : " . (keys %species) . "\n";

	foreach my $i (sort keys %species) {
		next if ($i eq "BLANK");
		printf ("%-32s $VBsp{$i} %4d\n", $i,$species{$i});
	}

	if ( $species{BLANK} ) {
		print "\nBLANK (no mosquitoes collected)   $species{BLANK}\n";
	}
}

sub parse_spreadsheet {

	my $key = "1Km3YHrBMtREyjPoR0TdOxCdxHjBEZzq0Sgfka_Db6Tw";
	my $curlcommand = "https://docs.google.com/spreadsheets/d/".$key."/export?exportFormat=tsv";
	my @f;

	# Get species metadata via curl and process
	#---------------------------------------------------------#

	open (FILE, "curl --silent $curlcommand |") or print "// WARNING failed to open spreadsheet\n";
	while (<FILE>) {
		chomp;
		@f = split/\t/;
		$VBsp{$f[1]} = $f[0];
#		print "$f[1] => $f[0]\n";
	}
	close FILE;
	return (%VBsp);
}


=pod

=head1 NAME - InterchangeValidator.pl

=head2 USAGE

This scripts checks a VectorBase PopBio interchange format csv file for syntax and data consistancy. Summary information is returned to the terminal

=head2 ARGUMENTS

B<InterchangeValidator.pl> arguments:


=head3 Input file

=over 4

=item -file, Select a comma-separated value file (csv) as input for validation

=back

=head3 Misc

=over 4

=item -verbose, Verbose mode toggle on extra command line output

=item -help, these help pages

=back

=head1 AUTHOR

Dan Lawson (daniel.lawson@imperial.ac.uk)

=cut
