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
my $rows;			# Count of rows in file
my $errors;  		# Number of errors found in file
my $warnings;		# Count of warnings
my $data_rows;		# Count of data rows
my $comments;		# Count of comment rows
my $blanks;		# Count of blank rows


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

# Assign error output filename
my $output = $file . "_output.txt";

##
## Parameters for validation checks
##

my $expected_columns = 21;		# Number of columns expected in the interchange format.

open (OUTPUT, "> $output");


# We are parsing csv files
my $csv = Text::CSV->new({ binary => 1 });              # create a new object

open (my $data , "< $file");
while (my $line = <$data>) {
	chomp $line;
	$rows++;

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
		## Min/Max values for latitude and longitude
		##
		if ( ($f[5] < -90) or ($f[5] > 90) ) {
			print OUTPUT "// WARNING: Latitude coordinate is out of range (-90 < Lat < 90), $f[5] at $rows \n";
			$warnings++;
			next;
		}
		if ( ($f[6] < -180) or ($f[6] > 180) ) {
			print OUTPUT "// WARNING: Longitude coordinate is out of range (-180 < Lat < 180), $f[6] at $rows\n";
			$warnings++;
			next;
		}

		##
		## Sample size must be an integer (warn if zero?)
		##
		if ( $f[16] eq "" ) {
			print OUTPUT "// ERROR: Collected sample size is absent $f[20] in line $rows \n";
			$errors++;
		}
		elsif ( $f[16] == 0 ) {
			print OUTPUT "// WARNING: Collected sample size is zero $f[20] in line $rows \n";
			$warnings++;
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
