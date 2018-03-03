#!/usr/bin/perl
#########################################################################################################
##  convert to interchange format                                                                      ##
##  Manatee County Mosquito Control District (MCMCD)                                                   ##
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
my $help;
my $prefix;
my $summary;
my $writesamples     	= 0;
my $writecollections 	= 0;
my %seen_species;
my %species_seen;
my %collection_history;
my %location_lat;
my %location_lon;
my %location_description;
my %location_count;
my @species;

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
my $outfile           = "${prefix}.pop";


my @columns;		# Column definitions

print "// Project: $project\n" if ($verbose);

# read species names and VBsp terms (why, is this used?)
my %species = &species_data;

# read location GPS coordinates and description from file
&location_data;


open (OUTPUT, "> $outfile");
print OUTPUT "collection,sample,collection_date_start,collection_date_end,trap_id,GPS_lat,GPS_lon,location_description,trap_type,attractant,trap_number,trap_duration,species,identification_method,sex,sample_count,assay_type,assay_method,assay_result,assay_units\n";

my $csv = Text::CSV->new({ binary => 1 } );              # create a new object

my $collection = 0;
my $sample     = 0;

## "",
## "Date",
## "Type",
## "Area",
## "Ref area/comp#",
## "Collection site latitude",
## "Collection site longitude",
## "Collection site locality",
## "variable",
## "Sample size",
## "Protocol REF",
## "Collection site county",
## "Collection site province",
## "Characteristics"

#Date,Type,Area,Ref area/comp#,Culex nigripalpus,Aedes taeniorhynchus,Psorophora ciliata,Psorophora columbiae,Psorophora ferox,Mansonia titillans,Coquillettidia perturbans,Anopheles crucians,Anopheles quadrimaculatus,Aedes vexans,Aedes aegypti,Aedes albopictus,Aedes infirmatus,Culex erraticus,Culex iolambdis,Culex pilosus,Culex (mel) spissipes,Culex quinquefasciatus,Culex salinarius,Anopheles albimanus,Anopheles atroparvus,Anopheles barberi,Anopheles bradleyi,Anopheles georgianus,Anopheles punctipennis,Anopheles perplexans,Anopheles walkeri,Aedes atlanticus,Aedes atlanticus-tormentor,Aedes bahamensis,Aedes canadensis,Aedes canadensis mathesoni,Aedes cinereus,Aedes dupreei,Aedes fulvus pallens,Aedes hendersoni,Aedes mitcheilae,Aedes sollicitans,Aedes sticticus,Aedes thelcter,Aedes thibaulti,Aedes tormentor,Aedes tortilis,Aedes triseriatus,Aedes spp.,Culex coronator,Culex restuans,Culex tarsalis,Culex territians,Culex spp.,Culex atratus,Culex mulrennani,Culex opisthopus,Culex peccator,Culiseta inornata,Culiseta melanura,Psorophora cyanescens,Psorophora discolor,Psorophora horrida,Psorophora howardii,Psorophora johnstoni,Psorophora pygmaea,Psorophora mathesoni,Mansonia dyari,Wyeomyia haynei,Wyeomyia mitchellii,Wyeomyia vanduzeei,Deinocerites cancer,Toxorhychites rutilus rutilus,Toxorhychites rutilus spetentrionalis,Orthopodomyia signifera,Orthopodomyia signifera alba,Uranotaenia lowii,Uranotaenia sapphirina,Anopheles spp.,Psorophora spp.,Mansonia spp.,Wyeomyia spp.,Toxyrhychites spp.,Orthopodomyia  spp.,Uranotaenia spp.,Culicines,Unidentified,Total,Ae.Total,Ps.Total,Mn.Total,Cx. Total,An.Total,Other total

my $total_count;
my $line_count;

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


		# Track collections by portmanteau of time/location/trap_type
		my $id = "$f[0]_$f[2]_$f[3]";
		print "ID :: $id\n" if ( $moreverbose );
		print  "// WARNING - seen this ID before, assign to existing collection \n" if ( ($collection_history{$id} ne "" ) && ($verbose) );

		if ( $summary ) {
			print "// Pushing \"$columns[4]\" to ID: $id\n" if ($verbose);
			push @{ $species_seen{$id} }, $f[8];
			print "// Array $id has size " . scalar @{ $species_seen{$id} } . "\n" if ($verbose);
			print "$id :: $species_seen{$id}[0]\n" if ($verbose);
			my $content = join " - ", @{ $species_seen{$id} };
			print "$id :: $content\n" if ($verbose);
		}
		else {

			## Species column counts
			## $f[4] to $f[28]

			for ( my $i = 4; $i <= 87; $i++ ) {

				# deal with total column
				next if ( ( $i == 87 ) && ( $f[$i] != 0 ) );

				# discard zero value cells
				if ( $f[$i] eq "" ) {
					print "// Confirmed absence, zero sample size, for $columns[$i]\n" if ($moreverbose);
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
					print OUTPUT "$f[0],$f[0],";

					# Trap ID
					print OUTPUT "$id,";

					# [WHERE] Collection location
					print OUTPUT "$location_lat{$f[2]},$location_lon{$f[2]},\"$location_description{$f[2]}\",";

					# [HOW] Trap type
					print OUTPUT "CDC_LIGHT,\"LIGHT,CO2\",1,1,";


					# [WHAT] Species and sample_count
					## Deal with edge case of no mosquitoes found in the collection
					if ( $i == 87 ) {
						if ( $f[$i] == 0 ) {
							print OUTPUT "BLANK,SPECIES_MORPHO,adult,female,0,,,,\n";
						}
						else {
							print OUTPUT "// LOOK It's the total column ($f[$i])\n";
						}
					}
					else {
						print OUTPUT "$columns[$i],SPECIES_MORPHO,adult,female,$f[$i],,,,\n";
						# add sample_count to total_count
						$total_count = $total_count + $f[$i];
						$seen_species{$columns[$i]} = $seen_species{$columns[$i]} + $f[$i];
					}
				}
			}
		}
	}
	else {
	 	print "// parsing line failed [$line]\n" if ($verbose);
	}

	# increment collection count
	$collection++;
}
close FILE;

close OUTPUT;

if ($verbose) {
	print "// Parsed $line_count rows and created " . $sample  . " samples from " . $collection . " collections\n";
	print "// Total number of mosquitoes: $total_count from " . keys(%seen_species) ." species\n";
}

if ( $summary ) {

	foreach my $i (sort keys %species_seen) {


		# print array contents (redundant list of species seen)

		printf ("%-36s :: [", $i);
		print scalar @{ $species_seen{$i} };
		print "] :: ";
		print join ',', @{ $species_seen{$i} };
		print "\n";

		# print distinct array contents (non-redundant list of species_seen)
		my @distinct = uniq( @{ $species_seen{$i} } );
		printf ("%-36s :: [", $i);
		print scalar @distinct;
		print "] :: ";
		print join ',', @distinct;
		print "\n";

		# print list of AWOL species based on the full enumerated list in %species

		printf ("%-36s :: [", $i);
		my @missing = grep { ! ({ $_, 0 } ~~ @distinct) } @species;
		print scalar @missing;
		print "] :: ";
		print join ", ", @missing;
		print "\n";

		print "\n";
	}
}



sub species_data {

my %species = (
 	'Aedes aegypti'				=> '0000518',
	'Aedes albopictus'				=> '0000522',
	'Aedes atlanticus'				=> '0000976',
	'Aedes atlanticus-tormentor'			=> '0003892',
	'Aedes bahamensis'				=> '0000895',
	'Aedes canadensis'				=> '0003846',
	'Aedes canadensis mathesoni'			=> '0003871',
	'Aedes cinereus'				=> '0000255',
	'Aedes dupreei'				=> '0001024',
	'Aedes fulvus pallens'			=> '0001039',
	'Aedes hendersoni'				=> '0001188',
	'Aedes infirmatus'				=> '0001059',
	'Aedes mitchellae'				=> '0001081',
	'Aedes sollicitans'				=> '0001138',
	'Aedes sp.'				=> '0000253',
	'Aedes sticticus'				=> '0001144',
	'Aedes taeniorhynchus'			=> '0001152',
	'Aedes thelcter'				=> '0001154',
	'Aedes thibaulti'				=> '0001156',
	'Aedes tormentor'				=> '0001157',
	'Aedes tortilis'				=> '0001158',
	'Aedes triseriatus'				=> '0001206',
	'Aedes vexans'				=> '0000372',
	#Anopheles
	'Anopheles albimanus'			=> '0000071',
	'Anopheles atropos'				=> '0000033',
	'Anopheles barberi'				=> '0000040',
	'Anopheles bradleyi'			=> '0000047',
	'Anopheles crucians'			=> '0000072',
	'Anopheles georgianus'			=> '0000094',
	'Anopheles perplexens'			=> '0003416',
	'Anopheles punctipennis'			=> '0003439',
	'Anopheles quadrimaculatus'			=> '0003441',
	'Anopheles spp.'				=> '0000015',
	'Anopheles walkeri'				=> '0003469',
	# Coquillettidia
	'Coquillettidia perturbans'			=> '0002347',
	# Culex
	'Culex atratus'				=> '0003008',
	'Culex coronator'				=> '0002544',
	'Culex erraticus'				=> '0003050',
	'Culex iolambdis'				=> '0003071',
	'Culex mulrennani'				=> '0003090',
	'Culex nigripalpus'				=> '0002622',
	'Culex opisthopus'				=> '0003891',
	'Culex peccator'				=> '0003103',
	'Culex pilosus'				=> '0003110'.
	'Culex quinquefasciatus'			=> '0002654',
	'Culex restuans'				=> '0002657',
	'Culex salinarius'				=> '0002661',
	'Culex spissipes'				=> '0003131',
	'Culex spp.'				=> '0002423',
	'Culex tarsalis'				=> '0002687',
	'Culex territans'				=> '0003218',
	# Culiseta
	'Culiseta inornata'				=> '0002409',
	'Culiseta melanura'				=> '0002381',
	'Culiseta spp.'				=> '0002373',
	# Deinocerites
	'Deinocerites cancer'			=> '0002287',
	# Mansonia
	'Mansonia dyari'				=> '0001257',
	'Mansonia spp.'				=> '0001252',
	'Mansonia titillans'			=> '0001265',
	# Orthopodomyia
	'Orthopodomyia signifera'			=> '0003681',
	'Orthopodomyia signifera alba'		=> '0001302',
	'Orthopodomyia spp.'			=> '0001301',
	# Psorophora
	'Psorophora ciliata'			=> '0001345',
	'Psorophora columbiae'			=> '0001307',
	'Psorophora cyanescens'			=> '0001326',
	'Psorophora discolor'			=> '0001310',
	'Psorophora ferox'				=> '0001328',
	'Psorophora horrida'			=> '0001331',
	'Psorophora howardii'			=> '0001348',
	'Psorophora johnstonii'			=> '0001332',
	'Psorophora mathesoni'			=> '0001336',
	'Psorophora pygmaea'			=> '0001317',
	'Psorophora spp.'				=> '0001304',
	# Toxorhynchites
	'Toxorhynchites rutilus rutilus'		=> '0001484',
	'Toxorhynchites rutilus septentrionalis'	=> '0001485',
	'Toxorhynchites spp.'			=> '0001454',
	# Uranotaenia
	'Uranotaenia lowii'				=> '0002089',
	'Uranotaenia sapphirina'			=> '0002128',
	'Uranotaenia spp.'				=> '0001927',
	# Wyeomyia
	'Wyeomyia haynei'				=> '0001857',
	'Wyeomyia mitchellii'			=> '0001869',
	'Wyeomyia spp.'				=> '0001772',
	'Wyeomyia vanduzeei'			=> '0001884',
	# Misc
	'Culicinae'                              	=> '0003653',
	'Culicidae'                             	=> '0003818',
);

print "// read species, " . keys(%species) . " species\n" if ( $verbose );

foreach my $j (sort keys %species) {
	 push @species, $j;
 }

return %species;
}

sub location_data {
	my $location_count;
	while (<DATA>) {
		chomp;
		my @f = split /\t/;
		$location_lat{$f[0]} = $f[1];
		$location_lon{$f[0]} = $f[2];
		$location_description{$f[0]} = $f[3];
		print "$f[0]\t$location_description{$f[0]}\t$location_lat{$f[0]}\t$location_lon{$f[0]}\n" if ($verbose);
		$location_count++;
	}
	print "// read location GPS coordinates from file: , $location_count locations\n" if ($verbose);
	return (%location_lon,%location_lat,%location_description);
}


sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}


__DATA__
B1	27.616367	-82.535300	Chapman Rd.
B1CH	27.600978	-82.529357	Church
B1F	27.603067	-82.535717	French
B1FC	27.590917	-82.548000	Frog Creek
B1GF	27.627967	-82.496833	Grassfarm Rd.
B1P	27.640483	-82.552667	Port
B2-1	27.575706	-82.577527	Ken Hubbard Rd.
B2-2	27.574542	-82.580910	Center Rd.
B2-3	27.581486	-82.587467	Terra Ceia Rd.
B2-4	27.564890	-82.587227	Horseshoe Loop Rd.
B2S	27.569469	-82.593177	Shiver
B3EP	27.532041	-82.628181	Emerson Point
B4	27.529167	-82.588117	Fielding
B4BP	27.531000	-82.582350	Palmetto Ball Park
B4HS	27.516317	-82.537735	Highland Shores
B4J	27.554018	-82.545481	Carrs
B4S	27.533281	-82.565712	Palmetto shop
C1	27.503120	-82.657390	Orbans
C2P	27.511090	-82.683495	Perico
C3C	27.468452	-82.684313	Cortez Village
C4	27.447151	-82.620520	Bolleteries Soccer field
D1	27.488870	-82.592140	Manatee HS
D1BP	27.475697	-82.613417	Colonial Baptist Church
D1ST	27.499288	-82.600544	St Stephens
D2BP1	27.489937	-82.533814	Lloyd Park
D2BP2	27.448930	-82.498664	Braden River Park
D2F	27.475085	-82.551283	Fulwood
D2S	27.469470	-82.524020	Samoset
D3	27.439449	-82.548928	Mullins
D3P	27.406663	-82.504054	Palm-Aire
D3PM	27.419700	-82.479820	Mote Ranch
D3T	27.426600	-82.524710	Tallevast
F1M	27.481256	-82.480420	Magnolia Manor
F2L	27.419517	-82.426200	Lakewood Ranch
F2SF	27.426485	-82.376623	Soccer Field
F3	27.484586	-82.401610	Rye Rd
F3B	27.438890	-82.438650	Braden Pines
F3BP	27.442420	-82.435780	Lakewood Ranch HS
F3SF	27.435822	-82.391094	Soccer Field
F3U	27.512610	-82.431900	Upper Manatee
F3WL	27.489400	-82.361667	Waterline Rd
F4	27.463570	-82.303790	Jamie's Way
F4A	27.412180	-82.287570	Apple
F4R	27.474551	-82.324001	Race Track
F4S	27.407790	-82.304000	South 675
F4SM	27.496665	-82.308603	Swift Mud
F4V	27.408050	-82.260329	Gopher Hill
F5-1	27.347367	-82.154433	Myakka
F5-2	27.346702	-82.158998	Myakka
F5-3	27.346935	-82.160157	Myakka
F5M	27.349733	-82.157581	Myakka
F5MH	27.286571	-82.230380	Mossy Hammock
F5S	27.347120	-82.243525	Singletary Rd
F5SB	27.272953	-82.116503	Sandy Baptist
F6	27.426958	-82.133258	Wauchula Rd.
F6B	27.470780	-82.169930	Bear Bay Rd.
G1BP	27.578355	-82.483695	Buffalo Creek Ball Park
G1C	27.639167	-82.417083	Prichard Rd.
G1E	27.574079	-82.504100	Erie Rd.
G1T	27.564989	-82.480332	Thousand Oaks
G2FBb	27.577557	-82.354793	Fox Brook back
G2FBf	27.549098	-82.353430	Fox Brook front
G2M	27.575050	-82.419667	Parrish
G2MS	27.574840	-82.419638	Parrish
G2P	27.564967	-82.424886	Parrish/Britt Rd.
G3B	27.591467	-82.199850	Bunker Hill Church
G3BR	27.603150	-82.098133	Bradley Rd.
G3CR	27.613375	-82.123478	Carlton Rd.
G3K	27.570935	-82.137841	Keentown
