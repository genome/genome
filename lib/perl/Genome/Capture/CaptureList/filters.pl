#!/usr/bin/env genome-perl
#Compare the first file with the rest of the files in terms of overlap in coordinates
use strict;
use warnings;
use Getopt::Std;
use DBI;
use strict;
use lib "/gscuser/jwallis/svn/perl_modules/test_project/jwallis";

my %opts = (g=>500, f=>0.8);
getopts('g:f:', \%opts);
die("
Usage:   filters.pl file_to_be_processed output_file_left output_file_filtered_out\n
Options:
         -g INT breakpoint buffer
         -f FLOAT  maximum fraction to suffer\n;
") unless (@ARGV);

my $centro_read = 0;
my %CenStart;
my %CenEnd;
my $assembly_hit_cen = 0;
my $breakdancer_hit_cen = 0;
my $BreakpointBuffer = $opts{g};
my $maxFraction = $opts{f};
my %satelliteToChrStatements;

open FILE_out_left, ">", $ARGV[1] or die $!;
close FILE_out_left;
open FILE_out_filtered, ">", $ARGV[2] or die $!;
close FILE_out_filtered;
open FILE_out_left, ">>", $ARGV[1] or die $!;
open FILE_out_filtered, ">>", $ARGV[2] or die $!;

open(CONFIG,"<$ARGV[0]") || die "unable to open $ARGV[0]\n";
while(<CONFIG>){
	my $line = $_;
	my ($input, )=split("\t", $_);
	my ($chr1,$outer_start_,$inner_start_,$chr2,$inner_end_,$outer_end_) = ($input =~ /(\S+)\.(\d+)\.(\d+)\.(\S+)\.(\d+)\.(\d+)/);
#	print "$chr1\t$outer_start_\t$inner_start_\t$chr2\t$inner_end_\t$outer_end_\n";
		my $hit_cen = 0;
		if($centro_read == 0){
			$centro_read = 1;
			my $chr_cen;
			#read data table
			foreach $chr_cen (1..22,"X","Y") {
				my $db = "ucsc";
				my $user = "mgg_admin";
				my $password = "c\@nc3r";
				my $dataBase = "DBI:mysql:$db:mysql2";
				my $dbh = DBI->connect($dataBase, $user, $password) ||
				    die "ERROR: Could not connect to database: $! \n";
    
			    my $table = "chr$chr_cen"."_rmsk";
    			my $query = "SELECT genoStart, genoEnd
              	FROM $table
              	WHERE (repClass = 'Satellite' || repClass = 'Low_complexity'|| repClass = 'Simple_repeat')
                    && genoEnd >= ? && genoStart <= ?
              	ORDER BY genoStart";

			    $satelliteToChrStatements{$chr_cen} = $dbh->prepare($query) ||
				die "Could not prepare statement '$query': $DBI::errstr \n";
			}
		}
		# make a region around breakpoints
		my %regions = ();
		${$regions{$outer_start_-$BreakpointBuffer}}{$inner_start_+$BreakpointBuffer} = 1;
   	    ${$regions{$inner_end_-$BreakpointBuffer}}{$outer_end_+$BreakpointBuffer} = 1;
		my $satelliteRegion = 0;
	    my %outputRegion = ();
	    my $chrStart;
	    my $chrStop;
	    my $start;
	    my $stop;
	    my %satelliteCoords;
	    my $totalSatellite;
	    my $fraction;
	    foreach $start ( keys %regions ) {
			foreach $stop ( keys %{$regions{$start}} ) {
			    %satelliteCoords = ();
    			$totalSatellite = 0;
    	
    			$satelliteToChrStatements{$chr1}->execute($start, $stop) ||
					die "Could not execute statement for repeat masker table with (chr$chr1, $start, $stop): $DBI::errstr \n";
			    while ( ($chrStart, $chrStop) =  $satelliteToChrStatements{$chr1}->fetchrow_array() ) {
			    	my $start_last = ($chrStart > $start) ? $chrStart : $start;
			    	my $stop_last = ($chrStop < $stop) ? $chrStop : $stop;
					foreach ($start_last..$stop_last) { $satelliteCoords{$_} = 1; }
   				}
    			foreach ($start..$stop) {
					if ( defined $satelliteCoords{$_}  ) { $totalSatellite++; }
    			}
    			$fraction = $totalSatellite/($stop - $start + 1);
    			if ( $fraction > $maxFraction ) { 
    				$satelliteRegion = 1; 
    				${$outputRegion{$start}}{$stop} = $fraction;     				
    				print FILE_out_filtered "$line";
	    			last;
	    		}# print "$chr1\t$start_\t$end_\t$start\t$stop\t$fraction\n"; last;}
			}
			last if($satelliteRegion == 1);
    	}
#		$assembly_hit_cen ++ if($satelliteRegion == 1 && $fin =~ /assembled/);
		print FILE_out_left "$line" if($satelliteRegion == 0);
		$breakdancer_hit_cen ++ if($satelliteRegion == 1);    	
#    	next if($satelliteRegion == 1);
}

close FILE_out_left;
close FILE_out_filtered;
print "hit satellite and repeat mask and low complexity: $breakdancer_hit_cen\n";
