
package Genome::Model::Tools::Capture::GeneCoverage;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# GeneCoverage - Compare tumor versus normal models to find somatic events
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/09/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Capture::GeneCoverage {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		chromosome	=> { is => 'Text', doc => "Chromosome name", is_optional => 0 },
		start	=> { is => 'Text', doc => "Chromosome start position", is_optional => 0 },
		stop	=> { is => 'Text', doc => "Chromosome end position " , is_optional => 0},
		reference	=> { is => 'Text', doc => "Size of reference/variant contigs to generate", is_optional => 1, default => "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fasta" },		
		coverage_file	=> { is => 'Text', doc => "Q20 coverage file in chrom, pos, cov format", is_optional => 0 },
		window_size	=> { is => 'Text', doc => "Window size for genome partitioning" , is_optional => 0, default => 10},
		target_file	=> { is => 'Text', doc => "Targets file in BED format" , is_optional => 1},
		output_file	=> { is => 'Text', doc => "Output file" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Build a TSV file for plotting gene coverage"                 
}

sub help_synopsis {
    return <<EOS
This command builds a TSV file for plotting gene coverage
EXAMPLE:	gmt capture gene-coverage ...
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	my $chrom = $self->chromosome;
	my $chr_start = $self->start;
	my $chr_stop = $self->stop;
	my $reference = $self->reference;
	my $window_size = $self->window_size;
	
	
	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
		print OUTFILE "chrom\twindow_start\twindow_stop\twindow_gc\tcoverage\n";	
	}


	my %target_positions = ();

	if($self->target_file)
	{
		warn "Loading Targets...\n";
		%target_positions = load_targets($self->target_file);
	}


	## Load coverage ##
	warn "Loading Q20 coverage...\n";
	my %coverage = load_coverage_file($self->coverage_file);

	print "Segmenting genome...\n";	
	
	my $window_number = 0;
	
	my $window_start = $chr_start;
	my $window_stop = 0;
	
	while($window_start <= $chr_stop)
	{
		$window_number++;
		$window_stop = $window_start + $window_size - 1;

		## Get sequence and GC content ##
		my $window_seq = `samtools faidx $reference $chrom:$window_start-$window_stop | grep -v \">\"`;
		$window_seq =~ s/[^A-Za-z]//g if($window_seq);
		my $window_gc = gc_content($window_seq);

		## Calculate average coverage depth ##

		my $cov_sum = my $cov_num = my $coverage = 0;

		my $target_flag = 0;
		
		for(my $position = $window_start; $position <= $window_stop; $position++)
		{
			my $key = join("\t", $chrom, $position);
			
			$cov_sum += $coverage{$key} if($coverage{$key});
			$cov_num++;

			$target_flag = 1 if($target_positions{$key});
		}
		
		if($cov_num)
		{
			$coverage = $cov_sum / $cov_num;
			$coverage = sprintf("%d", $coverage);
		}
		
		print "$window_number\t$chrom\t$window_start\t$window_stop\t$window_seq\t$window_gc\t$coverage\n";

		if($self->output_file)
		{
			print OUTFILE "$chrom\t$window_start\t$window_stop\t$window_gc\t$coverage";
			
			if($self->target_file)
			{
				print OUTFILE "\t$target_flag";
			}
			
			print OUTFILE "\n";
		}
		
		$window_start += $window_size;
	}


	close(OUTFILE) if($self->output_file);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}






#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub load_coverage_file
{
	my $FileName = shift(@_);
	my %coverage = ();
	
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $position, $q20_cov) = split(/\t/, $line);
		
		my $key = join("\t", $chrom, $position);
		
		$coverage{$key} = $q20_cov;
	}
	
	close($input);

	return(%coverage);
}






#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub load_targets
{
	my $FileName = shift(@_);
	my %coverage = ();
	
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $start, $stop) = split(/\t/, $line);

		for(my $position = $start; $position <= $stop; $position++)
		{
			my $key = join("\t", $chrom, $position);
			$coverage{$key} = 1;
		}

		
	}
	
	close($input);

	return(%coverage);
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub gc_content
{
	my $seq = shift(@_);
	
	$seq = uc($seq);
	
	my @bases = split(//, $seq);
	
	my $num_bases = my $numGC = 0;
	
	foreach my $base (@bases)
	{
		$num_bases++;
		
		$numGC++ if($base eq "G" || $base eq "C");

	}
	
	if($num_bases)
	{
		my $gc_content = $numGC / $num_bases * 100;
		$gc_content = sprintf("%d", $gc_content);
		return($gc_content);
	}
	
	return(0);
}


1;

