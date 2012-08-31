
package Genome::Model::Tools::Varscan::FilterSnps;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Varscan::FilterSnps	Process somatic pipeline output
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/29/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it


my $report_only = 0;
my %stats = ();

class Genome::Model::Tools::Varscan::FilterSnps {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		snp_file		=> { is => 'Text', doc => "File containing varscan calls, e.g. status.varscan.snp" , is_optional => 0, is_input => 1},
		indel_file	=> { is => 'Text', doc => "File containing varscan indel calls" , is_optional => 0, is_input => 1},
		output_file		=> { is => 'Text', doc => "Output file for filtered SNPs" , is_optional => 0, is_input => 1, is_output => 1},
		max_indel_proximity	=> { is => 'Text', doc => "If variant is this number of bp from an indel, filter it [3]" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Filter Varscan SNPs by indels and by clusters"                 
}

sub help_synopsis {
    return <<EOS
Filters Varscan SNP calls near indels and clusters of SNPs
EXAMPLE:	gmt varscan filter-snps ...
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

	## Get required parameters ##
	my $snp_file = $self->snp_file;
	my $indel_file = $self->indel_file;
	my $output_file = $self->output_file;
	my $max_indel_proximity = 3;
	$max_indel_proximity = $self->max_indel_proximity if($self->max_indel_proximity);

	if(!(-e $snp_file))
	{
		die "ERROR: SNP file $snp_file not found\n";
	}

	if(!(-e $indel_file))
	{
		die "ERROR: Indel file $indel_file not found\n";
	}

	## OPen outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	


	## Find SNP clusters ##
	my %snp_clusters = find_snp_clusters($snp_file);


	## Find indels ##
	my %indel_positions = load_indels($indel_file, $max_indel_proximity);

	## Parse the SNPs file ##

	my $input = new FileHandle ($snp_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);
		
		if(($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name"))
		{
			## Skip file header ##
			print OUTFILE "$line\n";
			
		}
		else
		{
			my $chrom = $lineContents[0];
			my $position = $lineContents[1];
			$stats{'total_snps'}++;

			if($snp_clusters{"$chrom\t$position"})
			{
				$stats{'cluster_snps'}++;		
			}
			elsif($indel_positions{"$chrom\t$position"})
			{
				$stats{'indel_snps'}++;
			}
			else
			{
				print OUTFILE "$line\n";
			}
		}
	}
	
	close($input);


	print "$stats{'total_snps'} total SNPs\n";
	print "$stats{'cluster_snps'} removed due to clustering\n";
	print "$stats{'indel_snps'} removed within $max_indel_proximity bp of indels\n";
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# Process results - filter variants by type and into high/low confidence
#
################################################################################################

sub find_snp_clusters
{
	my $variants_file = shift(@_);

	my %cluster_snps = ();

	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	
	my $window_start = my $window_stop = 0;

	my $chrom1 = my $chrom2 = my $chrom3 = "";
	my $pos1 = my $pos2 = my $pos3 = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);
		
		if(($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name"))
		{
			## Skip file header ##
		}
		else
		{
			my $chrom = $lineContents[0];
			my $position = $lineContents[1];
			
			## This SNP becomes number three ##
			$chrom3 = $chrom;
			$pos3 = $position;
			
			## Check for a SNP cluster ##

			if($chrom1 && $chrom2 && $chrom3 && $chrom1 eq $chrom2 && $chrom2 eq $chrom3)
			{
				my $this_window = $pos3 - $pos1 + 1;
				if($this_window <= 10)
				{
					$cluster_snps{"$chrom1\t$pos1"} = 1;
					$cluster_snps{"$chrom2\t$pos2"} = 1;
					$cluster_snps{"$chrom3\t$pos3"} = 1;
				}
			}
			
			## Move previous SNPs forward ##
			$chrom1 = $chrom2;
			$pos1 = $pos2;
			$chrom2 = $chrom3;
			$pos2 = $pos3;
		}
	}
	
	close($input);
	
	return(%cluster_snps);
	
}







################################################################################################
# Process results - filter variants by type and into high/low confidence
#
################################################################################################

sub load_indels
{
	my $variants_file = shift(@_);
	my $max_proximity = shift(@_);

	my %indels = ();

	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);
		
		if(($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name"))
		{
			## Skip file header ##
		}
		else
		{
			my $chrom = $lineContents[0];
			my $position = $lineContents[1];
			
			my $chr_start = my $chr_stop = $position;
			
			if($lineContents[2] && $lineContents[2] =~ /[0-9]/ && $lineContents[2] > $lineContents[1])
			{
				$chr_stop = $lineContents[2];
			}
			
			if($chrom && $chr_start && $chr_stop && $chr_stop >= $chr_start)
			{
				## Count the indel ##
				
				for($position = $chr_start - $max_proximity; $position <= $chr_stop + $max_proximity; $position++)
				{
					my $key = "$chrom\t$position";
					$indels{$key} = 1;
				}
			}
		}
	}
	
	close($input);
	
	return(%indels);
	
}



1;

