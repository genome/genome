
package Genome::Model::Tools::Capture::WiggleFromQ20;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# WiggleFromQ20 - Converts a pileup file to a three column file of chromosome, position, and bases with q>20.
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/10/2010 by D.K.
#	MODIFIED:	02/10/2010 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Capture::WiggleFromQ20 {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		normal_file	=> { is => 'Text', doc => "Q20 Coverage for Normal", is_optional => 0, is_input => 1 },
		tumor_file	=> { is => 'Text', doc => "Q20 Coverage for Tumor", is_optional => 0, is_input => 1 },
		regions_file	=> { is => 'Text', doc => "Tab-delimited list of regions", is_optional => 0, is_input => 1 },
		min_depth_normal	=> { is => 'Text', doc => "Minimum Q20 depth for Normal [6]", is_optional => 1, is_input => 1 },
		min_depth_tumor	=> { is => 'Text', doc => "Minimum Q20 depth for Tumor [8]", is_optional => 1, is_input => 1 },
		output_file     => { is => 'Text', doc => "Output file to receive per-base qual>min coverage", is_optional => 0, is_input => 1, is_output => 1 },
		gzip_after	=> { is => 'Text', doc => "If set to 1, compress the file after building [0]", is_optional => 1, is_input => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Builds a wiggle coverage file for a list of ROIs"                 
}

sub help_synopsis {
    return <<EOS
This command builds a wiggle coverage file for a list of ROIs
EXAMPLE:	gmt capture wiggle-from-q20 --normal-file [normal.pileup.q20] --tumor-file [tumor.pileup.q20] --regions-file [targets.tsv] --output-file [patient.wiggle]
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
	my $normal_file = $self->normal_file;
	my $tumor_file = $self->tumor_file;
	my $regions_file = $self->regions_file;
	my $output_file = $self->output_file;
	
	my $min_depth_normal = 6;
	my $min_depth_tumor = 8;
	
	$min_depth_normal = $self->min_depth_normal if($self->min_depth_normal);
	$min_depth_tumor = $self->min_depth_tumor if($self->min_depth_tumor);
	
	
	## Load the normal coverage ##
	
	print "Loading normal coverage...\n";
	my %normal_coverage = parse_q20_file($normal_file);
	print "Loading tumor coverage...\n";
	my %tumor_coverage = parse_q20_file($tumor_file);
	print "Parsing regions file...\n";	

	## Open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	
	my %stats = ();
	$stats{'bases'} = $stats{'covered'} = $stats{'not_covered'} = 0;
	
	## Parse the regions ##

	my $input = new FileHandle ($regions_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		(my $chrom, my $chr_start, my $chr_stop, my $region_name) = split(/\t/, $line);

		## Print wiggle header ##
		print OUTFILE "fixedStep chrom=chr$chrom start=$chr_start step=1\n";			

		for(my $position = $chr_start; $position <= $chr_stop; $position++)
		{
			my $key = "$chrom\t$position";
			$stats{'bases'}++;

			## Determine if coverage is met ##
			
			if($normal_coverage{$key} && $tumor_coverage{$key} && $normal_coverage{$key} >= $min_depth_normal && $tumor_coverage{$key} >= $min_depth_tumor)
			{				
				print OUTFILE "1\n";
				$stats{'covered'}++;
			}
			else
			{
				print OUTFILE "0\n";
				$stats{'not_covered'}++;
			}
		}
	}

	close($input);


	close(OUTFILE);

	print $stats{'bases'} . " bases in ROI\n";
	print $stats{'covered'} . " bases covered >= " . $min_depth_normal . "x in normal and >= " . $min_depth_tumor . "x in tumor\n";
	print $stats{'not_covered'} . " bases NOT covered >= " . $min_depth_normal . "x in normal and >= " . $min_depth_tumor . "x in tumor\n";

	if($self->gzip_after)
	{
		print "Compressing $output_file...\n";
		system("gzip $output_file"); 
	}
	
	return 1;
}








#############################################################
# parse_file - takes input file and parses it
#
#############################################################

sub parse_q20_file
{
	my $filename = shift(@_);
	my %coverage_by_position = ();

	## Parse the positions ##
	
	my $input = new FileHandle ($filename);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chrom, my $position, my $coverage) = split(/\t/, $line);

		$chrom =~ s/chr//;
		$chrom = "MT" if($chrom eq "M");

		my $key = "$chrom\t$position";
		$coverage_by_position{$key} = $coverage;
		
#		return(%coverage_by_position) if($lineCounter > 10000);

	}
	
	close($input);

	return(%coverage_by_position);
}




1;

