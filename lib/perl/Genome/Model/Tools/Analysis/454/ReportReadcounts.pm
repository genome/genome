
package Genome::Model::Tools::Analysis::454::ReportReadcounts;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ReportReadcounts - Report read counts at a list of positions across a set of 454 samples
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/25/2009 by D.K.
#	MODIFIED:	02/25/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

my $ref_seq = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa';
my $ref_index = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa.fai';

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::454::ReportReadcounts {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		samples_file	=> { is => 'Text', doc => "Tab-delimited file of sample and SFF file(s)" },
		variants_file	=> { is => 'Text', doc => "Tab-delimited file of variant positions" },		
		output_dir	=> { is => 'Text', doc => "Output directory" },
		aligner		=> { is => 'Text', doc => "Aligner that was used." },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Report read counts at a list of positions across a set of samples"                 
}

sub help_synopsis {
    return <<EOS
This command converts SAM files to BAM files
EXAMPLE:	gmt analysis 454 report-readcounts --samples-file data/samples.tsv --output-dir data --aligner ssaha2 --variants-file targets/myVariants.tsv
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
	my $samples_file = $self->samples_file;
	my $variants_file = $self->variants_file;
	my $output_dir = $self->output_dir;
	my $aligner = $self->aligner;

	if(!(-e $samples_file))
	{
		die "Error: Samples file not found!\n";
	}

	if(!(-e $variants_file))
	{
		die "Error: Variants file not found!\n";
	}


	my %samples = ();
	my %sample_readcounts = ();

	## Parse sample list ##

	my $input = new FileHandle ($samples_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my @lineContents = split(/\t/, $line);			
		my $sample_name = $lineContents[0];
		
		$samples{$sample_name} = 1;
		
		## Create sample output directories ##
		my $sample_output_dir = $output_dir . "/" . $sample_name;
		my $aligner_output_dir = $sample_output_dir . "/" . $aligner . "_out";
		my $aligner_scripts_dir = $sample_output_dir . "/" . $aligner . "_out/scripts";		
		my $aligner_output_file = "$aligner_output_dir/$sample_name.$aligner.sam";
		my $bam_file = "$aligner_output_dir/$sample_name.$aligner.bam";
		my $pileup_file = $bam_file . ".roi.pileup";
		my $cns_file = $bam_file . ".roi.pileup.cns";
		
		if(-e $pileup_file)
		{
			if(-e $cns_file)
			{
				print "$cns_file\n";
				
				my %readcounts = get_readcounts($variants_file, $cns_file);
				
				foreach my $key (keys %readcounts)
				{
					## Make new key with sample name ##
					
					my $sample_key = "$sample_name\t$key";
					
					$sample_readcounts{$sample_key} = $readcounts{$key};
				}				
			}

		}
		else
		{
			die "Pileup file $pileup_file not found!\n";
		}
		
	}
	
	close($input);


	## Print the header ##
	
	print "chrom\tposition";

	foreach my $sample_name (sort keys %samples)
	{
		print "\t$sample_name\t$sample_name\t$sample_name";
	}

	print "\n";

	## Parse the positions file one last time ##

	my %variant_positions = load_positions($variants_file);
	
	foreach my $key (sort keys %variant_positions)
	{
#		print "$key";

		foreach my $sample_name (sort keys %samples)
		{
			my $sample_key = "$sample_name\t$key";
			my $result = "N\t0\tN\t0";
			$result = $sample_readcounts{$sample_key} if($sample_readcounts{$sample_key});
			
			print "$sample_key\t$result\n";
			
#			print "\t$result";
		}
		
#		print "\n";
	}

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Execute - the main program logic
#
################################################################################################

sub get_readcounts
{
	(my $variants_file, my $cns_file) = @_;

	## Load variant positions##
	
	my %variant_positions = load_positions($variants_file);
	my %readcounts = ();

	my $input = new FileHandle ($cns_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;	
		
		my @lineContents = split(/\t/, $line);
		my $chrom = $lineContents[0];
		my $position = $lineContents[1];
		my $ref_base = $lineContents[2];
		my $coverage = $lineContents[4];
		my $reads1 = $lineContents[5];
		my $allele2 = $lineContents[8];
		my $reads2 = $lineContents[9];
		$allele2 = "" if(!$allele2);
		$reads2 = 0 if(!$reads2);

		my $key = "$chrom\t$position";
		
		## Get the read counts ##
		if($variant_positions{$key})
		{
			$readcounts{$key} = "$reads1\t$allele2\t$reads2";
		}
	}
	
	close($input);

	return(%readcounts);
}






################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_positions
{
	my $variants_file = shift(@_);

	my %positions = ();

	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;			
		
		(my $chrom, my $position) = split(/\t/, $line);
		my $key = "$chrom\t$position";
		
		$positions{$key} = $line;
	}
	
	close($input);

#	print "$lineCounter positions loaded\n";

	return(%positions);                       # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

