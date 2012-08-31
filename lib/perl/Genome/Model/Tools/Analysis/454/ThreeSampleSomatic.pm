
package Genome::Model::Tools::Analysis::454::ThreeSampleSomatic;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ThreeSampleSomatic - Runs the Varscan somatic pipeline on matched tumor-normal data
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/25/2009 by D.K.
#	MODIFIED:	02/25/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::454::ThreeSampleSomatic {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sample_info_file	=> { is => 'Text', doc => "Tab-delimited file of sample info from database. Col1=sample Col2=normal|tumor|relapse" },
		output_dir	=> { is => 'Text', doc => "Output directory for 454 data. Will create somatic_pipeline in each tumor sample dir" },
		readcounts_file		=> { is => 'Text', doc => "Path within sample directory to varscan readcounts file" },
		tier1_file		=> { is => 'Text', doc => "Path to tier 1 variants in relapse/tumor samples" },
		variant_file		=> { is => 'Text', doc => "Path to list of variants with annotation" },
		reference		=> { is => 'Text', doc => "Reference sequence [default=Hs36 ssaha2]", is_optional => 1 },
		skip_if_output_present	=> { is => 'Text', doc => "Skip if output present", is_optional => 1 },		
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compares tumor and relapse mutations from capture somatic pipeline"                 
}

sub help_synopsis {
    return <<EOS
This command runs the Varscan somatic pipeline on matched normal-tumor samples
EXAMPLE:	gmt analysis 454 three-sample-somatic --sample-info-file Sample-Info-From-DB.tsv --output-dir data --variant-file all-tier1-rna-snvs.tsv --readcounts-file varscan_somatic/varScan.output.snp.formatted.filter.novel.annotation
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
	my $sample_info_file = $self->sample_info_file;
	my $output_dir = $self->output_dir;
	my $readcounts_file = $self->readcounts_file;

	mkdir($output_dir . "/three_sample_somatic") if(!(-d "$output_dir/three_sample_somatic"));

	if(!(-e $sample_info_file))
	{
		die "Error: Samples file not found!\n";
	}

	my %patient_ids = ();
	my %patient_samples = ();

	my $input = new FileHandle ($sample_info_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my ($sample_name, $sample_type) = split(/\t/, $line);
		
		if($sample_type eq "normal" || $sample_type eq "tumor" || $sample_type eq "relapse")
		{
			my @sampleContents = split(/\-/, $sample_name);
			my $numContents = @sampleContents;
			my $patient_id = "";
			
			for(my $colCounter = 0; $colCounter < ($numContents - 1); $colCounter++)
			{
				$patient_id .= "-" if($patient_id);
				$patient_id .= "$sampleContents[$colCounter]";				
			}
			
			$patient_ids{$patient_id}++;
			my $key = "$patient_id\t$sample_type";
			
			warn "Warning: Multiple samples of type $sample_type for patient $patient_id; last one will be used\n" if($patient_samples{$key});

			$patient_samples{$key} = $sample_name;
			
			#print "$sample_name\t$patient_id\t$sample_type\n";
		}
		else
		{
			warn "Line $lineCounter: Unrecognized sample type $sample_type ignored\n";
		}

	}
	
	close($input);


	## Go through each patient ##

	my $output_file = $output_dir . "/three_sample_somatic/Tier1-Mutation-Readcounts.readcounts.tsv";
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

	foreach my $patient_id (keys %patient_ids)
	{
		my $normal_sample = my $tumor_sample = my $relapse_sample = "-";
		
		## Open the output file ##
#		my $output_file = $output_dir . "/three_sample_somatic/" . $patient_id . ".readcounts.tsv";
#		open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
		
		my $num_tumor_mutations = my $num_relapse_mutations = my $num_shared_mutations = 0;
		
		$normal_sample = $patient_samples{"$patient_id\tnormal"} if ($patient_samples{"$patient_id\tnormal"});
		$tumor_sample = $patient_samples{"$patient_id\ttumor"} if ($patient_samples{"$patient_id\ttumor"});
		$relapse_sample = $patient_samples{"$patient_id\trelapse"} if($patient_samples{"$patient_id\trelapse"});

		## Get the tumor mutations ##

		my %normal_readcounts = load_mutations($self, $normal_sample) if($normal_sample ne "-");
		my %tumor_readcounts = load_mutations($self, $tumor_sample) if($tumor_sample ne "-");
		my %relapse_readcounts = load_mutations($self, $relapse_sample) if($relapse_sample ne "-");

		## Get tier 1 calls for normal and relapse ##
		my %tumor_tier1_calls = load_tier1_calls($self, $tumor_sample) if($tumor_sample ne "-");
		my %relapse_tier1_calls = load_tier1_calls($self, $relapse_sample) if($relapse_sample ne "-");		

		## Parse the variants file

		my $input = new FileHandle ($self->variant_file);
		my $lineCounter = 0;
		
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;		
		
			my @lineContents = split(/\t/, $line);
			my ($chrom, $chr_start, $chr_stop, $ref, $var, $var_type, $gene_name) = split(/\t/, $line);

			my $key = join("\t", $chrom, $chr_start);

			## Get normal results ##
			
			my $normal_results = "-\t-\t-";
			
			if($normal_readcounts{$key})
			{
				$normal_results = get_rc_results($ref, $var, $normal_readcounts{$key});	
			}

			## Get tumor results ##
			
			my $tumor_results = "-\t-\t-";
			
			if($tumor_readcounts{$key})
			{
				$tumor_results = get_rc_results($ref, $var, $tumor_readcounts{$key});	
			}

			## Get relapse results ##
			
			my $relapse_results = "";
			## Only set defaults if relapse was found ##
			$relapse_results = "-\t-\t-" if($patient_samples{"$patient_id\trelapse"});
			
			if($relapse_readcounts{$key})
			{
				$relapse_results = get_rc_results($ref, $var, $relapse_readcounts{$key});	
			}

			if($tumor_tier1_calls{$key} || $relapse_tier1_calls{$key})
			{
				print OUTFILE join("\t", $line, $normal_sample, $tumor_sample, $relapse_sample, $normal_results, $tumor_results, $relapse_results) . "\n";
				print join("\t", $key, $normal_results, $tumor_results, $relapse_results) . "\n";				
			}

		}
		
		close($input);


#		close(OUTFILE);


#		print join("\t", $patient_id, $normal_sample, $tumor_sample, $relapse_sample, $num_tumor_mutations, $num_relapse_mutations, $num_shared_mutations) . "\n";
	}
	
	close(OUTFILE);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}





################################################################################################
# Get RC results - get the readcounts and var freq for this allele
#
################################################################################################

sub get_rc_results
{
	my ($ref_allele, $variant_allele, $rc_line) = @_;


	my $reads1 = my $reads2 = my $var_freq = "-";
	
	my @lineContents = split(/\t/, $rc_line);
	my $numContents = @lineContents;
	
	if($numContents > 4)
	{
		my $ref_counts = $lineContents[5];
		my @refCounts = split(/\:/, $ref_counts);
		$reads1 = $refCounts[1];
		$reads2 = 0;
	
		my $var_counts = "";
		my $other_alleles = "";
	
		for(my $colCounter = 6; $colCounter < $numContents; $colCounter++)
		{
			(my $allele) = split(/\:/, $lineContents[$colCounter]);
			if($allele eq $variant_allele)
			{
				$var_counts = $lineContents[$colCounter];
				my @varCounts = split(/\:/, $var_counts);
				$reads2 = $varCounts[1];
				$var_freq = sprintf("%.2f", ($reads2 / ($reads1 + $reads2) * 100)) . '%';
			}
		}		

		$var_freq = sprintf("%.2f", ($reads2 / ($reads1 + $reads2) * 100)) . '%';

		return(join("\t", $reads1, $reads2, $var_freq));
	}
	
	return("-\t-\t-");
}



################################################################################################
# Load mutations - find and parse the mutations file 
#
################################################################################################

sub load_mutations
{
	my $self = shift;
	my $sample_name = shift;

	my $output_dir = $self->output_dir;
	my $readcounts_file = $self->readcounts_file;

	my %mutations = ();

	## Determine path to file ##
	
	my $mutations_file = "$output_dir/$sample_name/$readcounts_file";
	
	if(-e $mutations_file)
	{
#		print "Parsing $mutations_file for $sample_name\n";
		my $input = new FileHandle ($mutations_file);
		my $lineCounter = 0;
		
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;		
		
			my @lineContents = split(/\t/, $line);
			my ($chrom, $chr_start, $chr_stop, $ref, $var, $var_type, $gene_name) = split(/\t/, $line);

			my $key = join("\t", $chrom, $chr_start);
			$mutations{$key} = $line;
		}
		
		close($input);


		
	}
	else
	{
#		warn "Warning: Mutations file $mutations_file not found\n";
	}



	return(%mutations);

}






################################################################################################
# Load mutations - find and parse the mutations file 
#
################################################################################################

sub load_tier1_calls
{
	my $self = shift;
	my $sample_name = shift;

	my $output_dir = $self->output_dir;
	my $tier1_file = $self->tier1_file;

	my %mutations = ();

	## Determine path to file ##
	
	my $mutations_file = "$output_dir/$sample_name/$tier1_file";
	
	if(-e $mutations_file)
	{
		my $input = new FileHandle ($mutations_file);
		my $lineCounter = 0;
		
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;		
		
			my @lineContents = split(/\t/, $line);
			my ($chrom, $chr_start, $chr_stop, $ref, $var, $var_type, $gene_name) = split(/\t/, $line);

			my $key = join("\t", $chrom, $chr_start);
			$mutations{$key} = $line;
		}
		
		close($input);


		
	}
	else
	{
#		warn "Warning: Mutations file $mutations_file not found\n";
	}



	return(%mutations);

}


1;

