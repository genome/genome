
package Genome::Model::Tools::Analysis::Mendelian::MpileupForVariants;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SearchRuns - Search the database for runs
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/01/2009 by D.K.
#	MODIFIED:	04/01/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it
use Genome::Model::Tools::Analysis::Helpers qw(
    code_to_genotype_returning_code
);

my $num_affected = my $affecteds_missing = my $unaffecteds_variant = my $affecteds_variant = my $affecteds_ambiguous = 0;
my %genotypes = ();

class Genome::Model::Tools::Analysis::Mendelian::MpileupForVariants {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file	=> { is => 'Text', doc => "List of variants to consider (annotation format)", is_optional => 0, is_input => 1},
		reference	=> { is => 'Text', doc => "Reference sequence to use", is_optional => 0, is_input => 1},
		bam_files	=> { is => 'Text', doc => "One or more BAM files for samples, comma-separated", is_optional => 0, is_input => 1},
		sample_names	=> { is => 'Text', doc => "Descriptive sample names for samples, comma-separated", is_optional => 0, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for QC result", is_optional => 1, is_input => 1},
		samtools_params	=> { is => 'Text', doc => "Parameters for SAMtools mpileup2cns", is_optional => 0, is_input => 1, default => '-q 1'},
		varscan_params	=> { is => 'Text', doc => "Parameters for VarScan mpileup2cns", is_optional => 0, is_input => 1, default => '--min-coverage 3 --min-var-freq 0.20 --p-value 0.10 --strand-filter 1'},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compiles a genotype table for all variants in a file"                 
}

sub help_synopsis {
    return <<EOS
This command compiles a genotype table for all variants in a file using SAMtools mpileup
EXAMPLE:	gmt analysis mendelian genotype-table --variant-file all.snvs --consensus-files myCons1,myCons2
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

	my $variant_file = $self->variant_file;

	my @bam_files = split(/\,/, $self->bam_files);
	my @sample_names = split(/\,/, $self->sample_names);
	
	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";		
	}

	my %stats = ();
	
		
	## Count the files of each type and print the header ##
	my $header = "chrom\tchr_start\tchr_stop\tref\tvar";
	foreach my $sample (@sample_names)
	{
		$header .= "\t" if($header);
		$header .= $sample . "\treads1\treads2\tfreq";
	}

	## Build BAM string for SAMtools mpileup command ##
	
	my $bam_string = "";
	
	foreach my $bam_file (@bam_files)
	{
		$bam_string .= "$bam_file ";
	}

	if($self->output_file)
	{
		print OUTFILE "$header\n";
	}

	## Print the variants ##

	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chrom, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);

		if(length($ref) > 1 || $var eq "-")
		{
			## Undo the adjustment made when formatting deletions for annotation.
			$chr_start--;
			$chr_stop--;
		}

		my $variant_type = "SNP";
		
		if(length($ref) > 1 || $var eq "-")
		{
			$variant_type = "DEL";
		}
		elsif(length($var) > 1 || $ref eq "-")
		{
			$variant_type = "INS";
		}

		$stats{'num_variants'}++;
					
		## Build the query ##
		
		my $query = $chrom . ":" . $chr_start . "-" . $chr_stop;
		
		my $cmd = "samtools mpileup -r " . $query . " -f " . $self->reference . " " . $self->samtools_params . " " . $bam_string . " | " . varscan_cmd() . "mpileup2cns " . $self->varscan_params;

		my $varscan_result = `$cmd | grep -v Chrom`;

		if($varscan_result)
		{
			print OUTFILE join("\t", $chrom, $chr_start, $chr_stop, $ref, $var) if($self->output_file);

			my @varscanContents = split(/\t/, $varscan_result);	
			my @varscan_calls = split(/\s+/, $varscanContents[10]);
			
			foreach my $varscan_call (@varscan_calls)
			{
				my ($cns, $cov, $reads1, $reads2, $freq, $pval) = split(/\:/, $varscan_call);
				print OUTFILE "\t" . join("\t", code_to_genotype_returning_code($cns), $reads1, $reads2, $freq) if($self->output_file);
			}
	
			print OUTFILE "\n" if($self->output_file);
	

		}

	
		

	
		
	}
	
	close($input);
	
	if($self->output_file)
	{
		close(OUTFILE);
	}
	
	print $stats{'num_variants'} . " variants\n";

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# The VarScan command
#
################################################################################################

sub varscan_cmd
{
#	return("java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.jar ");
	return("java -cp /gscuser/dkoboldt/Software/VarScan net.sf.varscan.VarScan ");
}

1;
