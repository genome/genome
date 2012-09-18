package Genome::Model::Tools::Vcf::RareDelTest;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MutationRate - Calculate the mutation rate (per megabase) given a list of mutations (e.g. tier1 SNVs) and a set of regions (e.g. coding space)
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	04/22/2011 by D.K.
#	MODIFIED:	04/22/2011 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it


## Pre-define a ranking system for VEP annotation, where higher = more severe ##


## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Vcf::RareDelTest {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		gene_file	=> { is => 'Text', doc => "Input rare-deleterious gene table" , is_optional => 0},
		variant_file	=> { is => 'Text', doc => "Sample Phenotype File with column named phenotype and values 0 or 1" , is_optional => 0},
		output_file	=> { is => 'Text', doc => "Output file for genes with FET p-value" , is_optional => 0},
		output_significant	=> { is => 'Text', doc => "Output file for genes with significant FET p-value" , is_optional => 0},
		output_details	=> { is => 'Text', doc => "Output file for details of significant genes" , is_optional => 1},
		p_value_threshold	=> { is => 'Text', doc => "Default p-value threshold to report details for genes" , is_optional => 0, default => 0.05},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Runs a Fisher's Exact Test on a rare deleterious table of genes"                 
}

sub help_synopsis {
    return <<EOS
This command produces a rare/deleterious table by gene from VCF and VEP files
EXAMPLE:	gmt vcf rare-del-test --gene-file my.deltable.tsv --variant-file my.deltable.variants.tsv --output-file my.deltable.fisher.tsv
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

        my $gene_file = $self->gene_file;
	my $variant_file = $self->variant_file;
	my $output_file = $self->output_file;
	my $output_significant = $self->output_significant;
	my $p_threshold = $self->p_value_threshold;

	## Open the output file ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	open(SIGNIFICANT, ">$output_significant") or die "Can't open outfile: $!\n";

	## Parse the file ##

	my $input = new FileHandle ($gene_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if($lineCounter == 1)
		{
			print OUTFILE "$line\tfet_p_value\n";
			print SIGNIFICANT "$line\tfet_p_value\n";
		}
		else
		{
			my ($gene, $rare_del_vars, $control_variants, $case_variants, $controls_without_var, $controls_with_var, $cases_without_var, $cases_with_var, $pct_controls, $pct_cases) = split(/\t/, $line);
			
#			if($lineCounter < 20)
#			{
				open(SCRIPT, ">temp.R");
				print SCRIPT qq{
deltable <- matrix(c($controls_without_var, $controls_with_var, $cases_without_var, $cases_with_var), nr=2, dimnames=list(c("Neut", "Delet"), c("Control", "Case")))
ftest <- fisher.test(deltable)
write(ftest\$p.value, file="temp.R.out", append=FALSE)
				};
				close(SCRIPT);

				system("R --no-save < temp.R 1>/dev/null 2>/dev/null");

				my $p_value = `cat temp.R.out`;
				chomp($p_value);
				$p_value = "NA" if(length($p_value) < 1);
				
				print OUTFILE "$line\t$p_value\n";
				warn join("\t", $lineCounter, $gene, $controls_with_var, $cases_with_var, $p_value) . "\n";

				if($p_value ne "NA" && $p_value < $p_threshold)
				{
					print SIGNIFICANT "$line\t$p_value\n";					
				}

#			}
		}
	}
	
	close($input);

	close(SIGNIFICANT);
	close(OUTFILE);

	foreach my $key (sort keys %stats)
	{
		print "$stats{$key}\t$key\n";
	}

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

1;
