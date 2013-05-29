package Genome::Model::Tools::Vcf::RareDelPower;     # rename this when you give the module file a different name <--

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



## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Vcf::RareDelPower {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		num_cases	=> { is => 'Text', doc => "Total number of cases in cohort" , is_optional => 0},
		num_controls	=> { is => 'Text', doc => "Total number of controls in cohort" , is_optional => 0},
		max_variant	=> { is => 'Text', doc => "Maximum number of cases-variant to test" , is_optional => 0},
		output_file	=> { is => 'Text', doc => "Output file for genes with FET p-value" , is_optional => 0},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compute possible FET p-values given a number of cases and controls"                 
}

sub help_synopsis {
    return <<EOS
This command computes possible FET p-values given a number of cases and controls
EXAMPLE:	gmt vcf rare-del-power --num-cases 100 --num-controls 100 --output-file fisher.power.tsv
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

	my $output_file = $self->output_file;
	my $num_cases = $self->num_cases;
	my $num_controls = $self->num_controls;
	my $max_variant = $self->max_variant;

	## Open the output file ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

	print OUTFILE "Cases variant (of $num_cases)\tControls_variant (of $num_controls)\n";

	for(my $controlCounter = 0; $controlCounter <= $max_variant; $controlCounter++)
	{
		print OUTFILE "\t$controlCounter";
	}
	print OUTFILE "\n";
	
	print "controls_wt\tcontrols_var\tcases_wt\tcases_var\tp_value\n";
	
	
	for(my $caseCounter = 0; $caseCounter <= $max_variant; $caseCounter++)
	{
		print OUTFILE "$caseCounter";
		for(my $controlCounter = 0; $controlCounter <= $max_variant; $controlCounter++)
		{
			my $cases_with_var = $caseCounter;
			my $cases_without_var = $num_cases - $cases_with_var;
			my $controls_with_var = $controlCounter;
			my $controls_without_var = $num_controls - $controls_with_var;
			
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
				
				print OUTFILE "\t$p_value";
#				print OUTFILE join("\t", $controls_without_var, $controls_with_var, $cases_without_var, $cases_with_var, $p_value) . "\n";
				print join("\t", $controls_without_var, $controls_with_var, $cases_without_var, $cases_with_var, $p_value) . "\n";
		}
		
		print OUTFILE "\n";
	}


				

	close(OUTFILE);

	foreach my $key (sort keys %stats)
	{
		print "$stats{$key}\t$key\n";
	}

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

1;
