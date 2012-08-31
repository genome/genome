
package Genome::Model::Tools::SnpArray::RecurrentEvents;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# RecurrentEvents - merges adjoining segments of similar copy number; distinguishes amplifications and deletions
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

my %stats = ();

class Genome::Model::Tools::SnpArray::RecurrentEvents {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sample_events_file	=> { is => 'Text', doc => "Segments with p-values from running CBS on data", is_optional => 0, is_input => 1 },		
		output_file	=> { is => 'Text', doc => "Output file for recurrent event counts", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merges adjoining segments of similar copy number"                 
}

sub help_synopsis {
    return <<EOS
This command merges merges adjoining CBS segments of similar copy number and distinguishes amplifications and deletions
EXAMPLE:	gmt snp-array recurrent-events --sample-events-file Sample-CopyNumberEvent-Files.tsv
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

sub execute
{                               # replace with real execution logic.
	my $self = shift;

	## Get required parameters ##
	my $sample_events_file = $self->sample_events_file;
	my $output_file = $self->output_file;

	my $num_samples = `cat $sample_events_file | wc -l`;
	chomp($num_samples);
	
	print "$num_samples samples\n";

	open(OUTFILE, ">$output_file") if($self->output_file);
	print OUTFILE "chrom_arm\tamps\t%cases_amp\tdels\t%cases_del\n" if($self->output_file);

	my %amplifications = run_recurrence($self, "large-scale", "amplification");
	my %deletions = run_recurrence($self, "large-scale", "deletion");
	
	for(my $chrCounter = 1; $chrCounter <= 23; $chrCounter++)
	{
		my $chrom = $chrCounter;
		$chrom = "X" if($chrCounter == 23);

		## Get the P arm ##
		my $arm = "p";
		my $p_amps = 0;
		$p_amps = $amplifications{$chrom . $arm} if($amplifications{$chrom . $arm});
		my $p_amp_freq = sprintf("%.2f", $p_amps / $num_samples * 100) . '%';
		
		my $p_dels = 0;
		$p_dels = $deletions{$chrom . $arm} if($deletions{$chrom . $arm});
		my $p_del_freq = sprintf("%.2f", $p_dels / $num_samples * 100) . '%';
		print join("\t", $chrom . $arm, $p_amps, $p_amp_freq, $p_dels, $p_del_freq) . "\n";
		
		## Get the Q arm ##
		$arm = "q";
		my $q_amps = 0;
		$q_amps = $amplifications{$chrom . $arm} if($amplifications{$chrom . $arm});
		my $q_amp_freq = sprintf("%.2f", $q_amps / $num_samples * 100) . '%';
		
		my $q_dels = 0;
		$q_dels = $deletions{$chrom . $arm} if($deletions{$chrom . $arm});
		my $q_del_freq = sprintf("%.2f", $q_dels / $num_samples * 100) . '%';

		print join("\t", $chrom . $arm, $q_amps, $q_amp_freq, $q_dels, $q_del_freq) . "\n";
		print OUTFILE join("\t", $chrom . $arm, $q_amps, $q_amp_freq, $q_dels, $q_del_freq) . "\n" if($self->output_file);
	}
	
	close(OUTFILE) if($self->output_file);
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub run_recurrence
{
	my $self = shift(@_);
	my $event_size = shift(@_);
	my $event_type = shift(@_);

	## Get required parameters ##
	my $sample_events_file = $self->sample_events_file;


	my %arm_event_counts = ();
	
	my $input = new FileHandle ($sample_events_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if(1) #$lineCounter <= 20)#($lineCounter > 0)#
		{
			my ($sample_name, $file_name) = split(/\t/, $line);

			my $events = load_events($file_name, $event_size, $event_type);

			## Save the positions of this sample's events ##

			my @events = split(/\n/, $events);
			my $num_events = @events;

			my %event_counted = ();
			
			foreach my $event (@events)
			{
				if(!$event_counted{$event})
				{
					$arm_event_counts{$event}++;
					$event_counted{$event} = 1;
				}

			}
			
#			print "$chrom\t$sample_name\t$num_events\n";				
		}


	}
	
	close($input);
	
	
	return(%arm_event_counts);

}




sub numerically
{
	$a <=> $b;
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_events
{
	my $FileName = shift(@_);
#	my $desired_chrom = shift(@_);
	my $desired_class = shift(@_);
	my $desired_type = shift(@_);
	my $events = "";

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if($lineCounter > 1)
		{
			my ($chrom, $chr_start, $chr_stop, $seg_mean, $num_segments, $num_markers, $p_value, $event_type, $event_size, $size_class, $chrom_arm) = split(/\t/, $_);

#			if($chrom eq $desired_chrom)
#			{
				if($size_class eq $desired_class)
				{
					if($event_type eq $desired_type)
					{
						$events .= $chrom_arm . "\n";
					}
				}
#			}
		}


	}
	
	close($input);	
	
	return($events);
}



1;


