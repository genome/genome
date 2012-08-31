
package Genome::Model::Tools::Capture::MergeAdaptedIndels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MergeAdaptedIndels - Merge Indel Calls from Varscan and Somatic Sniper
#					
#	AUTHOR:		Will Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	1/13/2009 by W.S.
#	MODIFIED:	1/13/2009 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Capture::MergeAdaptedIndels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		glf_file	=> { is => 'Text', doc => "Somatic Sniper Adapted Indel Input File", is_optional => 0, is_input => 1 },
		varscan_file	=> { is => 'Text', doc => "Varscan Adapted Indel Input File", is_optional => 0, is_input => 1 },
		gatk_file	=> { is => 'Text', doc => "GATK Adapted Indel Input File", is_optional => 1, is_input => 1 },
		gatk_unified_file	=> { is => 'Text', doc => "GATK Unified Genotyper Adapted Indel Input File", is_optional => 1, is_input => 1 },
		output_file	=> { is => 'Text', doc => "Merged Indel Output File" , is_optional => 0, is_input => 1, is_output => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merge Indel Calls from Varscan and Somatic Sniper"                 
}

sub help_synopsis {
    return <<EOS
This file was created to merge Indel Calls from Varscan and Somatic Sniper.
This requires inputs of ADAPTED files with chr pos pos ref var as first 5 columns.
EXAMPLE:	gmt capture merge-adapted-indels ...
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
	my $glf_file = $self->glf_file;
	my $varscan_file = $self->varscan_file;
	my $gatk_file;
	if ($self->gatk_file) {
		$gatk_file = $self->gatk_file;
	}
	my $gatk_unified_file;
	if ($self->gatk_unified_file) {
		$gatk_unified_file = $self->gatk_unified_file;
	}
	my $output_file = $self->output_file;

	my %stats = ();

	## Load the indel sets ##
	
	my %glf_indels = load_indels($glf_file);
	my %varscan_indels = load_indels($varscan_file);
	my %gatk_indels;
	my %gatk_IndelGenotyperV2_indels;
	my %gatk_unified_indels;
	if ($self->gatk_file && $self->gatk_unified_file) {
		%gatk_IndelGenotyperV2_indels = load_indels($gatk_file);
		%gatk_unified_indels = load_indels($gatk_unified_file);
		foreach my $keys (sort keys %gatk_unified_indels) {
			$gatk_indels{$keys} = $gatk_unified_indels{$keys};
		}
		foreach my $keys (sort keys %gatk_IndelGenotyperV2_indels) {
			$gatk_indels{$keys} = $gatk_IndelGenotyperV2_indels{$keys};
		}
	}
	elsif ($self->gatk_file) {
		%gatk_indels = load_indels($gatk_file);
	}
	elsif ($self->gatk_unified_file) {
		%gatk_indels = load_indels($gatk_unified_file);
	}

	## Build a list of all unique indel keys ##

	my %indel_keys = ();

	foreach my $key (keys %glf_indels)
	{
		$indel_keys{$key} = "sniper-only";
	}

	foreach my $key (keys %varscan_indels)
	{
		if($indel_keys{$key})
		{
			$indel_keys{$key} = "sniper-varscan-shared";
		}
		else
		{
			$indel_keys{$key} = "varscan-only";
		}
	}
	if ($self->gatk_file) {
		foreach my $key (keys %gatk_indels)
		{
			if($indel_keys{$key})
			{
				if($indel_keys{$key} eq "sniper-varscan-shared")
				{
					$indel_keys{$key} = "tri-shared";
				}
				elsif($indel_keys{$key} eq "varscan-only")
				{
					$indel_keys{$key} = "varscan-gatk-shared";
				}
				elsif($indel_keys{$key} eq "sniper-only")
				{
					$indel_keys{$key} = "sniper-gatk-shared";
				}
			}
			else
			{
				$indel_keys{$key} = "gatk-only";
			}
		}
	}

	## Open output files ##

	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	open(SHARED, ">$output_file.shared") or die "Can't open outfile: $!\n";
	open(SNIPER, ">$output_file.sniper-only") or die "Can't open outfile: $!\n";
	open(VARSCAN, ">$output_file.varscan-only") or die "Can't open outfile: $!\n";


	if ($self->gatk_file) {
		open(GATK, ">$output_file.gatk-only") or die "Can't open outfile: $!\n";
	}


	## Go thru all unique keys and handle them ##

	foreach my $indel_key (keys %indel_keys)
	{
		$stats{'total'}++;
		
		if($indel_keys{$indel_key} eq "tri-shared")
		{
			$stats{'tri-shared'}++;
			print SHARED $indel_key . "\t" . $varscan_indels{$indel_key} . "\t" . $glf_indels{$indel_key} . "\t" . $gatk_indels{$indel_key} . "\n";
			print OUTFILE $indel_key . "\t" . $varscan_indels{$indel_key} . "\t" . $glf_indels{$indel_key} . "\t" . $gatk_indels{$indel_key} . "\n";
		}
		if($indel_keys{$indel_key} eq "sniper-varscan-shared")
		{
			$stats{'sniper-varscan-shared'}++;
			print SHARED $indel_key . "\t" . $varscan_indels{$indel_key} . "\t" . $glf_indels{$indel_key} . "\n";
			print OUTFILE $indel_key . "\t" . $varscan_indels{$indel_key} . "\t" . $glf_indels{$indel_key} . "\n";
		}
		if($indel_keys{$indel_key} eq "sniper-gatk-shared")
		{
			$stats{'sniper-gatk-shared'}++;
			print SHARED $indel_key . "\t" . $gatk_indels{$indel_key} . "\t" . $glf_indels{$indel_key} . "\n";
			print OUTFILE $indel_key . "\t" . $gatk_indels{$indel_key} . "\t" . $glf_indels{$indel_key} . "\n";
		}
		if($indel_keys{$indel_key} eq "varscan-gatk-shared")
		{
			$stats{'varscan-gatk-shared'}++;
			print SHARED $indel_key . "\t" . $varscan_indels{$indel_key} . "\t" . $gatk_indels{$indel_key} . "\n";
			print OUTFILE $indel_key . "\t" . $varscan_indels{$indel_key} . "\t" . $gatk_indels{$indel_key} . "\n";
		}
		elsif($indel_keys{$indel_key} eq "varscan-only")
		{
			$stats{'varscan-only'}++;
			print VARSCAN $indel_key . "\t" . $varscan_indels{$indel_key} . "\n";
			print OUTFILE $indel_key . "\t" . $varscan_indels{$indel_key} . "\n";
		}
		elsif($indel_keys{$indel_key} eq "sniper-only")
		{
			$stats{'sniper-only'}++;
			print SNIPER $indel_key . "\t" . $glf_indels{$indel_key} . "\n";
			print OUTFILE $indel_key . "\t" . $glf_indels{$indel_key} . "\n";
		}
		elsif($indel_keys{$indel_key} eq "gatk-only")
		{
			$stats{'gatk-only'}++;
			print GATK $indel_key . "\t" . $gatk_indels{$indel_key} . "\n";
			print OUTFILE $indel_key . "\t" . $gatk_indels{$indel_key} . "\n";
		}
		else
		{
			## Where did this come from !? ##
		}
	}
	
	close(OUTFILE);
	close(SHARED);
	close(SNIPER);
	close(VARSCAN);
	close(GATK);
	
	## Sort the file ##
	
	my $cmd_obj = Genome::Model::Tools::Capture::SortByChrPos->create(
	    input_file => $output_file,
	    output_file => $output_file,
	);

	$cmd_obj->execute();
	
	## Reset stats if necessary ##
	$stats{'tri-shared'} = 0 if(!$stats{'tri-shared'});
	$stats{'sniper-varscan-shared'} = 0 if(!$stats{'sniper-varscan-shared'});
	$stats{'sniper-gatk-shared'} = 0 if(!$stats{'sniper-gatk-shared'});
	$stats{'varscan-gatk-shared'} = 0 if(!$stats{'varscan-gatk-shared'});
	$stats{'shared'} = 0 if(!$stats{'shared'});
	$stats{'varscan-only'} = 0 if(!$stats{'varscan-only'});
	$stats{'sniper-only'} = 0 if(!$stats{'sniper-only'});
	$stats{'gatk-only'} = 0 if(!$stats{'gatk-only'});

	print $stats{'total'} . " unique indels\n";
	if ($self->gatk_file) {
		print $stats{'tri-shared'} . " shared in gatk, sniper, and varscan\n";
		print $stats{'sniper-varscan-shared'} . " shared in sniper and varscan\n";
		print $stats{'sniper-gatk-shared'} . " shared in gatk and sniper\n";
		print $stats{'varscan-gatk-shared'} . " shared in gatk and varscan\n";
	}
	else {
		print $stats{'shared'} . " shared\n";		
	}
	print $stats{'varscan-only'} . " Varscan-only\n";
	print $stats{'sniper-only'} . " Sniper-only\n";
	if ($self->gatk_file) {
		print $stats{'gatk-only'} . " gatk-only\n";
	}

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_indels
{
	my $infile = shift(@_);
	
	my %indels = ();

	my $input = new FileHandle ($infile);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chrom, my $chr_start, my $chr_stop, my $allele1, my $allele2) = split(/\t/, $line);
		
		if(lc($chrom) ne "chrom")	## Skip header lines ##
		{
			$allele1 = "-" if($allele1 eq "0");
			$allele2 = "-" if($allele2 eq "0");
			
			my $indel_key = "$chrom\t$chr_start\t$chr_stop\t$allele1\t$allele2";
			
			my @lineContents = split(/\t/, $line);
			my $numContents = @lineContents;
			
			## Build rest of line ##
			my $rest_of_line = "";
			for(my $colCounter = 5; $colCounter < $numContents; $colCounter++)
			{
				$rest_of_line .= "\t" if($rest_of_line);
				$rest_of_line .= $lineContents[$colCounter];
			}
			
			## Save indel ##
			
			$indels{$indel_key} = $rest_of_line;
		}
	}
	
	close($input);
	
	return(%indels);
}


1;

