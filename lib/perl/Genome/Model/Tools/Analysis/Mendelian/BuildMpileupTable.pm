
package Genome::Model::Tools::Analysis::Mendelian::BuildMpileupTable;     # rename this when you give the module file a different name <--

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

my $num_affected = my $affecteds_missing = my $unaffecteds_variant = my $affecteds_variant = my $affecteds_ambiguous = 0;
my %genotypes = ();

class Genome::Model::Tools::Analysis::Mendelian::BuildMpileupTable {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sample_file	=> { is => 'Text', doc => "Tab-delimited file of family, sample, status, dir", is_optional => 0, is_input => 1},
		variant_file	=> { is => 'Text', doc => "Variants in sorted annotation format to build mpileup for", is_optional => 0, is_input => 1},
		family	=> { is => 'Text', doc => "Family for which to build mpileup", is_optional => 0, is_input => 1},
		reference	=> { is => 'Text', doc => "Path to the reference to use [defaults to build 37]", is_optional => 0, is_input => 1, default => '/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa'},
		output_file	=> { is => 'Text', doc => "Output directory to contain files", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Runs a Mendelian analysis of variants on a per-family basis"                 
}

sub help_synopsis {
    return <<EOS
This command runs a Mendelian analysis of variants on a per-family basis
EXAMPLE:	gmt analysis mendelian germline-comparison --sample-file Family-Sample-Status-Dir.tsv --output-dir mendelian_out
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

	my $sample_file = $self->sample_file;
	my $output_file = $self->output_file;
	my $variant_file = $self->variant_file;
	my $reference = $self->reference;
	
	my %stats = ();

	my $sample_list = my $bam_list = "";

	## Print the variants ##

	my $input = new FileHandle ($sample_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($family, $sample_name, $affected_status, $dir) = split(/\t/, $line);	

		if($family eq $self->family && $affected_status && $affected_status ne "control")
		{
			my $bam_file = `ls $dir/alignments/*.bam | head -1`;
			chomp($bam_file);

			$sample_list .= "\t" if($sample_list);
			$sample_list .= $sample_name;
			
			$bam_list .= " " if($bam_list);
			$bam_list .= $bam_file;
		}
	}
	
	close($input);
	
	print "$sample_list\n$bam_list\n";

	$input = new FileHandle ($variant_file);
	$lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);		

		## Determine variant type ##
		my $variant_type = "SNP";
		
		if($ref eq '-' || $ref eq '0' || length($var) > 1)
		{
			$variant_type = "INS";
		}
		elsif($var eq '-' || $var eq '0' || length($ref) > 1)
		{
			$variant_type = "DEL";
		}		

		if($variant_type eq "INS")
		{
			my $indel_size = length($var);	
			my $query = $chrom . ":" . ($chr_start - $indel_size - 1) . "-" . ($chr_stop + $indel_size + 1);
			my $cmd = "samtools mpileup -f $reference -q 10 -r $query $bam_list 2>/dev/null >>$output_file";
			print join("\t", $lineCounter, $line) . "\n";
			system($cmd);
		}
		elsif($variant_type eq "DEL")
		{
			my $indel_size = length($ref);	
			my $query = $chrom . ":" . ($chr_start - $indel_size - 1) . "-" . ($chr_stop + $indel_size + 1);
			my $cmd = "samtools mpileup -f $reference -q 10 -r $query $bam_list 2>/dev/null >>$output_file";
			print join("\t", $lineCounter, $line) . "\n";
			system($cmd);			
		}
		else
		{
			my $query = $chrom . ":" . $chr_start . "-" . $chr_stop;
			my $cmd = "samtools mpileup -f $reference -q 10 -r $query $bam_list 2>/dev/null >>$output_file";
			print join("\t", $lineCounter, $line) . "\n";
			system($cmd);			
		}
		
#		return(0) if($lineCounter > 10);
	}
	
	close($input);

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# sorting subroutine by chromosome and then position
#
################################################################################################

sub byChrPos
{
	my ($chrom_a, $pos_a) = split(/\t/, $a);
	my ($chrom_b, $pos_b) = split(/\t/, $b);

	$chrom_a = 23 if($chrom_a eq 'X');
	$chrom_a = 24 if($chrom_a eq 'Y');
	$chrom_a = 25 if($chrom_a eq 'MT');
	
	$chrom_b = 23 if($chrom_b eq 'X');
	$chrom_b = 24 if($chrom_b eq 'Y');
	$chrom_b = 25 if($chrom_b eq 'MT');

	$chrom_a =~ s/[^0-9]//g;
	$chrom_b =~ s/[^0-9]//g;
	
	$chrom_a <=> $chrom_b
	or
	$pos_a <=> $pos_b;
}


1;


