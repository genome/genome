
package Genome::Model::Tools::Analysis::Indels::CompileBowtieResults;     # rename this when you give the module file a different name <--

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

my %stats = ();

class Genome::Model::Tools::Analysis::Indels::CompileBowtieResults {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file	=> { is => 'Text', doc => "Indels in annotation format" },
		alignment_files	=> { is => 'Text', doc => "Alignments in Bowtie native format [comma-separated]" },
		output_file	=> { is => 'Text', doc => "Output of indels with alignment counts", is_optional => 1 },
		max_q2_run	=> { is => 'Text', doc => "Maximum number of consecutive Q=2 bases before read is discarded", is_optional => 1, default=>30 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compile results of Bowtie alignments to indel contigs"                 
}

sub help_synopsis {
    return <<EOS
This command compiles results of Bowtie alignments to indel contigs
EXAMPLE:	gmt analysis indels compile-bowtie-results --variant-file [indels.formatted.tsv] --alignment-files s_1_sequence.contigs.bowtie --output-file [results.bowtie.tsv]
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
	my $variant_file = $self->variant_file;
	my $alignment_files = $self->alignment_files;
	my $output_file = $self->output_file;
	my $max_q2_run = $self->max_q2_run;


	$stats{'num_indels'} = $stats{'covered_0x'} = $stats{'covered_1x'} = $stats{'covered_10x'} =  $stats{'covered_50x'} = $stats{'covered_100x'} = 0;
	$stats{'num_alignments'} = $stats{'num_alignments_q2_fail'} = $stats{'num_good_alignments'} = 0;

	if($output_file)
	{
		open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	}

	## Get the alignment files ##
	
	my %alignment_totals = ();
	
	my @alignment_files = split(/\,/, $alignment_files);
	
	foreach my $alignment_file (@alignment_files)
	{
		print "Parsing $alignment_file\n";

		## Parse file for alignment counts ##
		my %alignment_counts = parse_alignments($self, $alignment_file);

		## Add these counts to total ##
		
		foreach my $contig (keys %alignment_counts)
		{
			$alignment_totals{$contig} = 0 if(!$alignment_totals{$contig});
			$alignment_totals{$contig} += $alignment_counts{$contig};
		}
	}

		
	## Parse the variants file ##
	
	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		(my $chrom, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
		$chrom =~ s/[^0-9XYMNT\_]//g;
		$chr_start =~ s/[^0-9]//g if($chr_start);
		$chr_stop =~ s/[^0-9]//g if($chr_stop);
	
		if($chrom && $chr_start && $chr_stop)
		{
			$stats{'num_indels'}++;
			my $indel_name = my $indel_type = my $indel_size = my $allele = "";
			
			if($ref eq "0" || $ref eq "-")
			{
				$indel_type = "Ins";
				$indel_size = length($var);
                $allele = uc($var);

				## Build indel name ##			
				$indel_name = "$chrom:$chr_start-$chr_stop:$indel_type:$indel_size:$allele";
			}
			else
			{
				$indel_type = "Del";
				$indel_size = length($ref);
                $allele = uc($ref);

				## Build indel name ##			
				$indel_name = "$chrom:$chr_start-$chr_stop:$indel_type:$indel_size:$allele";
			}
			
			my $ref_contig_name = $indel_name . "_ref";
			my $var_contig_name = $indel_name . "_var"; 
			
			my $coverage = my $reads1 = my $reads2 = my $var_freq = 0;
			
			$reads1 = $alignment_totals{$ref_contig_name} if($alignment_totals{$ref_contig_name});
			$reads2 = $alignment_totals{$var_contig_name} if($alignment_totals{$var_contig_name});
			$coverage = $reads1 + $reads2;

			if($coverage)
			{
				$var_freq = $reads2 / $coverage * 100;
				$var_freq = sprintf("%.2f", $var_freq) . '%';
				
				$stats{'covered_1x'}++;
				$stats{'covered_10x'}++ if($coverage >= 10);
				$stats{'covered_50x'}++ if($coverage >= 50);
				$stats{'covered_100x'}++ if($coverage >= 100);
			}
			else
			{
				$stats{'covered_0x'}++;
			}
			
			if($output_file)
			{
				print OUTFILE join("\t", $line, $coverage, $reads1, $reads2, $var_freq) . "\n";
			}

			print "$chrom\t$chr_start\t$chr_stop\t$ref\t$var\t$indel_type-$indel_size\t";
			print "$coverage\t$reads1\t$reads2\t$var_freq\n";
			
		}
	}
	
	close($input);

	close(OUTFILE) if($output_file);
	
	print $stats{'num_alignments'} . " alignments parsed\n";
	print $stats{'num_alignments_q2_fail'} . " had >= $max_q2_run Q2 bases and were discarded\n";
	print $stats{'num_good_alignments'} . " good alignments were used\n";
	
	print $stats{'num_indels'} . " indels in file\n";	
	print $stats{'covered_0x'} . " had no coverage\n";
	print $stats{'covered_1x'} . " had at least 1x coverage\n";
	print $stats{'covered_10x'} . " had at least 10x coverage\n";
	print $stats{'covered_50x'} . " had at least 50x coverage\n";
	print $stats{'covered_100x'} . " had at least 100x coverage\n";

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Parse alignments - parse the bowtie alignments 
#
################################################################################################

sub parse_alignments
{
	my $self = shift(@_);
	my $FileName = shift(@_);
	my %alignment_counts = ();
	
	my $max_q2_run = $self->max_q2_run;	
	
	my $q2_string = '#' x $max_q2_run;
	
	## Parse the variants file ##
	
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);
		my $contig_name = $lineContents[2];
		my $read_quals = $lineContents[5];
		
		$stats{'num_alignments'}++;
		
		if($read_quals =~ $q2_string)
		{
			$stats{'num_alignments_q2_fail'}++;
		}
		else
		{
			$alignment_counts{$contig_name}++;
			$stats{'num_good_alignments'}++;
		}
	}
	
	close($input);
	
	return(%alignment_counts);
}


1;

