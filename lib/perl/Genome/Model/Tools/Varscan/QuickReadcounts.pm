
package Genome::Model::Tools::Varscan::QuickReadcounts;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Varscan::Somatic	Runs Varscan somatic pipeline on Normal/Tumor BAM files
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/29/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Varscan::QuickReadcounts {
	is => 'Genome::Model::Tools::Varscan',
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		bam_file	=> { is => 'Text', doc => "Path to BAM file", is_optional => 0 },
		variants_file	=> { is => 'Text', doc => "Path to variant positions file", is_optional => 0 },
		output_file	=> { is => 'Text', doc => "Path to output file" , is_optional => 0},
		reference        => { is => 'Text', doc => "Reference FASTA file for BAMs" , is_optional => 1, default_value => (Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa')},
		min_coverage	=> { is => 'Text', doc => "Minimum base coverage to report readcounts [8]" , is_optional => 1},
		min_base_qual	=> { is => 'Text', doc => "Minimum base quality to count a read [30]" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Uses BAM index to retrieve pileup and generates read counts"                 
}

sub help_synopsis {
    return <<EOS
Runs Varscan readcounts from indexed BAM output
EXAMPLE:	gmt varscan readcounts --bam-file [sample.bam] --variants-file [variants.tsv] --output-file readcounts.txt ...
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
	my $bam_file = $self->bam_file;
	my $reference = $self->reference;
	my $variants_file = $self->variants_file;
	my $output_file = $self->output_file;
	my $min_coverage = 4;
	$min_coverage = $self->min_coverage if($self->min_coverage);
	

	if(-e $bam_file)
	{
		if(-e $variants_file)
		{
			## parse the variants, generating the pileup files ##

			open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
			print OUTFILE "chrom\tposition\treads1\treads2\tvar_freq\n";

			if(-e "$output_file.pileup")
			{
				print "Removing existing file $output_file.pileup\n";
				system("rm -rf $output_file.pileup");
			}
			

			my %target_positions = ();

			my $query_string = "";

			my $input = new FileHandle ($variants_file);
			my $lineCounter = 0;
			
			while (<$input>)
			{
				chomp;
				my $line = $_;
				$lineCounter++;
				
				(my $chromosome, my $position) = split(/\t/, $line);
				
				$chromosome =~ s/[^0-9XYMT]//g;
				$position =~ s/[^0-9]//g;
				
				if($chromosome && $position)# && $lineCounter < 10)
				{
					$query_string .= " " if($query_string);
					$query_string .= $chromosome . ":" . $position . "-" . $position;
					## Run Command ##
				#	my $cmd = "samtools view -b -u -q 1 $bam_file $chromosome:$position-$position | samtools pileup -f $reference - >>$output_file.pileup";# - | java -classpath ~dkoboldt/Software/Varscan net.sf.varscan.Varscan pileup2snp --min-coverage $min_coverage";
				#	print "RUN $cmd\n";
				#	system("$cmd");
					
				}
			}
			
			close($input);

			## Run Pileup Command ##
			my $cmd = "samtools view -b -u -q 1 $bam_file $query_string | samtools mpileup -f $reference - >>$output_file.pileup";# - | $self->java_command_line("pileup2snp --min-coverage $min_coverage");
			print "RUN $cmd\n";
			system("$cmd");

			## Run Limit Command ##
			$cmd = $self->java_command_line("limit $output_file.pileup --positions-file $variants_file --output-file $output_file.pileup.roi");# - | pileup2snp --min-coverage $min_coverage";
			print "RUN $cmd\n";
			system("$cmd");

			## Run Varscan
			
			$cmd = $self->java_command_line("pileup2cns $output_file.pileup.roi --min-coverage $min_coverage >$output_file.pileup.roi.cns");
			print "RUN $cmd\n";
			system("$cmd");

			## Load the results ##
			
			print "Loading results...\n";
			my %varscan_results = load_results("$output_file.pileup.roi.cns");


			## Parse the input file again ##

			
			$input = new FileHandle ($variants_file);
			$lineCounter = 0;
			
			
			while (<$input>)
			{
				chomp;
				my $line = $_;
				$lineCounter++;
				
				(my $chromosome, my $position) = split(/\t/, $line);
				
				$chromosome =~ s/[^0-9XYMT]//g;
				$position =~ s/[^0-9]//g;
				
				if($chromosome && $position)# && $lineCounter < 10)
				{
					if($varscan_results{"$chromosome\t$position"})
					{
						(my $ref, my $cns, my $reads1, my $reads2, my $var_freq) = split(/\t/, $varscan_results{"$chromosome\t$position"});
						
						print OUTFILE "$chromosome\t$position\t$reads1\t$reads2\t$var_freq\n";
						print "$chromosome\t$position\t$reads1\t$reads2\t$var_freq\n";
					}
					else
					{
						print OUTFILE "$chromosome\t$position\t0\t0\t0\n";
						print "$chromosome\t$position\t0\t0\t0\n";
					}
				}
			}

		}
		
		## Prepare pileup commands ##
		
#		my $pileup = "samtools pileup -f $reference $bam_file";
#		my $cmd = "bash -c \"java -classpath ~dkoboldt/Software/Varscan net.sf.varscan.Varscan readcounts <\($pileup\) --variants-file $variants_file --output-file $output_file\"";
#		system($cmd);
	}
	else
	{
		die "Error: One of your BAM files doesn't exist!\n";
	}
	
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}






################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_results
{
	my $FileName = shift(@_);	

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	my %results = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chromosome, my $position, my $ref, my $cns, my $reads1, my $reads2, my $var_freq) = split(/\t/, $line);
		$results{"$chromosome\t$position"} = "$ref\t$cns\t$reads1\t$reads2\t$var_freq";
	}
	close($input);

	return(%results);
}


1;

