
package Genome::Model::Tools::Analysis::Solexa::ContigResults;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ContigResults - Get the variant contig results
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/20/2009 by D.K.
#	MODIFIED:	04/20/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

my $ref_dir = "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/";

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Solexa::ContigResults {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File of SNPs in chrom, pos, ref, var TSV format" },
		alignments_file => { is => 'Text', doc => "Alignments (Bowtie) to variant contigs" },
		output_file	=> { is => 'Text', doc => "Optional output file", is_optional => 1},
	],
};

#, is_optional => 1

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Build reference and variant contigs for alignment purposes"                 
}

sub help_synopsis {
    return <<EOS
This command builds reference and variant contigs for alignment purposes
EXAMPLE:	gmt analysis variant-contigs --variants-file test.snps
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
	my $variants_file = $self->variants_file;
	my $alignments_file = $self->alignments_file;
	my $output_file = $self->output_file;

	if(!(-e $variants_file))
	{
		die "Error: Variants file not found!\n";
	}

	if(!(-e $alignments_file))
	{
		die "Error: Alignments file not found!\n";
	}

	my %stats = ();
	my %read_counts = ();

	print "Parsing the alignments file...\n";

	## Parse the alignments file ##

 	my $input = new FileHandle ($alignments_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my @lineContents = split(/\t/, $line);
		my $align_strand = $lineContents[1];
		my $contig_name = $lineContents[2];
		if($contig_name)
		{		
			$read_counts{$contig_name}++;
		}
	}
	
	close($input);

#	foreach my $contig (keys %read_counts)
#	{
#		print "$read_counts{$contig}\t$contig\n";
#	}


	if($output_file)
	{
		open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
		print OUTFILE "chrom\tposition\tref\tvar\tref_reads\tvar_reads\n";
	}

	print "Parsing the variants file...\n";

	## Parse the variants file ##

 	$input = new FileHandle ($variants_file);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter > 1)
		{
			my @lineContents = split(/\t/, $line);
			
			if($lineContents[0] && $lineContents[0] ne "chromosome" && $lineContents[0] ne "chrom" && $lineContents[0] ne "ref_name")
			{
				if($lineContents[3] && ($lineContents[3] eq "INSERTION" || $lineContents[3] eq "DELETION"))
				{
					## HANDLE INDELS ##
					
					(my $chromosome, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $indel_context) = split(/\t/, $line);

					my $indel_name = $chromosome . '_' . $chr_start . '_' . $chr_stop . '_' . $indel_type . '_' . $indel_size;
					my $contig_name_ref = $indel_name . "_ref";
					my $contig_name_var = $indel_name . "_var";
	
					## Get the read counts for each contig ##
					
					my $reads_ref = $read_counts{$contig_name_ref};
					my $reads_var = $read_counts{$contig_name_var};
					$reads_ref = 0 if(!$reads_ref);
					$reads_var = 0 if(!$reads_var);
					
					$stats{'num_indels'}++;
					
					if(($reads_ref + $reads_var) >= 10)
					{
						$stats{'num_indels_cov'}++;

						$stats{'num_supported_1read'}++ if($reads_var >= 1);
						$stats{'num_supported_2reads'}++ if($reads_var >= 2);
					}
					

					
					print "$chromosome\t$chr_start\t$chr_stop\t$indel_type\t$indel_size\t$indel_context\t$reads_ref\t$reads_var\n";
					print OUTFILE "$chromosome\t$chr_start\t$chr_stop\t$indel_type\t$indel_size\t$indel_context\t$reads_ref\t$reads_var\n" if($output_file);
				}
				else
				{
					## HANDLE SNPS ##
					
					(my $chromosome, my $position, my $allele1, my $allele2) = split(/\t/, $line);
					my $contig_name_ref = $chromosome . "_" . $position . "_" . $allele1 . "_" . $allele2 . "_ref";
					my $contig_name_var = $chromosome . "_" . $position . "_" . $allele1 . "_" . $allele2 . "_var";
	
					## Get the read counts for each contig ##
					
					my $reads_ref = $read_counts{$contig_name_ref};
					my $reads_var = $read_counts{$contig_name_var};
					$reads_ref = 0 if(!$reads_ref);
					$reads_var = 0 if(!$reads_var);
					
					print "$chromosome\t$position\t$allele1\t$allele2\t$reads_ref\t$reads_var\n";
					print OUTFILE "$chromosome\t$position\t$allele1\t$allele2\t$reads_ref\t$reads_var\n" if($output_file);
				}
			}
		}
	}
	
	close($input);

	print "$stats{'num_indels'} indels attempted\n";
	print "$stats{'num_indels_cov'} had >10x solexa coverage\n";
	print "$stats{'num_supported_1read'} supported by 1+ solexa reads\n";
	print "$stats{'num_supported_2reads'} supported by 2+ solexa reads\n";
 
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

