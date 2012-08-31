
package Genome::Model::Tools::Analysis::SomaticPipeline::MergeIndelsWithAnnotation;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MergeIndelsWithAnnotation - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	10/23/2009 by D.K.
#	MODIFIED:	10/23/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my %stats = ();

class Genome::Model::Tools::Analysis::SomaticPipeline::MergeIndelsWithAnnotation {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File of variants in indel format", is_optional => 0 },
		annotation_file	=> { is => 'Text', doc => "Annotate-indel output file", is_optional => 0 },
		output_file     => { is => 'Text', doc => "Output file to receive merged data", is_optional => 0 },
		full_annotation => { is => 'Text', doc => "If set to 1, appends the full annotation to line", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merges indels with their annotations"                 
}

sub help_synopsis {
    return <<EOS
This command merges variant calls from the pipeline with their annotation information
EXAMPLE:	gmt analysis somatic-pipeline merge-snvs-with-annotation --variants-file [file] --annotation-file [file] --output-file [file]
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
	my $annotation_file = $self->annotation_file;
	my $output_file = $self->output_file;

	## Determine if variant file is in bed format ##
	my $is_bed = 0;

	if(substr($variants_file, length($variants_file) - 3, 3) eq "bed")
	{
		$is_bed = 1;
	}
	
	## Load the SNP calls ##

	my %snp_calls = load_snp_calls($variants_file);
	
	## Open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	open(TIER1, ">$output_file.tier1") or die "Can't open tier1 outfile: $!\n";
	
	## Parse the annotation file ##
	
	my $input = new FileHandle ($annotation_file);
	my $lineCounter = 0;
	
	$stats{'tier1'} = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);

		if(!($lineContents[0] =~ "chrom" || $lineContents[0] =~ "ref_name"))
		{
			my $chrom = $lineContents[0];
			$chrom =~ s/[^0-9XYMTNT\_]//g;
			my $chr_start = $lineContents[1];
			my $chr_stop = $lineContents[2];
			my $allele1 = $lineContents[3];
			my $allele2 = $lineContents[4];
			my $variant_type = $lineContents[5];
			my $gene_name = $lineContents[6];
			my $transcript_name = $lineContents[7];
			my $trv_type = $lineContents[13];
			my $c_position = $lineContents[14];
			my $aa_change = $lineContents[15];
			my $ucsc_cons = $lineContents[16];
			my $domain = $lineContents[17];

			if($allele1 eq "-" || $allele1 eq "0")
			{
				$allele1 = "-";
			}
			else
			{
				$allele2 = "-";
			}
			
			my $snp_key = "$chrom\t$chr_start\t$chr_stop\t$allele1\t$allele2" if($chrom && $chr_start && $chr_stop);
			
			if($variant_type ne "SNP" && $snp_key && $snp_calls{$snp_key})
			{
				my $annotation = "$variant_type\t$gene_name\t$transcript_name\t$trv_type\t$c_position\t$aa_change\t$ucsc_cons\t$domain";
				
				if($self->full_annotation)
				{
					$annotation = $line;
				}

#				print OUTFILE "$annotation\t$snp_calls{$snp_key}\n";
				my $newline = $snp_key . "\t";
				
				if($is_bed)
				{
					## For a bed file, do slightly alternate ##
					$newline = "$chrom\t$chr_start\t$chr_stop\t$allele1/$allele2\t";
				}
				
				$newline .= $annotation . "\t";
				$newline .= $snp_calls{$snp_key};

				print OUTFILE "$newline\n";

				#my $newline = "$chrom\t$position\t$position\t$allele1\t$allele2\tSNP\t$annotation\t$iupac_code\t$p_value";
				
				if($trv_type eq "missense" || $trv_type eq "nonsense" || $trv_type eq "nonstop" || $trv_type =~ "splice_site" || $trv_type =~ "frame" || $trv_type =~ "silent" || $trv_type =~ "rna")
				{
					print TIER1 "$newline\n";
#					print "$gene_name\t$trv_type\t$c_position $aa_change\t$ucsc_cons\t$snp_key\n";
					$stats{'tier1'}++;
					$stats{$trv_type}++;
				}

			}
		}
	}

	close($input);

	$stats{'missense'} = 0 if(!$stats{'missense'});
	$stats{'nonsense'} = 0 if(!$stats{'nonsense'});
	$stats{'nonstop'} = 0 if(!$stats{'nonstop'});
	$stats{'splice_site'} = 0 if(!$stats{'splice_site'});
	$stats{'in_frame_ins'} = 0 if(!$stats{'in_frame_ins'});
	$stats{'in_frame_del'} = 0 if(!$stats{'in_frame_del'});
	$stats{'frame_shift_ins'} = 0 if(!$stats{'frame_shift_ins'});
	$stats{'frame_shift_del'} = 0 if(!$stats{'frame_shift_del'});
	$stats{'silent'} = 0 if(!$stats{'silent'});

	print $stats{'tier1'} . " tier 1 mutations (";
	
	print $stats{'missense'} . " missense, " if($stats{'missense'});
	print $stats{'nonsense'} . " nonsense, " if($stats{'nonsense'});
	print $stats{'nonstop'} . " nonstop, " if($stats{'nonstop'});
	print $stats{'splice_site'} . " splice_site, " if($stats{'splice_site'});
	print $stats{'in_frame_ins'} . " in_frame_ins, " if($stats{'in_frame_ins'});
	print $stats{'in_frame_del'} . " in_frame_del, " if($stats{'in_frame_del'});
	print $stats{'frame_shift_ins'} . " frame_shift_ins, " if($stats{'frame_shift_ins'});
	print $stats{'frame_shift_del'} . " frame_shift_del, " if($stats{'frame_shift_del'});
	print $stats{'silent'} . " silent, " if($stats{'silent'});
	print $stats{'rna'} . " miRNA " if($stats{'rna'});
	
	print ")\n";	
		
	close(OUTFILE);
	close(TIER1);
}



#############################################################
# ParseBlocks - takes input file and parses it
#
#############################################################

sub load_snp_calls
{
	my $FileName = shift(@_);

	my %snps = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	my $is_bed = 0;

	if(substr($FileName, length($FileName) - 3, 3) eq "bed")
	{
		$is_bed = 1;
	}
	my @formatted = ();
	my $formatCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		my $chrom = $lineContents[0];
		
		if(!($lineContents[0] =~ "chrom" || $lineContents[0] =~ "ref_name"))
		{
			$chrom =~ s/[^0-9XYMTNT\_]//g;
			my $chr_start = $lineContents[1];
			my $chr_stop = my $allele1 = my $allele2 = "";
			my $ref = my $var = "";
			my $indel_type = my $indel_size = "";

			if($is_bed)
			{
				$chr_start = $lineContents[1];
				$chr_stop = $lineContents[2];
				($ref, $var) = split(/\//, $lineContents[3]);
				$ref =~ s/\*/\-/g;
				$var =~ s/\*/\-/g;
				if($ref eq '-')
				{
					$chr_stop++;
				}
				else
				{
					$chr_start++;
				}
				#print "$line\n";
				#print join("\t", $chrom, $chr_start, $chr_stop, $allele1, $allele2) . "\n\n";
			}
			elsif($lineContents[2] =~ /[0-9]/ && $lineContents[2] ne "0")
			{
				$chr_stop = $lineContents[2];
				$ref = $lineContents[3];
				$var = $lineContents[4];
			}
			else
			{
				$ref = $lineContents[2];
				$var = $lineContents[3];
			}

			## Correct alleles ##

			if($ref eq '-' || $var eq '-' || $ref eq "0" || $var eq "0")
			{
				$allele1 = $ref;
				$allele2 = $var;
				
				if($ref eq '-' || $ref eq "0")
				{
					$allele1 = '-';
					$indel_type = "INSERTION";
					$indel_size = length($var);
				}
				else
				{
					$allele2 = '-';
					$indel_type = "DELETION";
					$indel_size = length($ref);
				}
			}
			elsif(substr($var, 0, 1) eq '+')
			{
				$allele1 = "-";
				$allele2 = uc($var);
				$allele2 =~ s/[^ACGTN]//g;
				$indel_type = "INSERTION";
				$indel_size = length($allele2);
			}
			elsif(substr($var, 0, 1) eq '-')
			{
				$allele2 = "-";
				$allele1 = uc($var);
				$allele1 =~ s/[^ACGTN]//g;
				$indel_type = "DELETION";
				$indel_size = length($allele1);
			}
			else
			{
				warn "Unable to format $line\n";
				$chrom = $chr_start = $chr_stop = $allele1 = $allele2 = "";
			}

			## If no chr stop, calculate it ##
			if(!$chr_stop)
			{
				if($indel_type eq "INSERTION" || $indel_size == 1)
				{
					$chr_stop = $chr_start + 1;
				}
				else
				{
					$chr_stop = $chr_start + $indel_size;
				}
			}

			my $snp_key = "$chrom\t$chr_start\t$chr_stop\t$allele1\t$allele2";
		
			my $snp_call = "";

			for(my $colCounter = 4; $colCounter < $numContents; $colCounter++)
			{
				$snp_call .= "\t" if($colCounter > 4);
				$snp_call .= $lineContents[$colCounter];
			}
		
#			$snps{$snp_key} = $line;
			$snps{$snp_key} = $snp_call;
		}
		else
		{
			$snps{'header'} = $line;
		}
	}

	close($input);

	#print "$lineCounter SNP calls loaded\n";
	
	return(%snps);

	return 0;
}




1;

