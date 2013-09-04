
package Genome::Model::Tools::Analysis::454::CompileVariants;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# CompileVariants - Load 454 reads from a sample-SFF tab-delimited file
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/25/2009 by D.K.
#	MODIFIED:	02/25/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it
use Genome::Model::Tools::Analysis::Helpers qw(
    code_to_genotype
);

class Genome::Model::Tools::Analysis::454::CompileVariants {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		
		dbsnp_file	=> { is => 'Text', doc => "Path to dbSNP file from sample dir", is_optional => 1},
		annotation_file	=> { is => 'Text', doc => "Path to transcript annotation file from sample dir", is_optional => 1},
		samples_file	=> { is => 'Text', doc => "Tab-delimited file of sample and SFF file(s)" },
		output_dir	=> { is => 'Text', doc => "Output directory" },
		variants_file	=> { is => 'Text', doc => "Path to variants file from sample dir" },
		output_file	=> { is => 'Text', doc => "Output file" },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compile variant calls from 454 projects"                 
}

sub help_synopsis {
    return <<EOS
This command compiles variants from 454 projects
	
EXAMPLE: gmt analysis 454 compile-variants --samples-file data/samples.tsv --output-dir data --variants merged.snp --dbsnp merged.snp.known --output MyVariantCalls.tsv
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
	my $samples_file = $self->samples_file;
	my $output_dir = $self->output_dir;
	my $variants_file = $self->variants_file;
	my $dbsnp_file = $self->dbsnp_file;
	my $annotation_file = $self->annotation_file;	
	my $output_file = $self->output_file;
	
	## Parse the samples file ##
	
	my %all_samples = ();
	my %all_variants = ();
	my %all_annotation = ();
	my %all_dbsnp = ();
	my %sample_genotypes = ();
	my %all_sample_cns = ();
	
	## Open the output file ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: !$\n";
	
	
	my $input = new FileHandle ($samples_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my @lineContents = split(/\t/, $line);			
		my $sample_name = $lineContents[0];

		my $sample_dir = $output_dir . "/" . $sample_name;
		
		my $sample_variants_file = $sample_dir . "/" . $variants_file;
		my $sample_dbsnp_file = $sample_dir . "/" . $dbsnp_file;
		my $sample_annotation_file = $sample_dir . "/" . $annotation_file;
		
		my %sample_stats = ();
		
		if(-e $sample_variants_file)
		{
			print "$sample_name\n";
			
			## Load the dbSNP variants ##
			
			my %sample_dbsnps = load_dbsnps($sample_dbsnp_file) if(-e $sample_dbsnp_file);
			
			## Load the annotation ##
			
			my %sample_annotation = load_annotation($sample_annotation_file) if(-e $sample_annotation_file);

			## Load the CNS Calls ##
			
			my $sample_cns_file = $sample_dir . "/ssaha2_out/$sample_name.ssaha2.bam.pileup.cns.allVariants";
			my %sample_cns_calls = load_cns_calls($sample_cns_file);

			my $callCounter =0;
			foreach my $position_key (keys %sample_cns_calls)
			{
				$callCounter++;
				my $cns_call = $sample_cns_calls{$position_key};
				$all_sample_cns{$position_key . "\t" . $sample_name} = $cns_call;
			}
			

			## Parse the variants file ##

			my $input2 = new FileHandle ($sample_variants_file);
			my $lineCounter2 = 0;
			
			while (<$input2>)
			{
				chomp;
				my $line2 = $_;
				$lineCounter2++;		
			
				my ($chrom, $chr_start, $chr_stop, $ref, $var, $reads1, $reads2, $var_freq) = split(/\t/, $line2);				
				my $dbsnp_key = join("\t", $chrom, $chr_start);
				my $annotation_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
				$sample_stats{'num_variants'}++;
				
				$all_samples{$sample_name} = 1;
				
				## Calculate genotype ##
				
				my $genotype = "";
				
				my $numeric_freq = $var_freq;
				$numeric_freq =~ s/\%//;
				if($numeric_freq >= 80)
				{
					$genotype = $var . $var;
				}
				else
				{
					$genotype = $ref . $var;
				}
				
				## Get dbSNP RS number ##
				
				my $rs_number = "-";
				if($sample_dbsnps{$dbsnp_key})
				{
					$rs_number = $sample_dbsnps{$dbsnp_key};
					$sample_stats{'num_dbsnps'}++;
				}
				
				## Get Annotation ##
				
				my $variant_type = my $gene = my $transcript = my $trv_type = my $codon = my $aa_change = my $conservation = my $domain1 = my $domain2 = "-";
				
				if($sample_annotation{$annotation_key})
				{
					my @annotation = split(/\t/, $sample_annotation{$annotation_key});
					$variant_type = $annotation[5];
					$gene = $annotation[6];
					$transcript = $annotation[7];
					$trv_type = $annotation[13];
					$codon = $annotation[14] if($annotation[14]);
					$aa_change = $annotation[15] if($annotation[15]);
					$conservation = $annotation[16] if($annotation[16]);
					$domain1 = $annotation[17] if($annotation[17]);
					$domain2 = $annotation[18] if($annotation[18]);
				}
				
				my $annotation_string = join("\t", $variant_type, $gene, $transcript, $trv_type, $codon, $aa_change);#, $conservation, $domain1, $domain2);


				## Save this information if not done already ##
				
				$all_variants{$annotation_key}++;
				$all_annotation{$annotation_key} = $annotation_string;
				$all_dbsnp{$annotation_key} = $rs_number;
				$sample_genotypes{$annotation_key . "\t" . $sample_name} = $genotype;

				print join("\t", $annotation_key, $reads1, $reads2, $var_freq, $rs_number, $annotation_string) . "\n";
			}
			
			close($input2);
			
			print $sample_stats{'num_variants'} . " variants\n";
			print $sample_stats{'num_dbsnps'} . " dbSNP\n";
			print "$callCounter CNS calls loaded\n";

		}
		else
		{
			warn "No variants file ($variants_file) found in $sample_dir\n";
		}
	}
	
	close($input);


	## Go through each of the samples, printing them across the top ##

	## Header ##

	print OUTFILE "chrom\tchr_start\tchr_stop\tref\tvar\trs_number\ttype\tgene\ttranscript\ttrv_type\tcodon\taa_change\tsamples_ref\tsamples_het\tsamples_hom\tallele_freq";

	foreach my $sample (sort keys %all_samples)
	{
		print OUTFILE "\t$sample";
	}
	
	print OUTFILE "\n";

	## Variants ##
	
	foreach my $key (sort byChrPos keys %all_variants)
	{
		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $key);

		## Build sample genotypes ##
		
		my $sample_gt_string = "";
		
		my $num_alleles = my $num_variant_alleles = 0;
		my $num_samples_ref = my $num_samples_het = my $num_samples_hom = 0;
		
		foreach my $sample (sort keys %all_samples)
		{
			my $genotype = "--";
			
			my $cns_key = join("\t", $chrom, $chr_start, $sample);
			
			if($all_sample_cns{$cns_key})
			{
				$genotype = $all_sample_cns{$cns_key};

				if($genotype eq ($ref . $ref))
				{
					$genotype = "..";
					$num_samples_ref++;
				}

			}
			
			$genotype = $sample_genotypes{$key . "\t" . $sample} if($sample_genotypes{$key . "\t" . $sample});
		
			if($sample_genotypes{$key . "\t" . $sample})
			{
				my ($a1, $a2) = split(//, $genotype);
				
				$num_variant_alleles++ if($a1 ne $ref);
				$num_variant_alleles++ if($a2 ne $ref);

				if($a1 ne $ref && $a2 ne $ref)
				{
					$num_samples_hom++;
				}
				elsif($a1 ne $ref || $a2 ne $ref)
				{
					$num_samples_het++;
				}
			}


			## Count these alleles if either reference or variant genotype was recorded ##

			if($genotype && $genotype ne "--")
			{
				$num_alleles += 2;				
			}

			$sample_gt_string .= "\t" if($sample_gt_string);
			$sample_gt_string .= $genotype;
		}
		

		## Calculate allele frequency ##
		
		my $allele_frequency = sprintf("%.2f", ($num_variant_alleles / $num_alleles * 100)) . '%' if($num_alleles);
		
		## Build sample genotype string ##

		print OUTFILE join("\t", $key, $all_dbsnp{$key}, $all_annotation{$key}, $num_samples_ref, $num_samples_het, $num_samples_hom, $allele_frequency, $sample_gt_string) . "\n";
	}


	close(OUTFILE);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# load_dbsnps
#
################################################################################################

sub load_dbsnps
{                               # replace with real execution logic.
	my $filename = shift(@_);
	my %dbsnps = ();

	my $input = new FileHandle ($filename);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);
		my @lineContents = split(/\t/, $line);			
		my $numContents = @lineContents;
		my $dbsnp_rs = $lineContents[$numContents - 1];
		
		my $key = join("\t", $chrom, $chr_start);
		$dbsnps{$key} = $dbsnp_rs;

	}
	
	close($input);


	return(%dbsnps);

}



################################################################################################
# load_dbsnps
#
################################################################################################

sub load_annotation
{                               # replace with real execution logic.
	my $filename = shift(@_);
	my %annotation = ();

	my $input = new FileHandle ($filename);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);		
		my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
		$annotation{$key} = $line;

	}
	
	close($input);


	return(%annotation);

}




################################################################################################
# load_dbsnps
#
################################################################################################

sub load_cns_calls
{                               # replace with real execution logic.
	my $filename = shift(@_);
	my %calls = ();

	my $input = new FileHandle ($filename);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my ($chrom, $position, $ref, $cns) = split(/\t/, $line);		
		my $key = join("\t", $chrom, $position);
		$calls{$key} = code_to_genotype($cns);

	}
	
	close($input);


	return(%calls);

}


sub byChrPos
{
    (my $chrom_a, my $pos_a) = split(/\t/, $a);
    (my $chrom_b, my $pos_b) = split(/\t/, $b);

	$chrom_a =~ s/X/23/;
	$chrom_a =~ s/Y/24/;
	$chrom_a =~ s/M/25/;

	$chrom_b =~ s/X/23/;
	$chrom_b =~ s/Y/24/;
	$chrom_b =~ s/M/25/;

    $chrom_a <=> $chrom_b
    or
    $pos_a <=> $pos_b;
    
#    $chrom_a = 23 if($chrom_a =~ 'X');
#    $chrom_a = 24 if($chrom_a =~ 'Y');
    
}


1;

