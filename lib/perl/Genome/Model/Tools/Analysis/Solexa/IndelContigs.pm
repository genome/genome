
package Genome::Model::Tools::Analysis::Solexa::IndelContigs;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# IndelContigs - Given a set of SNP variants, build reference and variant contigs for alignment purposes
#			Note: if there are several variants near one another, many near-matching contigs will be generated, and
#			Bowtie will not place reads that look like repeats.  Thus, it may be best to run this analysis in batches.
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
use Bio::DB::Fasta;

my $ref_dir = "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/";

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Solexa::IndelContigs {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File of indels in chrom, pos, type, size, context TSV format" },
		output_dir	=> { is => 'Text', doc => "Output dir for FASTA and variant files [contig_dir]", is_optional => 1 },
		flank_size	=> { is => 'Text', doc => "Number of flanking bases to include [48]", is_optional => 1},
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
	my $output_dir = "contig_dir";
	$output_dir = $self->output_dir if($self->output_dir);

	mkdir($output_dir) if(!(-d $output_dir));

	my $flank_size = 48;
	$flank_size = $self->flank_size if($self->flank_size);

	if(!(-e $variants_file))
	{
		die "Error: Variants file not found!\n";
	}

	## Build reference db ##
	my $reference_db;

	if($ref_dir && -d $ref_dir)
	{
	    ## Build the Bio::DB::Fasta index ##
	     $reference_db = Bio::DB::Fasta->new($ref_dir) or die "!$\n";
	}	

	## Variants by chrom ##
	open(ALLCONTIGS, ">$output_dir/contigs.fasta");
	
	my %variants_by_chrom = ();
	
 	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter >= 1)
		{
			(my $chromosome, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $indel_context) = split(/\t/, $line);
			
			if($chromosome && $chromosome ne "chrom" && $chromosome ne "chromosome" && $chromosome ne "ref_name")
			{
				print "Saving $chromosome $chr_start\t$chr_stop\t$indel_type\t$indel_size\n";

				my $indel_name = $chromosome . '_' . $chr_start . '_' . $chr_stop . '_' . $indel_type . '_' . $indel_size;

				my $flank_5 = uc($reference_db->seq($chromosome, ($chr_start - $flank_size) => ($chr_start - 1)));
				my $flank_3 = uc($reference_db->seq($chromosome, ($chr_stop + 1) => ($chr_stop + 1 + $flank_size)));

				my @tempArray = split(/[\[\/\]]/, $indel_context);
				my $allele1 = $tempArray[1];
				my $allele2 = $tempArray[2];
				$allele1 =~ s/[^ACGTN]//g;
				$allele2 =~ s/[^ACGTN]//g;
				
				my $ref_contig = $flank_5 . $allele1 . $flank_3;
				my $var_contig = $flank_5 . $allele2 . $flank_3;				

				## Adjust for indel content ##
				
				my $desired_flank_size = 2 * $flank_size + 1;
				$ref_contig = substr($ref_contig, 0, $desired_flank_size);
				$var_contig = substr($var_contig, 0, $desired_flank_size);
				
				print ALLCONTIGS ">$indel_name" . "_ref\n";
				print ALLCONTIGS "$ref_contig\n";

				print ALLCONTIGS ">$indel_name" . "_var\n";
				print ALLCONTIGS "$var_contig\n";

#				print "$indel_context\n$flank_5" . '[' . "$allele1/$allele2]$flank_3\n";

#				$variants_by_chrom{$chromosome} .= "$line\n";
				$variants_by_chrom{$chromosome . "\t" . $chr_start} = "$indel_type\t$indel_size\t$indel_context\n";
			}
			else
			{

			}
		}
	}
	
	close($input);
	print "$lineCounter variants\n";
	close(ALLCONTIGS);

	exit(0);

	my $contig_number = 0;
	my %remaining_variants = ();
	my %included_variants = ();
	my $num_included = 0;
	my $num_remaining = $lineCounter;

	while($num_remaining)
	{
		$contig_number++;
		
		open(OUTFILE, ">$output_dir/contigs-$contig_number.indels");
		open(OUTFASTA, ">$output_dir/contigs-$contig_number.fasta");
		
		$num_included = $num_remaining = 0;
		
		foreach my $chrom_key (keys %variants_by_chrom)
		{
			(my $chromosome, my $position) = split(/\t/, $chrom_key);
	
			my $check_flag = 0;
			for(my $check = $position - $flank_size; $check <= $position + $flank_size; $check++)
			{
				if($included_variants{$chromosome . "\t" . $check})
				{
					$check_flag = 1;
				}
			}
			
			if(!$check_flag)
			
			{
				(my $allele1, my $allele2) = split(/\t/, $variants_by_chrom{$chrom_key});
				$allele1 =~ s/[^ACGTN]//;
				$allele2 =~ s/[^ACGTN]//;
				my $contig_name = $chromosome . "_" . $position . "_" . $allele1 . "_" . $allele2;
				my $flank_5 = uc($reference_db->seq($chromosome, ($position - 1 - $flank_size) => ($position - 1)));
				my $flank_3 = uc($reference_db->seq($chromosome, ($position + 1) => ($position + 1 + $flank_size)));
				my $ref_contig = $flank_5 . $allele1 . $flank_3;
				my $var_contig = $flank_5 . $allele2 . $flank_3;
				print OUTFILE "$chrom_key\t$allele1\t$allele2\n";

				print OUTFASTA ">" . $contig_name . "_ref\n";
				print OUTFASTA $ref_contig . "\n";
				print OUTFASTA ">" . $contig_name . "_var\n";
				print OUTFASTA $var_contig . "\n";

				
				
				$num_included++;
				$included_variants{$chrom_key} = 1;
			}
			else
			{
				$remaining_variants{$chrom_key} = $variants_by_chrom{$chrom_key};
				$num_remaining++;
			}
		}
		
		close(OUTFASTA);
		close(OUTFILE);
		%variants_by_chrom = ();
		%variants_by_chrom = %remaining_variants;
		%remaining_variants = %included_variants = ();
		## Calc percent remaining ##
#		my $pct_remain = $num_remaining / ($num_included + $num_remaining) * 100;
#		$pct_remain = sprintf("%.2f", $pct_remain) . '%';
#		print "$num_included included, $num_remaining remain ($pct_remain)\n";

		$num_remaining = 0;	# Hack
	}

	print "$contig_number contig files created\n";
#	exit(0);

	foreach my $chrom (sort keys %variants_by_chrom)
	{
		my @chrom_variants = split(/\n/, $variants_by_chrom{$chrom});
		my $num_chrom_variants = @chrom_variants;

		my %variant_positions = ();
		my $num_included = 0;

		foreach my $variant (@chrom_variants)
		{
			(my $chromosome, my $position, my $allele1, my $allele2) = split(/\t/, $variant);

			## Check to see if there's a variant already registered within 50 bp ##
			my $check_flag = 0;
			for(my $check = $position - $flank_size; $check <= $position + $flank_size; $check++)
			{
				if($variant_positions{$check})
				{
					$check_flag = 1;
				}
			}
			
			if(!$check_flag)
			{
				$num_included++;
				$variant_positions{$position} = 1;
			}
			
		}
		
		print "$chrom\t$num_chrom_variants\t$num_included\n";
	}

	## Open output file ##
	
#	open(OUTFASTA, ">$output_fasta") or die "Can't open outfile: $!\n";

 
 #	close(OUTFASTA);

 
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

