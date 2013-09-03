package Genome::Model::Tools::Germline::VcfMaker;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# GermlinePipelineFinisher - Generate Vcf File
#					
#	AUTHOR:		Will Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	04/01/2011 by W.S.
#	MODIFIED:	04/01/2011 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it
use Genome::Model::Tools::Capture::Helpers 'iupac_to_base';

class Genome::Model::Tools::Germline::VcfMaker {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => "List of mutations in Vcf format",
        },
#        maf_file => { 
#            is => 'Text',
#            doc => "List of mutations in MAF format",
#            is_input => 1,
#            is_optional => 1,
#        },
#	data_dir => {
#            is => 'Text',
#            doc => "Base Data Directory i.e. /gscmnt/sata424/info/medseq/Freimer-Boehnke/Analysis-1033Samples/ ",
#            is_optional => 0
#	},

	bam_file	=> { is => 'Text', doc => "Bam File" , is_optional => 0, is_input => 1},
	build_id	=> { is => 'Text', doc => "Build Id" , is_optional => 0, is_input => 1},
	variant_file	=> { is => 'Text', doc => "Tier 1 SNV File" , is_optional => 0, is_input => 1},
	dbsnp_file	=> { is => 'Text', doc => "dbsnp File" , is_optional => 0, is_input => 1},
	snv_filtered_file	=> { is => 'Text', doc => "Strandfilter Output File" , is_optional => 0, is_input => 1},
	snv_failfiltered_file	=> { is => 'Text', doc => "Strandfilter Failed Output File" , is_optional => 0, is_input => 1},
	snv_annotation_file	=> { is => 'Text', doc => "Annotation File" , is_optional => 0, is_input => 1},
	indel_file	=> { is => 'Text', doc => "Tier 1 Indel File" , is_optional => 0, is_input => 1},
	indel_filtered_file	=> { is => 'Text', doc => "Strandfilter Output File" , is_optional => 0, is_input => 1},
	indel_failfiltered_file	=> { is => 'Text', doc => "Strandfilter Failed Output File" , is_optional => 0, is_input => 1},
	indel_annotation_file	=> { is => 'Text', doc => "Annotation File" , is_optional => 0, is_input => 1},
	samtools_file	=> { is => 'Text', doc => "Annotation File" , is_optional => 0, is_input => 1},
	varscan_file	=> { is => 'Text', doc => "Annotation File" , is_optional => 0, is_input => 1},
	output_file	=> { is => 'Text', doc => "MAF File" , is_optional => 0, is_output => 1, is_input => 1},
	project_name	=> { is => 'Text', doc => "Name of the project i.e. ASMS" , is_optional => 1, default => "Germline Project", is_input => 1},
	center 		=> { is => 'Text', doc => "Genome center name" , is_optional => 1, default => "genome.wustl.edu", is_input => 1},
	build 		=> { is => 'Text', doc => "Reference genome build" , is_optional => 1, example_values => ["36"], is_input => 1},
	sequence_phase	=> { is => 'Text', doc => "Sequencing phase" , is_optional => 1, default => "4", is_input => 1},
	sequence_source	=> { is => 'Text', doc => "Sequence source" , is_optional => 1, default => "Capture", is_input => 1},
	sequencer	=> { is => 'Text', doc => "Sequencing platform name" , is_optional => 1, default => "Illumina_GAIIx_or_Hiseq", is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Generate Vcf File"
}

sub help_synopsis {
    return <<EOS
Generate Vcf File
EXAMPLE:	gmt germline vcf-maker
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

	my $build_id = $self->build_id;
	my $build = Genome::Model::Build->get(build_id => $build_id);
	my $sample_name = $build->subject_name;

	my $bam_file = $self->bam_file;

	my $samtools_file = $self->samtools_file;
	my $varscan_file = $self->varscan_file;

	my $variant_file = $self->variant_file;
	my $dbsnp_file = $self->dbsnp_file;
	my $snv_filtered_file = $self->snv_filtered_file;
	my $snv_failfiltered_file = $self->snv_failfiltered_file;
	my $snv_annotation_file = $self->snv_annotation_file;
	my $indel_file = $self->indel_file;
	my $indel_filtered_file = $self->indel_filtered_file;
	my $indel_failfiltered_file = $self->indel_failfiltered_file;
	my $indel_annotation_file = $self->indel_annotation_file;

	my $project_name = $self->project_name;
	my $center = $self->center;
	my $ref_build = $self->build;
	my $sequence_phase = $self->sequence_phase;
	my $sequence_source = $self->sequence_source;
	my $sequencer = $self->sequencer;

	my $quality = 0;

	## Open the outfile ##
	open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $date = $year.$mon.$mday;
	print OUTFILE "##fileformat=VCFv4.0"."\n";
	print OUTFILE "##fileDate=$date"."\n";
#	print OUTFILE "##source=myImputationProgramV3.1"."\n";
	print OUTFILE "##reference=Build$ref_build"."\n";
#	print OUTFILE "##phasing=partial"."\n";
	print OUTFILE "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"."\n";   #####DP combined depth across samples, e.g. DP=154
#	print OUTFILE "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">"."\n"; #####AF allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
	print OUTFILE "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">"."\n"; #####DB dbSNP membership
	print OUTFILE "##INFO=<ID=VALIDATED,Number=0,Type=Flag,Description=\"Validation Status\">"."\n"; #####VALIDATED validated by follow-up experiment
#	print OUTFILE "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">"."\n"; #####H2 membership in hapmap2
#	print OUTFILE "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"."\n";
#	print OUTFILE "##INFO=<ID=BQ,Number=1,Type=String,Description=\"RMS base quality at this position\">"."\n"; #BQ RMS base quality at this position
	print OUTFILE "##INFO=<ID=AB,Number=1,Type=Float,Description=\"Allele Balance for hets (ref/(ref+alt))\">"."\n";

	print OUTFILE "##INFO=<ID=Dels,Number=1,Type=Float,Description=\"Fraction of Reads Containing Spanning Deletions\">"."\n"; #not set for varscan output
	print OUTFILE "##INFO=<ID=HRun,Number=1,Type=Integer,Description=\"Largest Contiguous Homopolymer Run of Variant Allele In Either Direction\">"."\n"; #not set for varscan output
	print OUTFILE "##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description=\"Consistency of the site with two (and only two) segregating haplotypes\">"."\n"; #not set for varscan output
	print OUTFILE "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">"."\n"; #not set for varscan output #MQ RMS mapping quality, e.g. MQ=52
	print OUTFILE "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">"."\n"; #not set for varscan output
	print OUTFILE "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">"."\n"; #not set for varscan output

	print OUTFILE "##INFO=<ID=EID,Number=1,Type=Float,Description=\"Entrez ID\">"."\n";
	print OUTFILE "##INFO=<ID=MT,Number=0,Type=Float,Description=\"Mutation Type\">"."\n";


	print OUTFILE "##FILTER=<ID=Filter,Description=\"Failed the FilterFalsePositives or the IndelFilter\">"."\n";

	print OUTFILE "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"."\n";
	print OUTFILE "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"."\n";
	print OUTFILE "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"."\n";
	print OUTFILE "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"."\n";

##FORMAT=<ID=GL,Number=3,Type=Float,Description="Log-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic">

	print OUTFILE "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_name\n";

	## Load dbSNPs ##
	
	my %dbsnp_rs_ids = load_dbsnps($dbsnp_file);

	## Load strandfilter ##
		
	my %strandfilter_lines = load_strandfilter($snv_filtered_file, $snv_failfiltered_file);

	## Load the SNVs ##
		
	my %snvs = load_mutations($variant_file);
	my %snv_annotation = load_annotation($snv_annotation_file);

	## Load the homo vars ##

	my %homo_var_calls = %{&load_homozygous($varscan_file, $samtools_file)};

	foreach my $key (sort byChrPos keys %snvs) {

		my ($chromosome, $chr_start, $chr_stop, $ref, $var, $var_pct, $quality, $depth, $VCF) = split(/\t/, $key);
		if ($VCF == 1 && $strandfilter_lines{$key} eq 'PASS') {
			print OUTFILE "$snvs{$key}\n";
			next;
		}
		my $snv = $snvs{$key};
		my $strand = "+";

		$var_pct =~ s/%//;
		my $ref_pct = (100 - $var_pct);
		
		my $AB = "AB=$ref_pct";
		my $DP = "DP=$depth";

		my @temp = split("\t", $snv);
		$key = "$chromosome\t$chr_start\t$chr_stop\t$ref\t$var";
		if($snv_annotation{$key}) {
			my @annotation = split(/\t/, $snv_annotation{$key});
			my $gene_name = $annotation[6];
			my $tumor_gt_allele1 = $annotation[3];
			my $tumor_gt_allele2 = $annotation[4];

			## Get the gene ID ##
			my $gene_id = 0;

			my @ea = Genome::Site::TGI::EntityAlias->get(alias => "$gene_name", alias_source => "HUGO", entity_type_name => "gene sequence tag");
					
			if(@ea) {
				my @tags = Genome::Site::TGI::SequenceTag->get(stag_id => [ map {$_->entity_id} @ea ]);
				if(@tags) {
					$gene_id = $tags[0]->ref_id;						
				}
			}

			my $trv_type = $annotation[13];
			my $mutation_type = trv_to_mutation_type($trv_type);

			my $MT = "MT=$mutation_type";
			##Get Strandfilter Status
			my $strandfilter_status = $strandfilter_lines{$key};

			## Get dbSNP Status
			my $db = 0;
			my $dbsnp_rs = ".";

			if($dbsnp_rs_ids{$key}) {
				$dbsnp_rs = $dbsnp_rs_ids{$key};
				$db = "DB";
			}

			my $EID = "EID=$gene_id";

			my $info_col = "$AB;$DP;$EID;$MT";
			if ($db) {
				$info_col .= ";$db";
			}

			my $GT = "0/1";
			if (defined $homo_var_calls{$key}) {
				$GT = "1/1";
			}
			my $format_col = "GT\:DP"."\t"."$GT\:$depth"; #GT:DP:GL:GQ	1/1:21:-65.27,-6.34,-0.03:63.09
			print OUTFILE "$chromosome\t$chr_start\t$dbsnp_rs\t$tumor_gt_allele1\t$tumor_gt_allele2\t$quality\t$strandfilter_status\t$info_col\t$format_col\n";

=cut
$chr_stop
$gene_name
$center
$strand
$sequence_phase
$sequence_source
$bam_file
$sequencer
@annotation
=cut
	
		}
		else {
			warn "No annotation for $key in $snv_annotation_file!\n";
		}
	}

	## Write indels to file ##
	my %indels_written = ();
	
	## Load the Indels ##
	my %indels = load_mutations($indel_file);
	my %indel_annotation = load_annotation($indel_annotation_file);
	
	## Load strandfilter ##
	%strandfilter_lines = load_strandfilter($indel_filtered_file, $indel_failfiltered_file);

	foreach my $key (sort byChrPos keys %indels) {
		my ($chromosome, $chr_start, $chr_stop, $ref, $var, $var_pct, $quality, $depth, $VCF) = split(/\t/, $key);
		if ($VCF == 1 && $strandfilter_lines{$key} eq 'PASS') {
			print OUTFILE "$indels{$key}\n";
			next;
		}
		my $indel = $indels{$key};
		my $strand = "+";

		$var_pct =~ s/%//;
		my $ref_pct = (100 - $var_pct);
		
		my $AB = "AB=$ref_pct";
		my $DP;
		if ($depth > 0) {
			$DP = "DP=$depth";
		}
		else {
			$DP = "DP=NA";
		}

		my $variant_type = "Unknown";
				
		if($ref eq "0" || $ref eq "-" || length($var) > 1) {
			$variant_type = "INS";
		}
		else {
			$variant_type = "DEL";
		}
				
		my @temp = split("\t", $indel);
		$key = "$chromosome\t$chr_start\t$chr_stop\t$ref\t$var";
		if($indel_annotation{$key}) {
			my @annotation = split(/\t/, $indel_annotation{$key});

			my $tumor_gt_allele1 = $annotation[3];
			my $tumor_gt_allele2 = $annotation[4];

			my $gene_name = $annotation[6];
	
			## Get the gene ID ##
			my $gene_id = 0;
	
			my @ea = Genome::Site::TGI::EntityAlias->get(alias => "$gene_name", alias_source => "HUGO", entity_type_name => "gene sequence tag");
						
			if(@ea) {
				my @tags = Genome::Site::TGI::SequenceTag->get(stag_id => [ map {$_->entity_id} @ea ]);
				if(@tags) {
					$gene_id = $tags[0]->ref_id;						
				}
			}
		
			my $trv_type = $annotation[13];
			my $mutation_type = trv_to_mutation_type($trv_type); ####NOT USED YET
			my $MT = "MT=$mutation_type";

			##Get Strandfilter Status
			my $strandfilter_status = $strandfilter_lines{$key};

			## Get dbSNP Status 	
			my $dbsnp_rs = ".";
			my $dbsnp_status = "unknown";

			my $indel_key = "$chromosome\t$chr_start\t$chr_stop\t$variant_type";
			$indels_written{$indel_key} = 1;

			my $EID = "EID=$gene_id";

			my $info_col = "$AB;$DP;$EID;$MT";
			my $format_col = "";#"GT\:DP"."\t"."$GT\:$depth"; #GT:DP:GL:GQ	1/1:21:-65.27,-6.34,-0.03:63.09
			print OUTFILE "$chromosome\t$chr_start\t$dbsnp_rs\t$tumor_gt_allele1\t$tumor_gt_allele2\t$quality\t$strandfilter_status\t$info_col\t$format_col\n";
		}
		else {
			warn "No annotation for $key in $indel_annotation_file!\n";
		}
	}

	return 1;
}


################################################################################################
# SUBS
#
################################################################################################

sub load_mutations
{  
	my $variant_file = shift(@_);
	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	my %mutations = ();

	while (my $line = <$input>) {
		chomp($line);
		$lineCounter++;
		my $VCF = 0;
		(my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var, my @other) = split(/\t/, $line);
		my $var_pct;
		my $germ_likelihood;
		my $depth;
		if ($other[2] =~ m/%/) {
			$var_pct = $other[2];
			if ($other[7] == 0) {
				print "$other[7]\n$line\n";
				$germ_likelihood = 256;
			}
			else {
				$germ_likelihood = (-10 * log10($other[7]));
			}
			$depth = ($other[0] + $other[1]);
		}
		elsif ($other[0] =~ m/OBS_COUNTS/) {
			(my @titles) = split(/:/, $other[0]);
			(my @freqs) = split(/\//, $titles[1]);
			$var_pct = ($freqs[0] / $freqs[2]);
			$depth = $freqs[2];
		}
		elsif ($other[2] =~ m/HaplotypeScore/) {
			$VCF = 1;
			$germ_likelihood = $other[0];
			($var_pct) = $other[2] =~ m/AB=(\d+\.\d+);/;
			($depth) = $other[2] =~ m/DP=(\d+\.\d+);/;
		}
		else {
			if ($ref eq '-' || $var eq '-') {
				$germ_likelihood = $other[1];
				my $refcount = $other[7];
				my $varcount = $other[6];
				$depth = ($varcount + $refcount);
				if ($depth > 0) { 
					$var_pct = ($varcount / $depth);
				}
				else {
					$var_pct = 0;
				}
			}
			else {
				$germ_likelihood = $other[1];
				my @freqs = split(//, $other[4]);
				my $refcount = my $varcount= 0;
				foreach my $base (@freqs) {
					if ($base eq '.' || $base eq ',') {
						$refcount++;
					}
					elsif ($base =~ m/$var/i) {
						$varcount++;
					}
				}
				$depth = ($varcount + $refcount);
				if ($depth > 0) { 
					$var_pct = ($varcount / $depth);
				}
				else {
					$var_pct = 0;
				}
			}
		}


	
		$var = Genome::Model::Tools::Capture::iupac_to_base($ref, $var);
	
		my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var, $var_pct, $germ_likelihood, $depth, $VCF);
		
		$mutations{$key} = $line;		
	}
	
	close($input);


	return(%mutations);
}


sub load_annotation
{  
	my $variant_file = shift(@_);
	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	my %annotation = ();

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);

		my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var);
		
		$annotation{$key} = $line;		
	}
	
	close($input);


	return(%annotation);
}

sub load_dbsnps
{  
	my $variant_file = shift(@_);
	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

#	print "Parsing $variant_file\n";

	my %mutations = ();

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		my $dbsnp_rs_id = $lineContents[$numContents - 1];
		my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var);
		
		$mutations{$key} = $dbsnp_rs_id;
	}
	
	close($input);

	return(%mutations);
}

sub load_strandfilter
{  
	my $strandfilter_file = shift(@_);
	my $strandfilter_junk_file = shift(@_);
	my $strandfilter = new FileHandle ($strandfilter_file);
	my $strandfilter_junk = new FileHandle ($strandfilter_junk_file);

	my $lineCounter = 0;
	my $lineCounter2 = 0;
	my %mutations = ();

	while (my $line = <$strandfilter>)
	{
		chomp($line);
		$lineCounter++;
		
		(my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
		my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var);
		$mutations{$key} = "PASS";
	}
	close($strandfilter);
	while (my $line = <$strandfilter_junk>)
	{
		chomp($line);
		$lineCounter2++;
		
		(my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var, my @else) = split(/\t/, $line);
		my $fail_reason = pop(@else);

		my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var);
		$mutations{$key} = "$fail_reason";
	}
	close($strandfilter_junk);

	return(%mutations);
}

sub load_homozygous
{  
	my (@files) = @_;
	my %homo_var_calls;
	foreach my $file (@files) {
		my $variant_file = new FileHandle ($file);
		while (my $line = <$variant_file>) {
			chomp($line);
			my ($Chrom,$Position,$Ref,$Cons) = split(/\t/, $line);
			my $variant = "$Chrom\t$Position\t$Position\t$Ref\t$Cons";
			if ($Cons =~ m/A|C|T|G/) {
				$homo_var_calls{$variant}++;
			}
		}
	}
	return \%homo_var_calls;
}

#############################################################
# ParseBlocks - takes input file and parses it
#
#############################################################

sub trv_to_mutation_type
{
	my $trv_type = shift(@_);
	
	return("Missense_Mutation") if($trv_type eq "missense");	
	return("Nonsense_Mutation") if($trv_type eq "nonsense" || $trv_type eq "nonstop");	
	return("Silent") if($trv_type eq "silent");		
	return("Splice_Site_SNP") if($trv_type eq "splice_site");
	return("Splice_Site_Indel") if($trv_type eq "splice_site_del");		
	return("Splice_Site_Indel") if($trv_type eq "splice_site_ins");		
	return("Frame_Shift_Del") if($trv_type eq "frame_shift_del");		
	return("Frame_Shift_Ins") if($trv_type eq "frame_shift_ins");		
	return("In_Frame_Del") if($trv_type eq "in_frame_del");		
	return("In_Frame_Ins") if($trv_type eq "in_frame_ins");		
	return("RNA") if($trv_type eq "rna");		

	warn "Unknown mutation type $trv_type\n";
	return("Unknown");
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub byChrPos
{
	my ($chrom_a, $pos_a) = split(/\t/, $a);
	my ($chrom_b, $pos_b) = split(/\t/, $b);
	
	$chrom_a =~ s/X/23/;
	$chrom_a =~ s/Y/24/;
	$chrom_a =~ s/MT/25/;
	$chrom_a =~ s/[^0-9]//g;

	$chrom_b =~ s/X/23/;
	$chrom_b =~ s/Y/24/;
	$chrom_b =~ s/MT/25/;
	$chrom_b =~ s/[^0-9]//g;

	$chrom_a <=> $chrom_a
	or
	$pos_a <=> $pos_b;
}

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

1;












