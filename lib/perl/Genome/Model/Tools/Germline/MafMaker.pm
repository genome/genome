
package Genome::Model::Tools::Germline::MafMaker;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# GermlinePipelineMafMaker - Generate MAF File after germline pipeline is done
#					
#	AUTHOR:		Will Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	02/01/2011 by W.S.
#	MODIFIED:	02/01/2011 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it
use Genome::Model::Tools::Capture::Helpers qw(
    byChrPos
    iupac_to_base
);

## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Germline::MafMaker {
    is => 'Command',                       

    has => [                                # specify the command's single-value properties (parameters) <--- 
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
    "Generate MAF File after germline pipeline is done"                 
}

sub help_synopsis {
    return <<EOS
Generate MAF File, Get dbsnp output, and strandfilter -- for GERMLINE events
EXAMPLE:	gmt capture germline-pipeline-maf-maker
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

	my $bam_file = $self->bam_file;
	my $variant_file = $self->variant_file;
	my $dbsnp_file = $self->dbsnp_file;
	my $snv_filtered_file = $self->snv_filtered_file;
	my $snv_failfiltered_file = $self->snv_failfiltered_file;
	my $snv_annotation_file = $self->snv_annotation_file;
	my $indel_file = $self->indel_file;
	my $indel_filtered_file = $self->indel_filtered_file;
	my $indel_failfiltered_file = $self->indel_failfiltered_file;
	my $indel_annotation_file = $self->indel_annotation_file;
	my $output_file = $self->output_file;

	my $project_name = $self->project_name;
	my $center = $self->center;
	my $ref_build = $self->build;
	my $sequence_phase = $self->sequence_phase;
	my $sequence_source = $self->sequence_source;
	my $sequencer = $self->sequencer;

	my $build_id = $self->build_id;
	my $build = Genome::Model::Build->get(build_id => $build_id);
	my $sample_name = $build->subject_name;

	## Open the outfile ##
	open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
	print OUTFILE join("\t", "Hugo_Symbol","Entrez_Gene_Id","Center","NCBI_Build","Chromosome","Start_position","End_position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2", "dbSNP_RS","dbSNP_Val_Status","Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode","Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2","Tumor_Validation_Allele1","Tumor_Validation_Allele2","Match_Norm_Validation_Allele1","Match_Norm_Validation_Allele2", "Verification_Status","Validation_Status","Mutation_Status","Validation_Method","Sequencing_Phase","Sequence_Source","Score","BAM_file","Sequencer","chromosome_name_WU","start_WU","stop_WU","reference_WU","variant_WU", "type_WU","gene_name_WU","transcript_name_WU","transcript_species_WU","transcript_source_WU","transcript_version_WU","strand_WU","transcript_status_WU","trv_type_WU","c_position_WU","amino_acid_change_WU","ucsc_cons_WU", "domain_WU","all_domains_WU","deletion_substructures_WU","transcript_error_WU") . "\n";

	## Load dbSNPs ##
	
	my %dbsnp_rs_ids = load_dbsnps($dbsnp_file);

	## Load strandfilter ##
		
	my %strandfilter_lines = load_strandfilter($snv_filtered_file, $snv_failfiltered_file);
	
	## Load the SNVs ##
		
	my %snvs = load_mutations($variant_file);
	my %snv_annotation = load_annotation($snv_annotation_file);

	foreach my $key (sort byChrPos keys %snvs) {
		$stats{'tier1_snvs'}++;
			
		my ($chromosome, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $key);
		my $snv = $snvs{$key};
		my $strand = "+";
			
		my @temp = split("\t", $snv);
			
		if($snv_annotation{$key}) {
			$stats{'tier1_snvs_with_annotation'}++;
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

			##Get Strandfilter Status
			my $strandfilter_status = $strandfilter_lines{$key};

			## Get dbSNP Status 	
			my $dbsnp_rs = "novel";
			my $dbsnp_status = "unknown";

			if($dbsnp_rs_ids{$key}) {
				$dbsnp_rs = $dbsnp_rs_ids{$key};
				$dbsnp_status = "unknown";
			}

			print OUTFILE join("\t", $gene_name,$gene_id,$center,$ref_build,$chromosome,$chr_start,$chr_stop,$strand,$mutation_type,"SNP",$ref,$tumor_gt_allele1,$tumor_gt_allele2,$dbsnp_rs,$dbsnp_status,$sample_name,$sample_name,$ref,$ref,"","","","",$strandfilter_status,"Unknown","Germline",$sequence_phase,$sequence_source,"","1",$bam_file,$sequencer, @annotation) . "\n";

			$stats{'tier1_snvs_written'}++;
	
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
		$stats{'tier1_indels'}++;
			
		my ($chromosome, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $key);
		my $indel = $indels{$key};
		my $strand = "+";
				
		my $variant_type = "Unknown";
				
		if($ref eq "0" || $ref eq "-" || length($var) > 1) {
			$variant_type = "INS";
		}
		else {
			$variant_type = "DEL";
		}
				
		my @temp = split("\t", $indel);
				
		if($indel_annotation{$key}) {
			$stats{'tier1_indels_with_annotation'}++;
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
			my $mutation_type = trv_to_mutation_type($trv_type);

			##Get Strandfilter Status
			my $strandfilter_status = $strandfilter_lines{$key};

			## Get dbSNP Status 	
			my $dbsnp_rs = "novel(indel)";
			my $dbsnp_status = "unknown";

			my $indel_key = "$chromosome\t$chr_start\t$chr_stop\t$variant_type";
			$indels_written{$indel_key} = 1;

			print OUTFILE join("\t", $gene_name,$gene_id,$center,$ref_build,$chromosome,$chr_start,$chr_stop,$strand,$mutation_type,$variant_type,$ref,$tumor_gt_allele1,$tumor_gt_allele2,$dbsnp_rs,$dbsnp_status,$sample_name,$sample_name,$ref,$ref,"","","","",$strandfilter_status,"Unknown","Germline",$sequence_phase,$sequence_source,"","1",$bam_file,$sequencer, @annotation) . "\n";				
			$stats{'tier1_indels_written'}++;
		}
		else {
			warn "No annotation for $key in $indel_annotation_file!\n";
		}
	}

	$stats{'tier1_snvs_not_included'} = 0 if(!$stats{'tier1_snvs_not_included'});

	$stats{'tier1_snvs'} = 0 if(!$stats{'tier1_snvs'});
	$stats{'tier1_snvs_not_included'} = 0 if(!$stats{'tier1_snvs_not_included'});
	$stats{'tier1_snvs_with_annotation'} = 0 if(!$stats{'tier1_snvs_with_annotation'});
	$stats{'tier1_snvs_written'} = 0 if(!$stats{'tier1_snvs_written'});
	$stats{'tier1_indels'} = 0 if(!$stats{'tier1_indels'});
	$stats{'tier1_indels_not_included'} = 0 if(!$stats{'tier1_indels_not_included'});
	$stats{'tier1_indels_with_annotation'} = 0 if(!$stats{'tier1_indels_with_annotation'});
	$stats{'tier1_indels_written'} = 0 if(!$stats{'tier1_indels_written'});	

	print $stats{'tier1_snvs'} . " tier 1 SNVs\n";
	print $stats{'tier1_snvs_not_included'} . " not included in target list\n";
	print $stats{'tier1_snvs_with_annotation'} . " met criteria and had annotation\n";
	print $stats{'tier1_snvs_written'} . " were written to MAF file\n\n";

	print $stats{'tier1_indels'} . " tier 1 Indels\n";
	print $stats{'tier1_indels_not_included'} . " not included in target list\n";
	print $stats{'tier1_indels_with_annotation'} . " met criteria and had annotation\n";
	print $stats{'tier1_indels_written'} . " were written to MAF file\n";

}

1;


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

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
	
		$var = iupac_to_base($ref, $var);
	
		my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var);
		
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

	print "$lineCounter dbSNPs loaded\n";

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
		$mutations{$key} = "Strandfilter_Passed";
	}
	close($strandfilter);
	while (my $line = <$strandfilter_junk>)
	{
		chomp($line);
		$lineCounter2++;
		
		(my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
		my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var);
		$mutations{$key} = "Strandfilter_Failed";
	}
	close($strandfilter_junk);

	print "$lineCounter strandfilter_passed\n";
	print "$lineCounter2 strandfilter_failed\n";

	return(%mutations);
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


1;
