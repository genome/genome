package Genome::Model::Tools::Capture::ParsingScript;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ParsingScript - Take glfSomatic and Varscan outputs and merge with Annotator Output
#					
#	AUTHOR:		Will Schierding
#
#	CREATED:	1/09/2009
#
#	NOTES:	
#			
#####################################################################################################################################

#__STANDARD PERL PACKAGES
   use warnings;
   use strict;
   use FileHandle;
   use Getopt::Long;

#__SPECIAL GENOME CENTER PACKAGES
   use Genome::Model::Tools::Capture::PipelineParser;
#   use MG::IO::Parse::Cosmic;
   use Genome;

class Genome::Model::Tools::Capture::ParsingScript {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		annotation_file	=> { is => 'Text', doc => "Annotation Output File", is_optional => 0 },
		merged_snv	=> { is => 'Text', doc => "Overlapping Calls from glfSomatic and Varscan", is_optional => 0 },
		varscan_unique_snv	=> { is => 'Text', doc => "Unique Calls from Varscan", is_optional => 0 },
		glf_unique_snv	=> { is => 'Text', doc => "Unique Calls from glfSomatic", is_optional => 0 },
		merged_indel	=> { is => 'Text', doc => "Overlapping Calls from glfSomatic and Varscan", is_optional => 0 },
		varscan_unique_indel	=> { is => 'Text', doc => "Unique Calls from Varscan", is_optional => 0 },
		glf_unique_indel	=> { is => 'Text', doc => "Unique Calls from glfSomatic", is_optional => 0 },
		out_file	=> { is => 'Text', doc => "Output file with Merged Calls and Software Status" , is_optional => 0},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Take glfSomatic and Varscan outputs and merge with Annotator Output"                 
}

sub help_synopsis {
    return <<EOS
Take glfSomatic and Varscan outputs and merge with Annotator Output
EXAMPLE:	gmt capture build-models ...
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
	my $annotation_file = $self->annotation_file;
	my $merged_snv = $self->merged_snv;
	my $varscan_unique_snv = $self->varscan_unique_snv;
	my $glf_unique_snv = $self->glf_unique_snv;
	my $merged_indel = $self->merged_indel;
	my $varscan_unique_indel = $self->varscan_unique_indel;
	my $glf_unique_indel = $self->glf_unique_indel;
	my $out_file = $self->out_file;

# Open Input
unless (open(ANNOT_FILE_MERGED_SNV,"<$merged_snv")) {
   die "Could not open input file '$merged_snv' for writing";
  }
unless (open(ANNOT_FILE_GLF_SNV,"<$glf_unique_snv")) {
   die "Could not open input file '$glf_unique_snv' for writing";
  }
unless (open(ANNOT_FILE_VARSCAN_SNV,"<$varscan_unique_snv")) {
   die "Could not open input file '$varscan_unique_snv' for writing";
  }
unless (open(ANNOT_FILE_MERGED_INDEL,"<$merged_indel")) {
   die "Could not open input file '$merged_indel' for writing";
  }
unless (open(ANNOT_FILE_GLF_INDEL,"<$glf_unique_indel")) {
   die "Could not open input file '$glf_unique_indel' for writing";
  }
unless (open(ANNOT_FILE_VARSCAN_INDEL,"<$varscan_unique_indel")) {
   die "Could not open input file '$varscan_unique_indel' for writing";
  }

my $header = 1;
if ($header) {
   my $head = <ANNOT_FILE_MERGED_SNV>;
#  $head = <ANNOT_FILE_GLF_SNV>;
   $head = <ANNOT_FILE_VARSCAN_SNV>;
   $head = <ANNOT_FILE_MERGED_INDEL>;
#  $head = <ANNOT_FILE_GLF_INDEL>;
   $head = <ANNOT_FILE_VARSCAN_INDEL>;
}

my %columns_merged_snv;
while( my $additions = <ANNOT_FILE_MERGED_SNV> ) {
   my($chrom,$position,$ref,$var,$normal_reads1,$normal_reads2,$normal_var_freq,$normal_gt,$tumor_reads1,$tumor_reads2,$tumor_var_freq,$tumor_gt,$somatic_status,$pValue) = split(/\t/, $additions);
   my $Chromosome = $chrom;
   my $Start_position = $position;
   my $End_position = $position;
   my $Reference_Allele = $ref;
   my $Tumor_Seq_Allele1 = $var;
   my $merger = "$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1";
   $columns_merged_snv{$merger} = "$normal_reads1\t$normal_reads2\t$normal_var_freq\t$normal_gt\t$tumor_reads1\t$tumor_reads2\t$tumor_var_freq\t$tumor_gt\t$somatic_status\t$pValue";
}

close ANNOT_FILE_MERGED_SNV;

my %columns_glf_snv;
while( my $additions = <ANNOT_FILE_GLF_SNV> ) {
   my ($Chromosome,$position,$reference,$tumor_consensus_call,$somatic_score,$tumor_consensus_quality_score,$tumor_snp_quality_score,$Root_Mean_Square_mapping_quality,$Tumor_Reads,$Normal_reads) = split(/\t/, $additions);
#   my $Chromosome = $Chromosome;
   my $Start_position = $position;
   my $End_position = $position;
   my $Reference_Allele = $reference;
   my $Tumor_Seq_Allele1 = $tumor_consensus_call;
   my $merger = "$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1";
   $columns_glf_snv{$merger} = "$somatic_score\t$tumor_consensus_quality_score\t$tumor_snp_quality_score\t$Root_Mean_Square_mapping_quality\t$Tumor_Reads\t$Normal_reads";
}

close ANNOT_FILE_GLF_SNV;

my %columns_varscan_snv;
while( my $additions = <ANNOT_FILE_VARSCAN_SNV> ) {
   my($chrom,$position,$ref,$var,$normal_reads1,$normal_reads2,$normal_var_freq,$normal_gt,$tumor_reads1,$tumor_reads2,$tumor_var_freq,$tumor_gt,$somatic_status,$pValue) = split(/\t/, $additions);
   my $Chromosome = $chrom;
   my $Start_position = $position;
   my $End_position = $position;
   my $Reference_Allele = $ref;
   my $Tumor_Seq_Allele1 = $var;
   my $merger = "$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1";
   $columns_varscan_snv{$merger} = "$normal_reads1\t$normal_reads2\t$normal_var_freq\t$normal_gt\t$tumor_reads1\t$tumor_reads2\t$tumor_var_freq\t$tumor_gt\t$somatic_status\t$pValue";
}

close ANNOT_FILE_VARSCAN_SNV;

my %columns_merged_indel;
while( my $additions = <ANNOT_FILE_MERGED_INDEL> ) {
   my($chrom,$position,$ref,$var,$normal_reads1,$normal_reads2,$normal_var_freq,$normal_gt,$tumor_reads1,$tumor_reads2,$tumor_var_freq,$tumor_gt,$somatic_status,$pValue) = split(/\t/, $additions);
   my $Chromosome = $chrom;
   my $Start_position = $position;
   my $End_position = $position;
   my $Reference_Allele = $ref;
   my $Tumor_Seq_Allele1 = $var;
   my $merger = "$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1";
   $columns_merged_indel{$merger} = "$normal_reads1\t$normal_reads2\t$normal_var_freq\t$normal_gt\t$tumor_reads1\t$tumor_reads2\t$tumor_var_freq\t$tumor_gt\t$somatic_status\t$pValue";
}

close ANNOT_FILE_MERGED_INDEL;

my %columns_glf_indel;
while( my $additions = <ANNOT_FILE_GLF_INDEL> ) {
   my ($CHR,$POS,$always_a_star,$SOMATIC_SCORE,$INDEL1_SEQUENCE,$INDEL2_SEQUENCE,$INDEL1_LENGTH,$INDEL2_LENGTH,$tumor_something_something,$tumor_consensus_quality,$tumor_reference_quality,$tumor_max_mapping_quality,$tumor_number_of_reads_at_position,$tumor_reads_1,$tumor_reads_2,$tumor_reads_ambi,$tumor_reads_anti,$tumor_min_likelihood,$tumor_likelihood1,$tumor_likelihood2,$tumor_likelihood3,$normal_something_something,$normal_consensus_quality,$normal_reference_quality,$normal_max_mapping_quality,$normal_number_of_reads_at_position,$normal_reads_1,$normal_reads_2,$normal_reads_ambi,$normal_reads_anti,$normal_min_likelihood,$normal_likelihood1,$normal_likelihood2,$normal_likelihood3,$LIBS_INDEL1,$LIBS_INDEL2) = split(/\t/, $additions);

   my $Chromosome = $CHR;
   my $Start_position = $POS;
   my $End_position = $POS;
############# THIS MAY NOT BE RIGHT, NEXT TWO LINES####################
   my $Reference_Allele = $INDEL1_SEQUENCE;
   my $Tumor_Seq_Allele1 = $INDEL2_SEQUENCE;
   my $merger = "$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1";
   $columns_glf_indel{$merger} = "$SOMATIC_SCORE\t$INDEL1_LENGTH\t$INDEL2_LENGTH\t$tumor_something_something\t$tumor_consensus_quality\t$tumor_reference_quality\t$tumor_max_mapping_quality\t$tumor_number_of_reads_at_position\t$tumor_reads_1\t$tumor_reads_2\t$tumor_reads_ambi\t$tumor_reads_anti\t$tumor_min_likelihood\t$tumor_likelihood1\t$tumor_likelihood2\t$tumor_likelihood3\t$normal_something_something\t$normal_consensus_quality\t$normal_reference_quality\t$normal_max_mapping_quality\t$normal_number_of_reads_at_position\t$normal_reads_1\t$normal_reads_2\t$normal_reads_ambi\t$normal_reads_anti\t$normal_min_likelihood\t$normal_likelihood1\t$normal_likelihood2\t$normal_likelihood3\t$LIBS_INDEL1\t$LIBS_INDEL2";
}


close ANNOT_FILE_GLF_INDEL;

my %columns_varscan_indel;
while( my $additions = <ANNOT_FILE_VARSCAN_INDEL> ) {
   my($chrom,$position,$ref,$var,$normal_reads1,$normal_reads2,$normal_var_freq,$normal_gt,$tumor_reads1,$tumor_reads2,$tumor_var_freq,$tumor_gt,$somatic_status,$pValue) = split(/\t/, $additions);
   my $Chromosome = $chrom;
   my $Start_position = $position;
   my $End_position = $position;
   my $Reference_Allele = $ref;
   my $Tumor_Seq_Allele1 = $var;
   my $merger = "$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1";
   $columns_varscan_indel{$merger} = "$normal_reads1\t$normal_reads2\t$normal_var_freq\t$normal_gt\t$tumor_reads1\t$tumor_reads2\t$tumor_var_freq\t$tumor_gt\t$somatic_status\t$pValue";
}

close ANNOT_FILE_VARSCAN_INDEL;
################
#  PROCESSING  #
################

#__SET SOME PARSING PARAMETERS -- UNSURE OF MEANING OF ORIGINAL COMMENTS (MCW)
   my %parse_args = (

   #__PROCESS EVERYTHING INTO A SINGLE STRUCTURE (AND THEN PROCESSED)
      'all' => 1,

#   #__HAVE EVERYTHING PROCESSED INTO A SINGLE STRUCTURE (AND NOT PROCESSED)
#      'no_process' => 1,
   );

my $fh1 = new FileHandle;
my $fh2 = new FileHandle;
# open "file.annotated"
   unless ($fh1->open (qq{$annotation_file})) {
      die "Could not open mutation project file '$annotation_file' for reading";
   }
# output file MAF
   unless (open($fh2,">$out_file")) {
      die "Could not open mutation project file '$out_file' for reading";
   }

my $parser = Genome::Model::Tools::Capture::PipelineParser->new();
my $annotation = $parser->Parse ($fh1, $annotation_file, %parse_args);
foreach my $hugo (keys %{$annotation}) {
# gene name = $hugo
    foreach my $line_num (keys %{$annotation->{$hugo}}) {
        print STDOUT ".";   #report that we are starting a sample (For commandline user feedback)

        my ($line, $aa_change,$transcript,$mstatus,$Variant_Type,$Chromosome,$Start_position,$End_position,$Reference_Allele,$Tumor_Seq_Allele1,$source,$genome,$strand,$trv_type,$c_position,$ucsc_cons,$domain) =
            (
                $annotation->{$hugo}->{$line_num}->{file_line},
                $annotation->{$hugo}->{$line_num}->{AA_CHANGE},
                $annotation->{$hugo}->{$line_num}->{TRANSCRIPT},
                $annotation->{$hugo}->{$line_num}->{MUTATION_STATUS},
                $annotation->{$hugo}->{$line_num}->{VARIANT_TYPE},
                $annotation->{$hugo}->{$line_num}->{CHROMOSOME},
                $annotation->{$hugo}->{$line_num}->{START_POSITION},
                $annotation->{$hugo}->{$line_num}->{END_POSITION},
                $annotation->{$hugo}->{$line_num}->{REFERENCE_ALLELE},
                $annotation->{$hugo}->{$line_num}->{TUMOR_SEQ_ALLELE1},
		$annotation->{$hugo}->{$line_num}->{SOURCE},
		$annotation->{$hugo}->{$line_num}->{GENOME},
		$annotation->{$hugo}->{$line_num}->{STRAND},
		$annotation->{$hugo}->{$line_num}->{TYPE},
		$annotation->{$hugo}->{$line_num}->{C_POSITION},
		$annotation->{$hugo}->{$line_num}->{UCSC},
		$annotation->{$hugo}->{$line_num}->{DOMAIN},
            );
#1,2,3,4,5,6,7,8,14,15,16,17,18,19 
	my $merger = "$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1";

	if ($columns_merged_snv{$merger}) {
	    my $software = "Varscan,glfSomatic";
	    $fh2->print("$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1\t$aa_change\t$transcript\t$mstatus\t$c_position\t$ucsc_cons\t$domain\t$software\t$columns_merged_snv{$merger}\n");
	}
	elsif ($columns_glf_snv{$merger}) {
	    my $software = "glfSomatic";
	    $fh2->print("$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1\t$aa_change\t$transcript\t$mstatus\t$c_position\t$ucsc_cons\t$domain\t$software\t$columns_glf_snv{$merger}\n");
	}
	elsif ($columns_varscan_snv{$merger}) {
	    my $software = "Varscan";
	    $fh2->print("$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1\t$aa_change\t$transcript\t$mstatus\t$c_position\t$ucsc_cons\t$domain\t$software\t$columns_varscan_snv{$merger}\n");
	}
	elsif ($columns_merged_indel{$merger}) {
	    my $software = "Varscan,glfSomatic";
	    $fh2->print("$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1\t$aa_change\t$transcript\t$mstatus\t$c_position\t$ucsc_cons\t$domain\t$software\t$columns_merged_indel{$merger}\n");
	}
	elsif ($columns_glf_indel{$merger}) {
	    my $software = "glfSomatic";
	    $fh2->print("$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1\t$aa_change\t$transcript\t$mstatus\t$c_position\t$ucsc_cons\t$domain\t$software\t$columns_glf_indel{$merger}\n");
	}
	elsif ($columns_varscan_indel{$merger}) {
	    my $software = "Varscan";
	    $fh2->print("$Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_Allele1\t$aa_change\t$transcript\t$mstatus\t$c_position\t$ucsc_cons\t$domain\t$software\t$columns_varscan_indel{$merger}\n");
	}
    }
}
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


1;

