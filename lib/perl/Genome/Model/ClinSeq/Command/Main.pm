package Genome::Model::ClinSeq::Command::Main;

#Written by Malachi Griffith

#Load modules
use strict;
use warnings;
use Genome; 
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use File::Basename;
use Genome::Model::ClinSeq::Util qw(:all);

my $script_dir;
use Cwd 'abs_path';
BEGIN{
  $script_dir = abs_path(File::Basename::dirname(__FILE__) . "/../original-scripts");
  $script_dir .= "/";
}

class Genome::Model::ClinSeq::Command::Main {
  is => 'Command::V2',
  has_input => [
      build => { 
                    is => 'Genome::Model::Build::ClinSeq',
                    id_by => 'build_id',
                    doc => 'Used to pass in the current ID of a Clinseq build (not used when running clinseq.pl directly)',
                  },
      wgs_build => {
                    is => 'Genome::Model::Build::SomaticVariation',
                    doc => 'Whole genome sequence (WGS) somatic variation build',
                    is_optional => 1,
                  },
      exome_build => {
                    is => 'Genome::Model::Build::SomaticVariation',
                    doc => 'Exome capture sequence somatic variation build',
                    is_optional => 1,
                   },
      tumor_rnaseq_build => {
                    is => 'Genome::Model::Build::RnaSeq',
                    doc => 'RNA-seq build for the tumor sample',
                    is_optional => 1,
                   },
      normal_rnaseq_build => {       
                    is => 'Genome::Model::Build::RnaSeq',
                    doc => 'RNA-seq build for the normal sample',
                    is_optional => 1,
                   },
      working_dir => {
                    is => 'Text',
                    doc => 'Directory where a patient subdir will be created',
                  }, 
      common_name => {
                    is => 'Text',
                    doc => "Patient's common name (will be used for the name of a results dir and labeling purposes)",
                  },

  ],
  has_param => [                  
      verbose => {
                    is => 'Number',
                    doc => 'To display more output, set to 1',
                    default_value => 0,
                    valid_values => [0,1],
                  },
  ],
  has_output => [
    # for SummarizeTier1SnvSupport
    wgs_positions_file          => { is => 'FilesystemPath', is_optional => 1 },
    exome_positions_file        => { is => 'FilesystemPath', is_optional => 1 },
    wgs_exome_positions_file    => { is => 'FilesystemPath', is_optional => 1 },
    tumor_fpkm_file             => { is => 'FilesystemPath', is_optional => 1 },

  ],
  doc => "This script attempts to automate the process of running the 'clinseq' pipeline",
};

sub help_synopsis {
    return <<EOS
    genome model clin-seq main --wgs-som-var-data-set='2882504846'  --exome-som-var-data-set='2882505032'  --tumor-rna-seq-data-set='2880794613'  --working-dir=/gscmnt/sata132/techd/mgriffit/hgs/  --common-name='ALL1'
EOS
}

sub help_detail {
    return <<EOS
 This script attempts to automate the process of running the 'clinseq' pipeline

 Input data (one or more of the following)
 1.) Whole genome somatic variation model id
 2.) Whole exome somatic variation model id
 3.) RNA-seq model id
 4.) Whole genome germline variation model id

 Big picture goals.
 1.) Summarize somatic and germline variation for a single tumor/normal pair leveraging WGS and/or Exome data
 2.) Summarize RNA expression events relative to whatever comparison are available for the index patient
 3.) Specifically identify the events most likely to be useful in a clinical context (clinically actionable events) - Missense mutations, amplifications, over-expressed genes
 4.) Generate summary statistics files and figures that will serve as the input for a clincal genomics report to be generated downstream

 See the ClinSeq README.txt for details
EOS
}

sub execute {
    my $self = shift;
    
    
    #Before executing, change the environment variable for R_LIBS to be ''
    #This will force R to use it own local notion of library paths instead of the /gsc/ versions
    #This should work for R installed on the machine /usr/bin/R  OR  a standalone version of R installed by a local user. e.g. /gscmnt/gc2142/techd/tools/R/R-2.14.0/bin/R
    local $ENV{R_LIBS}='';

    #Make sure the right libraries are used (in case someone runs with a perl -I statement).
    my $prefix = UR::Util->used_libs_perl5lib_prefix;
    local $ENV{PERL5LIB} = $prefix . ':' . $ENV{PERL5LIB};

    # redirect STDOUT to STDERR
    open(OLD, ">&STDOUT"); #Save stdout
    open(STDOUT,">&STDERR"); #Redirect stdout to go to stderr
    
    my $result = eval { $self->_execute() };
    my $exception_saved = $@;

    # restore STDOUT
    open(STDOUT, ">&OLD"); #Return stdout to its usual state
    close(OLD);
    if ($exception_saved) {
        die "Exception: $exception_saved";
    }
    if (!$result){
        die "Bad return value from ClinSeq::Command::Main";
    } 
    
    return $result;
}

sub _execute {
  my $self = shift;
  my $clinseq_build_id = $self->build_id; #Build ID of the current clinseq run...
  my $clinseq_build = Genome::Model::Build->get($clinseq_build_id);
  my $wgs_som_var_data_set = $self->wgs_build;
  my $exome_som_var_data_set = $self->exome_build;
  my $tumor_rna_seq_data_set = $self->tumor_rnaseq_build;
  my $normal_rna_seq_data_set = $self->normal_rnaseq_build;
  my $working_dir = $self->working_dir;
  my $common_name = $self->common_name;
  my $verbose = $self->verbose;

  #Get build directories for the three datatypes: $data_paths->{'wgs'}->*, $data_paths->{'exome'}->*, $data_paths->{'tumor_rnaseq'}->*
  my $step = 0;
  $step++; print MAGENTA, "\n\nStep $step. Getting data paths from 'genome' for specified model ids\n", RESET;
  my ($data_paths, $builds) = &getDataDirsAndBuilds('-wgs_som_var_data_set'=>$wgs_som_var_data_set, '-exome_som_var_data_set'=>$exome_som_var_data_set, '-tumor_rna_seq_data_set'=>$tumor_rna_seq_data_set, '-normal_rna_seq_data_set'=>$normal_rna_seq_data_set);

  #Option to remove MT chr snvs/indels
  my $filter_mt = 1;

  #Set flags for each datatype
  my $wgs = exists $builds->{wgs} || 0;
  my $exome = exists $builds->{exome} || 0;
  my $tumor_rnaseq = exists $builds->{tumor_rnaseq} || 0;
  my $normal_rnaseq = exists $builds->{normal_rnaseq} || 0;

  #Check the working dir
  $working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");

  #Determine annotation version used in the analysis by examining input builds of the clinseq model - watch out for conflicting versions!
  $step++; print MAGENTA, "\n\nStep $step. Determining annotation build used in the analysis by examining input builds of the clinseq model", RESET;
  my $annotation_build_id = $self->getAnnotationBuild('-clinseq_build_id'=>$clinseq_build_id);
  print BLUE, "\n\tAnnotation build id = $annotation_build_id", RESET;

  #Set the patient subdir
  $step++; print MAGENTA, "\n\nStep $step. Checking/creating the working dir for this patient", RESET;
  my $patient_dir = $working_dir . '/' . $common_name . '/';

  #Get Entrez and Ensembl data for gene name mappings
  my $entrez_ensembl_data = &loadEntrezEnsemblData();

  #Reference annotations - Extra annotation files not currently part of the APIPE system...
  #TODO: Create official versions of these data on allocated disk
  my $clinseq_annotations_dir = "/gscmnt/sata132/techd/mgriffit/reference_annotations/";

  #Directory of gene lists for various purposes
  my $gene_symbol_lists_dir = $clinseq_annotations_dir . "GeneSymbolLists/";
  $gene_symbol_lists_dir = &checkDir('-dir'=>$gene_symbol_lists_dir, '-clear'=>"no");

  #Import a set of gene symbol lists (these files must be gene symbols in the first column, .txt extension, tab-delimited if multiple columns, one symbol per field, no header)
  #Different sets of genes list could be used for different purposes
  #Fix gene names as they are being imported
  $step++; print MAGENTA, "\n\nStep $step. Importing gene symbol lists (from $gene_symbol_lists_dir)", RESET;
  my $symbol_list_names = &importSymbolListNames('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-verbose'=>$verbose);
  my $master_list = $symbol_list_names->{master_list};
  my @symbol_list_names = sort {$master_list->{$a}->{order} <=> $master_list->{$b}->{order}} keys %{$master_list};
  my $gene_symbol_lists = &importGeneSymbolLists('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-symbol_list_names'=>\@symbol_list_names, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);

  #Create a hash for storing output files as they are created
  my %out_paths;
  my $out_paths = \%out_paths;

  #Create a summarized file of SNVs for: WGS, exome, and WGS+exome merged
  #Grab the gene name used in the 'annotation.top' file, but grab the AA changes from the '.annotation' file
  #Fix the gene name if neccessary...
  $step++; print MAGENTA, "\n\nStep $step. Summarizing SNVs and Indels", RESET;
  if ($wgs || $exome){
      $self->importSNVs('-data_paths'=>$data_paths, '-out_paths'=>$out_paths, '-patient_dir'=>$patient_dir, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>$verbose, '-filter_mt'=>$filter_mt);
  }

  #TODO: More comprehensive processing of SNVs and InDels
  #Import SNVs and Indels in a more complete form
  #Make copies of Tier1,2,3 files
  #Make a master list of all distinct SNV/Indels - add a column that classifies them by Tier - Exclude Mt positions
  #For the master list of SNVs (also for INDELS) get the BAM read counts for all positions in tumor and normal
  #Add dbSNP annotations to the SNVs/InDELs
  #Add 1000 genomes annotations to the SNVs/InDELs

  #Run RNA-seq analysis on the RNA-seq data (if available)
  my $rnaseq_dir = $patient_dir . 'rnaseq/';
  if ($tumor_rnaseq){
    my $tumor_rnaseq_dir = $rnaseq_dir . 'tumor/';

    #Perform QC, splice site, and junction expression analysis using the Tophat output
    $step++; print MAGENTA, "\n\nStep $step. Summarizing RNA-seq tophat alignment results, splice sites and junction expression - Tumor", RESET;

    my $tumor_rnaseq_build = $builds->{tumor_rnaseq};
    my $tumor_rnaseq_build_dir = $tumor_rnaseq_build->data_directory;
    my $junctions_dir = $tumor_rnaseq_build_dir . "/junctions/";
    if (-e $junctions_dir){
      my $results_dir = &createNewDir('-path'=>$tumor_rnaseq_dir, '-new_dir_name'=>'tophat_junctions_absolute', '-silent'=>1);
      my $cp_cmd = "cp -r $junctions_dir" . "* $results_dir";
      Genome::Sys->shellcmd(cmd=>$cp_cmd);
    }
  }
  if ($normal_rnaseq){
    my $normal_rnaseq_dir = &createNewDir('-path'=>$rnaseq_dir, '-new_dir_name'=>'normal', '-silent'=>1);

    #Perform QC, splice site, and junction expression analysis using the Tophat output
    $step++; print MAGENTA, "\n\nStep $step. Summarizing RNA-seq tophat alignment results, splice sites and junction expression - Normal", RESET;

    my $normal_rnaseq_build = $builds->{normal_rnaseq};
    my $normal_rnaseq_build_dir = $normal_rnaseq_build->data_directory;
    my $junctions_dir = $normal_rnaseq_build_dir . "/junctions/";
    if (-e $junctions_dir){
      my $results_dir = &createNewDir('-path'=>$normal_rnaseq_dir, '-new_dir_name'=>'tophat_junctions_absolute', '-silent'=>1);
      my $cp_cmd = "cp -r $junctions_dir" . "* $results_dir";
      Genome::Sys->shellcmd(cmd=>$cp_cmd);
    }
  }

  #Annotate gene lists to deal with commonly asked questions like: is each gene a kinase?
  #Read in file, get gene name column, fix gene name, compare to list, set answer to 1/0, overwrite old file
  #Repeat this process for each gene symbol list defined
  #TODO: Get this out of Main and do it without having to know about this list of files to be annotated in a hash
  #TODO: - produce a separate .cancer.anno.tsv file with the annotation appended and the sum of count for each row as an extra column
  $step++; print MAGENTA, "\n\nStep $step. Annotating gene files", RESET;
  &annotateGeneFiles('-gene_symbol_lists'=>$gene_symbol_lists, '-out_paths'=>$out_paths, '-verbose'=>$verbose);

  #TODO: Replace this with use of DGIdb command line tool that performs queries against more sources with better cancer relevance filtering
  #Create drugDB interaction files
  #Perform druggable genes analysis on each list (filtered, kinase-only, inhibitor-only, antineoplastic-only)
  $step++; print MAGENTA, "\n\nStep $step. Intersecting gene lists with druggable genes of various categories", RESET;
  &drugDbIntersections('-script_dir'=>$script_dir, '-out_paths'=>$out_paths, '-verbose'=>$verbose);

  #For each of the following: WGS SNVs, Exome SNVs, and WGS+Exome SNVs, do the following:
  #Get BAM readcounts for WGS (tumor/normal), Exome (tumor/normal), RNAseq (tumor), RNAseq (normal) - as available of course
  $self->wgs_positions_file($out_paths->{'wgs'}->{'snv'}->{path});
  $self->exome_positions_file($out_paths->{'exome'}->{'snv'}->{path});
  $self->wgs_exome_positions_file($out_paths->{'wgs_exome'}->{'snv'}->{path});
  
  #$self->tumor_fpkm_file($out_paths->{'tumor_rnaseq_cufflinks_absolute'}->{'isoforms.merged.fpkm.expsort.tsv'}->{path});
  
  # TODO: switch to relaying the annotation build ID instead. It is currently only needed by SummarizeTier1SnvSupport
  # TODO: refactor that and then remove this as an output from main.  Then the subroutine that determines it can also be removed
  # TODO: actually, all SummarizeTier1SnvSupport does with it is pass it to GetBamReadCounts which no longer accepts it! 

  #Store outputs needed for summarize tier1 snv support
  for my $p (qw/wgs_positions_file exome_positions_file wgs_exome_positions_file tumor_fpkm_file/) {
        no warnings;
        $self->status_message("$p set to " . $self->$p . "\n");
  }

  #print Dumper $out_paths;
  print "\n\nPROCESSING COMPLETE\n\n";

  return(1);
}


###############################################################################################################################
#Determine the ensembl version used by examining the underlying input builds of the clinseq build                             #
###############################################################################################################################
sub getAnnotationBuild{
  my $self = shift;
  my %args = @_;
  my $clinseq_build_id = $args{'-clinseq_build_id'};

  my $clinseq_build = Genome::Model::Build->get($clinseq_build_id);
  my $clinseq_model_id = $clinseq_build->model->id;

  my ($wgs_somvar_build, $exome_somvar_build, $tumor_rnaseq_build, $normal_rnaseq_build, $wgs_normal_refalign_build, $wgs_tumor_refalign_build, $exome_normal_refalign_build, $exome_tumor_refalign_build);
  $wgs_somvar_build = $clinseq_build->wgs_build;
  $exome_somvar_build = $clinseq_build->exome_build;
  $tumor_rnaseq_build = $clinseq_build->tumor_rnaseq_build;
  $normal_rnaseq_build = $clinseq_build->normal_rnaseq_build;
  $wgs_normal_refalign_build = $wgs_somvar_build->normal_build if ($wgs_somvar_build);
  $wgs_tumor_refalign_build = $wgs_somvar_build->tumor_build if ($wgs_somvar_build);
  $exome_normal_refalign_build = $exome_somvar_build->normal_build if ($exome_somvar_build);
  $exome_tumor_refalign_build = $exome_somvar_build->tumor_build if ($exome_somvar_build);
  my @builds = ($wgs_normal_refalign_build, $wgs_tumor_refalign_build, $wgs_somvar_build, $exome_normal_refalign_build, $exome_tumor_refalign_build, $exome_somvar_build, $tumor_rnaseq_build, $normal_rnaseq_build);

  my %annotation_refs;
  for my $build (@builds){
    next unless $build;
    my $m = $build->model;
    my $model_name = $m->name;
    my $pp_id = $m->processing_profile_id;
    my $pp = Genome::ProcessingProfile->get($pp_id);
    my $pp_type = $pp->type_name;
    if ($m->can("annotation_reference_build")){
      my $ab = $m->annotation_reference_build;
      if ($ab){
        my $ab_name = $ab->name;
        $annotation_refs{$ab_name}=$ab->id;
      }else{
        print YELLOW, "\n\nUndefined annotation build name for model!\n$model_name\n\n", RESET;
      }
    }elsif ($m->can("annotation_build")){
      my $ab = $m->annotation_build;
      if ($ab){
        my $ab_name = $ab->name;
        $annotation_refs{$ab_name}=$ab->id;
      }else{
        print YELLOW, "\n\nUndefined annotation build name for model!\n$model_name\n\n", RESET;
      }
    }
  }   

  my $ar_count = keys %annotation_refs;
  my $annotation_build_id;
  if ($ar_count == 0){
    print RED, "\n\nUnable to determine a single annotation reference name from the input models!\n\n", RESET;
    die;
  }elsif($ar_count == 1){
    foreach my $ar_name (keys %annotation_refs){
      $annotation_build_id = $annotation_refs{$ar_name};
    }
  }else{
    print RED, "\n\nFound conflicting annotation versions being used in the input models of this ClinSeq model: $clinseq_model_id\n\n", RESET;
    die;
  }

  return ($annotation_build_id);
}


###############################################################################################################################
#Get build directories for the three datatypes                                                                                #
###############################################################################################################################
sub getDataDirsAndBuilds{
  my %args = @_; #Contains a hash of wgs/exome/rna model or build ids and names for these

  my %data_paths;
  my %builds;
  
  my %arg_dt = (
    -wgs_som_var_data_set => 'wgs',
    -exome_som_var_data_set => 'exome',
    -tumor_rna_seq_data_set => 'tumor_rnaseq',
    -normal_rna_seq_data_set => 'normal_rnaseq',
  );

  for my $arg_name (keys %arg_dt) {
    # arg name is one of those defined in the hash above - if that particular arg name was not used in that call, skip
    # arg value is the actual model/build number
    my $arg_value = $args{$arg_name};
    next unless $arg_value;

    # this is the key to use in the data_path hash
    my $dt = $arg_dt{$arg_name};
    
    # this is used in error messages
    my $type = $dt;
    $type =~ s/wgs/WGS/;
    $type =~ s/rna_seq/RNA-seq/g;
    $type =~ s/_/ /;

    # identify the build and model
    my $build = $arg_value; #Genome::Model::Build->get($arg_value);
    my $model;
    if ($build) {
        # yay: build directly specified
        $model = $build->model;
    }
    else {
        print RED, "\n\nA $type ID was specified, but no model or build with that ID could be found!\n\n", RESET;
        exit 1;
    }

    $builds{$dt} = $build;
  
    # the build directory root
    my $root = $data_paths{$dt}{root} = $build->data_directory . '/';
   
    # record paths to essential data based on data type ($dt)
    if ($dt =~ /rna/i) {
        # tumor rna and normal rna
        my $reference_build = $model->reference_sequence_build;
        $data_paths{$dt}{reference_fasta_path} = $reference_build->full_consensus_path('fa');

        my $alignment_result = $build->alignment_result;
        $data_paths{$dt}{bam} = $alignment_result->bam_file;
    
        $data_paths{$dt}{alignments} = $root."alignments/";
        $data_paths{$dt}{coverage} = $root."coverage/";
        $data_paths{$dt}{expression} = $root."expression/";
        $data_paths{$dt}{logs} = $root."logs/";
        $data_paths{$dt}{reports} = $root."reports/";
    }
    else {
        # wgs and exome
        my $reference_build = $build->reference_sequence_build;
        $data_paths{$dt}{reference_fasta_path} = $reference_build->full_consensus_path('fa');

        $data_paths{$dt}{tumor_bam} = $build->tumor_bam;
        $data_paths{$dt}{normal_bam} = $build->normal_bam;
        
        $data_paths{$dt}{effects} = $root."effects/";
        $data_paths{$dt}{logs} = $root."logs/";
        $data_paths{$dt}{loh} = $root."loh/";
        $data_paths{$dt}{novel} = $root."novel/";
        $data_paths{$dt}{reports} = $root."reports/";
        $data_paths{$dt}{variants} = $root."variants/";
    }
  }

  #print Dumper \%data_paths;
  return(\%data_paths, \%builds);
}


###################################################################################################################
#Summarize SNVs/Indels                                                                                            #
###################################################################################################################
sub importSNVs{
  my $self = shift;
  my %args = @_;
  my $data_paths = $args{'-data_paths'};
  my $out_paths = $args{'-out_paths'};
  my $patient_dir = $args{'-patient_dir'};
  my $entrez_ensembl_data = $args{'-entrez_ensembl_data'};
  my $filter_mt = $args{'-filter_mt'};
  my $verbose = $args{'-verbose'};

  #Create SNV/INDEL dirs for each data type. e.g.: 'snv', 'snv/wgs/', 'snv/exome/', 'snv/wgs_exome/'
  my $snv_dir = &createNewDir('-path'=>$patient_dir, '-new_dir_name'=>'snv', '-silent'=>1);
  my $indel_dir = &createNewDir('-path'=>$patient_dir, '-new_dir_name'=>'indel', '-silent'=>1);

  #Define variant effect type filters
  #TODO: Allow different filters to be used as a parameter
  my $snv_filter = "missense|nonsense|splice_site|splice_region|rna";
  my $indel_filter = "in_frame_del|in_frame_ins|frame_shift_del|frame_shift_ins|splice_site_ins|splice_site_del|rna";

  #Define the dataset: WGS SNV, WGS indel, Exome SNV, Exome indel
  my %dataset;
  if ($self->wgs_build){
    my $snv_wgs_dir = &createNewDir('-path'=>$snv_dir, '-new_dir_name'=>'wgs', '-silent'=>1);
    my $indel_wgs_dir = &createNewDir('-path'=>$indel_dir, '-new_dir_name'=>'wgs', '-silent'=>1);
    my $effects_dir = $data_paths->{'wgs'}->{'effects'};
    $dataset{'1'}{data_type} = "wgs";
    $dataset{'1'}{var_type} = "snv";
    $dataset{'1'}{effects_dir} = $effects_dir;
    $dataset{'1'}{t1_hq_annotated} = "snvs.hq.tier1.v1.annotated";
    $dataset{'1'}{t1_hq_annotated_top} = "snvs.hq.tier1.v1.annotated.top";
    $dataset{'1'}{compact_file} = "$snv_wgs_dir"."snvs.hq.tier1.v1.annotated.compact.tsv";
    $dataset{'1'}{aa_effect_filter} = $snv_filter;
    $dataset{'1'}{target_dir} = $snv_wgs_dir;

    $dataset{'2'}{data_type} = "wgs";
    $dataset{'2'}{var_type} = "indel";
    $dataset{'2'}{effects_dir} = $effects_dir;
    $dataset{'2'}{t1_hq_annotated} = "indels.hq.tier1.v1.annotated";
    $dataset{'2'}{t1_hq_annotated_top} = "indels.hq.tier1.v1.annotated.top";
    $dataset{'2'}{compact_file} = "$indel_wgs_dir"."indels.hq.tier1.v1.annotated.compact.tsv";
    $dataset{'2'}{aa_effect_filter} = $indel_filter;
    $dataset{'2'}{target_dir} = $indel_wgs_dir;
  }
  if ($self->exome_build){
    my $snv_exome_dir = &createNewDir('-path'=>$snv_dir, '-new_dir_name'=>'exome', '-silent'=>1);  
    my $indel_exome_dir = &createNewDir('-path'=>$indel_dir, '-new_dir_name'=>'exome', '-silent'=>1);
    my $effects_dir = $data_paths->{'exome'}->{'effects'};
    $dataset{'3'}{data_type} = "exome";
    $dataset{'3'}{var_type} = "snv";
    $dataset{'3'}{effects_dir} = $effects_dir;
    $dataset{'3'}{t1_hq_annotated} = "snvs.hq.tier1.v1.annotated";
    $dataset{'3'}{t1_hq_annotated_top} = "snvs.hq.tier1.v1.annotated.top";
    $dataset{'3'}{compact_file} = "$snv_exome_dir"."snvs.hq.tier1.v1.annotated.compact.tsv";
    $dataset{'3'}{aa_effect_filter} = $snv_filter;
    $dataset{'3'}{target_dir} = $snv_exome_dir;

    $dataset{'4'}{data_type} = "exome";
    $dataset{'4'}{var_type} = "indel";
    $dataset{'4'}{effects_dir} = $effects_dir;
    $dataset{'4'}{t1_hq_annotated} = "indels.hq.tier1.v1.annotated";
    $dataset{'4'}{t1_hq_annotated_top} = "indels.hq.tier1.v1.annotated.top";
    $dataset{'4'}{compact_file} = "$indel_exome_dir"."indels.hq.tier1.v1.annotated.compact.tsv";
    $dataset{'4'}{aa_effect_filter} = $indel_filter;
    $dataset{'4'}{target_dir} = $indel_exome_dir;
  }

  my %data_merge;
  foreach my $ds (sort {$a <=> $b} keys %dataset){
    my %data_out;

    #Make a copy of the high quality .annotated and .annotated.top files
    my $data_type = $dataset{$ds}{data_type};
    my $var_type = $dataset{$ds}{var_type};
    my $effects_dir = $dataset{$ds}{effects_dir};
    my $t1_hq_annotated = $dataset{$ds}{t1_hq_annotated};
    my $t1_hq_annotated_top = $dataset{$ds}{t1_hq_annotated_top};
    my $compact_file = $dataset{$ds}{compact_file};
    my $aa_effect_filter = $dataset{$ds}{aa_effect_filter};
    my $target_dir = $dataset{$ds}{target_dir};

    my $new_annotated_file = "$target_dir$t1_hq_annotated".".tsv";
    my $new_annotated_top_file = "$target_dir$t1_hq_annotated_top".".tsv";
    my $cp_cmd1 = "cp $effects_dir$t1_hq_annotated $new_annotated_file";
    my $cp_cmd2 = "cp $effects_dir$t1_hq_annotated_top $new_annotated_top_file";
    if ($verbose){print YELLOW, "\n\n$cp_cmd1", RESET;}
    Genome::Sys->shellcmd(cmd => $cp_cmd1);
    if ($verbose){print YELLOW, "\n\n$cp_cmd2", RESET;}
    Genome::Sys->shellcmd(cmd => $cp_cmd2);

    #Get a column count on the file and use this to determine the correct header for the annotated variant file
    my $col_count = 0;
    open (TMP, $new_annotated_file) || die "\n\nCould not open variant annotation file: $new_annotated_file\n\n";
    while(<TMP>){
      chomp($_);
      my @line = split("\t", $_);
      $col_count = scalar(@line);
    }
    close(TMP);

    my @input_headers; 
    if ($col_count == 24){
      #New header
      #chr start stop ref_base var_base var_type gene_name transcript_id species transcript_source transcript_version strand transcript_status var_effect_type coding_pos aa_change ucsc_cons domain all_domains deletion_substructures transcript_error default_gene_name gene_name_source ensembl_gene_id  (24)
      @input_headers = qw (chr start stop ref_base var_base var_type gene_name transcript_id species transcript_source transcript_version strand transcript_status var_effect_type coding_pos aa_change ucsc_cons domain all_domains deletion_substructures transcript_error default_gene_name gene_name_source ensembl_gene_id);
    }else{
      $self->error_message("Unexpected column count ($col_count) found in SNV/INDEL file - ClinSeq is not compatible with old annotation format...");
      die();
    }

    #Get AA changes from full .annotated file
    my %aa_changes;
    if ($verbose){print YELLOW, "\n\nReading: $new_annotated_file", RESET;}
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
      headers => \@input_headers,
      input => "$new_annotated_file",
      separator => "\t",
    );
    while (my $data = $reader->next) {
      my $coord = $data->{chr} .':'. $data->{start} .'-'. $data->{stop};
      $data->{coord} = $coord;
      #Apply the AA effect filter
      unless ($data->{var_effect_type} =~ /$aa_effect_filter/){
        next();
      }
      #Apply the MT/chrM filter
      if ($filter_mt){
        my $chr = $data->{chr};
        if ($chr =~ /^MT$|^chrMT$|^M$|^chrM$/i){
          next();
        }
      }
      $aa_changes{$coord}{$data->{aa_change}}=1;
    }

    #Get compact SNV info from the '.top' file but grab the complete list of AA changes from the '.annotated' file
    if ($verbose){print YELLOW, "\n\nReading: $new_annotated_top_file", RESET;}
    $reader = Genome::Utility::IO::SeparatedValueReader->create(
      headers => \@input_headers,
      input => "$new_annotated_top_file",
      separator => "\t",
    );

    while (my $data = $reader->next){
      my $coord = $data->{chr} .':'. $data->{start} .'-'. $data->{stop};
      #Apply the AA effect filter
      unless ($data->{var_effect_type} =~ /$aa_effect_filter/){
        next();
      }
      #Apply the MT/chrM filter
      if ($filter_mt){
        my $chr = $data->{chr};
        if ($chr =~ /^MT$|^chrMT$|^M$|^chrM$/i){
          next();
        }
      }
      my %aa = %{$aa_changes{$coord}};
      my $aa_string = join(",", sort keys %aa);
      $data_out{$coord}{gene_name} = $data->{gene_name};
      $data_merge{$var_type}{$coord}{gene_name} = $data->{gene_name};
      $data_out{$coord}{ensembl_gene_id} = $data->{ensembl_gene_id};
      $data_merge{$var_type}{$coord}{ensembl_gene_id} = $data->{ensembl_gene_id};
      $data_out{$coord}{aa_changes} = $aa_string;
      $data_merge{$var_type}{$coord}{aa_changes} = $aa_string;
      $data_out{$coord}{ref_base} = $data->{ref_base};
      $data_merge{$var_type}{$coord}{ref_base} = $data->{ref_base};
      $data_out{$coord}{var_base} = $data->{var_base};
      $data_merge{$var_type}{$coord}{var_base} = $data->{var_base};

      #Make note of the datatype (wgs or exome) this variant was called by...
      if ($data_type eq "wgs"){
        $data_merge{$var_type}{$coord}{wgs} = 1;
      }
      if ($data_type eq "exome"){
        $data_merge{$var_type}{$coord}{exome} = 1;
      }

      #Attempt to fix the gene name:
      my $fixed_gene_name = &fixGeneName('-gene'=>$data->{gene_name}, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
      $data_out{$coord}{mapped_gene_name} = $fixed_gene_name;
      $data_merge{$var_type}{$coord}{mapped_gene_name} = $fixed_gene_name;
    }

    #Print out the resulting list, sorting on fixed gene name
    open (OUT, ">$compact_file") || die "\n\nCould not open output file: $compact_file\n\n";
    print OUT "coord\tgene_name\tmapped_gene_name\tensembl_gene_id\taa_changes\tref_base\tvar_base\n";
    foreach my $coord (sort {$data_out{$a}->{mapped_gene_name} cmp $data_out{$b}->{mapped_gene_name}} keys %data_out){
      print OUT "$coord\t$data_out{$coord}{gene_name}\t$data_out{$coord}{mapped_gene_name}\t$data_out{$coord}{ensembl_gene_id}\t$data_out{$coord}{aa_changes}\t$data_out{$coord}{ref_base}\t$data_out{$coord}{var_base}\n";
    }
    close(OUT);

    #Store the path for this output file
    $out_paths->{$data_type}->{$var_type}->{path} = $compact_file;
  }

  #If both WGS and Exome data were present, print out a data merge for SNVs and Indels
  if ($self->wgs_build && $self->exome_build){
    my $snv_wgs_exome_dir = &createNewDir('-path'=>$snv_dir, '-new_dir_name'=>'wgs_exome', '-silent'=>1);
    my $indel_wgs_exome_dir = &createNewDir('-path'=>$indel_dir, '-new_dir_name'=>'wgs_exome', '-silent'=>1);
    my $snv_merge_file = "$snv_wgs_exome_dir"."snvs.hq.tier1.v1.annotated.compact.tsv";
    my $indel_merge_file = "$indel_wgs_exome_dir"."indels.hq.tier1.v1.annotated.compact.tsv";

    open (OUT, ">$snv_merge_file") || die "\n\nCould not open output file: $snv_merge_file\n\n";
    print OUT "coord\tgene_name\tmapped_gene_name\tensembl_gene_id\taa_changes\tref_base\tvar_base\twgs_called\texome_called\n";
    my %data_out = %{$data_merge{'snv'}};
    foreach my $coord (sort {$data_out{$a}->{mapped_gene_name} cmp $data_out{$b}->{mapped_gene_name}} keys %data_out){
      my $wgs_called = 0;
      if (defined($data_out{$coord}{wgs})){ $wgs_called = 1; }
      my $exome_called = 0;
      if (defined($data_out{$coord}{exome})){ $exome_called = 1; }
      print OUT "$coord\t$data_out{$coord}{gene_name}\t$data_out{$coord}{mapped_gene_name}\t$data_out{$coord}{ensembl_gene_id}\t$data_out{$coord}{aa_changes}\t$data_out{$coord}{ref_base}\t$data_out{$coord}{var_base}\t$wgs_called\t$exome_called\n";
    }
    close(OUT);
    $out_paths->{'wgs_exome'}->{'snv'}->{path} = $snv_merge_file;

    open (OUT, ">$indel_merge_file") || die "\n\nCould not open output file: $indel_merge_file\n\n";
    print OUT "coord\tgene_name\tmapped_gene_name\tensembl_gene_id\taa_changes\tref_base\tvar_base\twgs_called\texome_called\n";
    %data_out = %{$data_merge{'indel'}};
    foreach my $coord (sort {$data_out{$a}->{mapped_gene_name} cmp $data_out{$b}->{mapped_gene_name}} keys %data_out){
      my $wgs_called = 0;
      if (defined($data_out{$coord}{wgs})){ $wgs_called = 1; }
      my $exome_called = 0;
      if (defined($data_out{$coord}{exome})){ $exome_called = 1; }
      print OUT "$coord\t$data_out{$coord}{gene_name}\t$data_out{$coord}{mapped_gene_name}\t$data_out{$coord}{ensembl_gene_id}\t$data_out{$coord}{aa_changes}\t$data_out{$coord}{ref_base}\t$data_out{$coord}{var_base}\t$wgs_called\t$exome_called\n";
    }
    close(OUT);
    $out_paths->{'wgs_exome'}->{'indel'}->{path} = $indel_merge_file;
  }

  return();
}


###################################################################################################################################
#Annotate gene lists to deal with commonly asked questions like: is each gene a kinase?                                           #
#Read in file, get gene name column, fix gene name, compare to list, set answer to 1/0, overwrite old file                        #
#Repeat this process for each gene symbol list defined                                                                            #
###################################################################################################################################
sub annotateGeneFiles{
  my %args = @_;
  my $gene_symbol_lists = $args{'-gene_symbol_lists'};
  my $out_paths = $args{'-out_paths'};
  my $verbose = $args{'-verbose'};

  foreach my $type (keys %{$out_paths}){
    my $sub_types = $out_paths->{$type};
    foreach my $sub_type (keys %{$sub_types}){
      #Store the file input data for this file
      my $path = $sub_types->{$sub_type}->{'path'};
      if ($verbose){print "\n\tProcessing: $path";}
      my $new_path = $path.".tmp";
      open (INDATA, "$path") || die "\n\nCould not open input datafile: $path\n\n";
      my %data;
      my %cols;
      my $header = 1;
      my $header_line = '';
      my $l = 0;
      while(<INDATA>){
        $l++;
        chomp($_);
        my $record = $_;
        my @line = split("\t", $_);
        if ($header == 1){
          my $c = 0;
          $header_line = $_;
          foreach my $colname (@line){
            $cols{$colname}{position} = $c;
            $c++;
          }
          $header = 0;
          unless ($cols{'mapped_gene_name'}){
            print RED, "\n\nFile has no 'mapped_gene_name' column: $path\n\n", RESET;
            exit 1;
          }
          next();
        }
        $data{$l}{record} = $record;
        $data{$l}{gene_name} = $line[$cols{'mapped_gene_name'}{position}];
      }
      close(INDATA);

      #Figure out the gene matches to the gene symbol lists
      #Test each gene name in this column against those in the list and add a column with the match status (i.e. is is a kinase, cancer gene, etc.)
      foreach my $l (keys %data){
        my $gene_name = $data{$l}{gene_name};
        foreach my $gene_symbol_type (keys %{$gene_symbol_lists}){
          my $gene_symbols = $gene_symbol_lists->{$gene_symbol_type}->{symbols};
          if ($gene_symbols->{$gene_name}){
            $data{$l}{$gene_symbol_type} = 1;
          }else{
            $data{$l}{$gene_symbol_type} = 0;
          }
        }
      }
      #Print out a new file contain the extra columns
      open (OUTDATA, ">$new_path") || die "\n\nCould not open output datafile: $new_path\n\n";
      my @gene_symbol_list_names = sort {$gene_symbol_lists->{$a}->{order} <=> $gene_symbol_lists->{$b}->{order}} keys %{$gene_symbol_lists};
      my $gene_symbol_list_name_string = join("\t", @gene_symbol_list_names);
      print OUTDATA "$header_line\t$gene_symbol_list_name_string\n";
      foreach my $l (sort {$a <=> $b} keys %data){
        my @tmp;
        foreach my $gene_symbol_list_name (@gene_symbol_list_names){
          push (@tmp, $data{$l}{$gene_symbol_list_name});
        }
        my $new_cols_string = join("\t", @tmp);
        print OUTDATA "$data{$l}{record}\t$new_cols_string\n";
      }
      close(OUTDATA);

      #Replace the original file with the new file
      my $mv_cmd = "mv $new_path $path";
      if ($verbose){print YELLOW, "\n\t\t $mv_cmd", RESET;}
      Genome::Sys->shellcmd(cmd => $mv_cmd);
    }
  }
  return();
}


###################################################################################################################################
#Create drugDB interaction files                                                                                                  #
###################################################################################################################################
sub drugDbIntersections{
  my %args = @_;
  my $script_dir = $args{'-script_dir'};
  my $out_paths = $args{'-out_paths'};
  my $verbose = $args{'-verbose'};

  my $drugdb_script = "$script_dir"."summary/identifyDruggableGenes.pl";

  my $drugbank_interactions_dir = "/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/DrugBank/query_files/";
  my %filter_options;
  $filter_options{'3'}{name} = ".default";   #Default filter
  $filter_options{'4'}{name} = ".antineo";   #Anti-neoplastic only
  $filter_options{'5'}{name} = ".inhibitor"; #Inhibitor only
  $filter_options{'6'}{name} = ".kinase";    #Kinases only

  foreach my $type (keys %{$out_paths}){
    my $sub_types = $out_paths->{$type};
    foreach my $sub_type (keys %{$sub_types}){
      #Store the file input data for this file
      my $path = $sub_types->{$sub_type}->{'path'};
      my $name_col = &getColumnPosition('-path'=>$path, '-column_name'=>'mapped_gene_name');

      #Note that this function returns the 0-based column position - The script below assumes 1 based
      $name_col += 1;

      #Get file path with the file extension removed:
      my $fb = &getFilePathBase('-path'=>$path);
      
      my $dgidb_dir = $fb->{$path}->{base_dir} . "dgidb/";
      my $drugbank_dir = $dgidb_dir . "drugbank/";
      
      unless (-e $dgidb_dir && -d $dgidb_dir){
        mkdir ($dgidb_dir);
      }
      unless (-e $drugbank_dir && -d $drugbank_dir){
        mkdir ($drugbank_dir);
      }

      #Run with each filtering option
      foreach my $filter (sort {$a <=> $b} keys %filter_options){
        my $filter_name = $filter_options{$filter}{name};
        my $out = $drugbank_dir . $fb->{$path}->{file_base} . "$filter_name" . $fb->{$path}->{extension};

        if (-e $out){
          if ($verbose){print YELLOW, "\n\tFile already exists - skipping ($out)", RESET;} 
        }else{
          my $cmd = "$drugdb_script --candidates_file=$path  --name_col_1=$name_col  --interactions_file=$drugbank_interactions_dir/DrugBank_WashU_INTERACTIONS.filtered."."$filter".".tsv  --name_col_2=12 > $out";
          if ($verbose){print YELLOW, "\n\t$cmd", RESET;}
          Genome::Sys->shellcmd(cmd => "$cmd");
        }
      }
    }
  }

  my $santa_monica_interactions_dir = "/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/SantaMonicaLung/";

  foreach my $type (keys %{$out_paths}){
    my $sub_types = $out_paths->{$type};
    foreach my $sub_type (keys %{$sub_types}){
      #Store the file input data for this file
      my $path = $sub_types->{$sub_type}->{'path'};
      my $name_col = &getColumnPosition('-path'=>$path, '-column_name'=>'mapped_gene_name');

      #Note that this function returns the 0-based column position - The script below assumes 1 based
      $name_col += 1;

      #Get file path with the file extension removed:
      my $fb = &getFilePathBase('-path'=>$path);
      
      my $dgidb_dir = $fb->{$path}->{base_dir} . "dgidb/";
      my $santa_monica_dir = $dgidb_dir . "santa_monica_lung/";
      
      unless (-e $dgidb_dir && -d $dgidb_dir){
        mkdir ($dgidb_dir);
      }
      unless (-e $santa_monica_dir && -d $santa_monica_dir){
        mkdir ($santa_monica_dir);
      }
      my $filter_name = ".default";
      my $out = $santa_monica_dir . $fb->{$path}->{file_base} . "$filter_name" . $fb->{$path}->{extension};
      if (-e $out){
        if ($verbose){print YELLOW, "\n\tFile already exists - skipping ($out)", RESET;} 
      }else{
        my $cmd = "$drugdb_script --candidates_file=$path  --name_col_1=$name_col  --interactions_file=$santa_monica_interactions_dir"."SantaMonicaLungCancerDrugDatabase.tsv  --name_col_2=1 > $out";
        if ($verbose){print YELLOW, "\n\t$cmd", RESET;}
        Genome::Sys->shellcmd(cmd => "$cmd");
      }
    }
  }

  return();
}



1;
