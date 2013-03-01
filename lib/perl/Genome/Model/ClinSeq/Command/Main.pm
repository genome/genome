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

  #Set flags for each datatype
  my $wgs = exists $builds->{wgs} || 0;
  my $exome = exists $builds->{exome} || 0;
  my $tumor_rnaseq = exists $builds->{tumor_rnaseq} || 0;
  my $normal_rnaseq = exists $builds->{normal_rnaseq} || 0;

  #Check the working dir
  $working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");

  #Create a hash for storing output files as they are created
  my %out_paths;
  my $out_paths = \%out_paths;

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
