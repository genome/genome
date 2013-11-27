package Genome::Model::ClinSeq::Command::RunCnView;
use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::RunCnView {
    is => 'Command::V2',
    has_input => [
        build  => { 
            is => 'Genome::Model::Build::SomaticVariation',
            is_many => 0,
            shell_args_position => 1,
            require_user_verify => 0,
            doc => 'somatic variation model to be used for CnView runs'
        },
        cancer_annotation_db => {
            is => 'Genome::Db',
            example_values => ['tgi/cancer-annotation/human/build37-20130401.1'],
            doc => 'data set of cancer annotation to use for analysis',
        },
        outdir => { 
            is => 'FilesystemPath',
            doc => 'Directory where output files will be written', 
        },
        cnv_hmm_file => {
            is => 'FilesystemPath',
            doc => 'CNVhmm results file from clonality analysis (or elsewhere)'
        },
        test => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Use for shorter tests',
        }
    ],
    has_output => [
        gene_amp_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
        gene_del_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
        gene_ampdel_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
    ],
    doc => 'run gmt copy-number cn-view in various ways and summarize results',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq run-cn-view --outdir=/tmp/cnv/ --cnv-hmm-file=/gscmnt/gc1401/info/model_data/2889933976/build132760359/AML109/clonality/cnaseq.cnvhmm 129399487 

EOS
}

sub help_detail {
    return <<EOS
Run gmt copy-number cn-view in various ways and summarize results

(put more content here)
EOS
}

sub execute {
  my $self = shift;
  my $build = $self->build;
  my $outdir = $self->outdir;
  my $cnv_hmm_file = $self->cnv_hmm_file;
  $outdir .= "/" unless ($outdir =~ /\/$/);
$DB::single = 1;
  my @cnv_symbol_lists;
  if ($self->test){
    @cnv_symbol_lists = qw (Kinase_dGene);
  }else{
    @cnv_symbol_lists = qw (All Kinase_dGene CancerGeneCensusPlus_Sanger AntineoplasticTargets_DrugBank);
  }

  my $cancer_annotation_db = $self->cancer_annotation_db;

  my $gene_symbol_lists_dir = $cancer_annotation_db->data_directory . "/GeneSymbolLists/";
  #my $gene_symbol_lists_dir = "/gscmnt/sata132/techd/mgriffit/reference_annotations/GeneSymbolLists";

  die $self->error_message("outdir does not exist") unless (-e $outdir && -d $outdir);
  die $self->error_message("segments file does not exist: $cnv_hmm_file") unless (-e $cnv_hmm_file);
  die $self->error_message("could not resolve annotation build id") unless ($build->annotation_build->id);
  die $self->error_message("could not find gene symbol lists dir: $gene_symbol_lists_dir") unless (-e $gene_symbol_lists_dir);
  my $annotation_build_id = $build->annotation_build->id;
 
  my $variants_dir = $build->data_directory . "/variants/";
  my $cnv_data_file = $variants_dir . "cnvs.hq";
  die $self->error_message("could not find cnvs.hq file somatic build here: $variants_dir") unless (-e $cnv_data_file);

  #Create main CNV dir: 'cnv'
  my $cnview_dir = $outdir . "cnview/";

  #Create a copy of the cnvs.hq file for later convenience
  my $new_cnv_data_file = $outdir . "cnvs.hq";
  Genome::Sys->copy_file($cnv_data_file, $new_cnv_data_file) unless (-e $new_cnv_data_file);

  #For each list of gene symbols, run the CNView analysis
  foreach my $symbol_list_name (@cnv_symbol_lists){
    if ($symbol_list_name eq "All"){
      my $cnview_cmd = Genome::Model::Tools::CopyNumber::CnView->create(annotation_build=>$build, cnv_file=>$cnv_data_file, segments_file=>$cnv_hmm_file, output_dir=>$cnview_dir, name=>$symbol_list_name, cancer_annotation_db => $cancer_annotation_db);
      $cnview_cmd->execute();

      #Copy these files to the top CNV dir
      my $new_dir = "$cnview_dir"."CNView_"."$symbol_list_name/";

      my @suffixes = qw (_genes _genes.amp _genes.del _genes.ampdel _transcripts _transcripts.amp _transcripts.del _transcripts.ampdel);
      foreach my $suffix (@suffixes){
        my $path1 = "$new_dir"."CNView_"."$symbol_list_name"."$suffix".".tsv";
        my $path2 = "$cnview_dir"."cnv."."$symbol_list_name"."$suffix".".tsv";
        Genome::Sys->copy_file($path1, $path2);
        my $dataname = "cnv".$suffix;
      }
    }else{
      my $gene_targets_file = "$gene_symbol_lists_dir/$symbol_list_name".".txt";
      my $cnview_cmd = Genome::Model::Tools::CopyNumber::CnView->create(annotation_build=>$build, cnv_file=>$cnv_data_file, segments_file=>$cnv_hmm_file, output_dir=>$cnview_dir, gene_targets_file=>$gene_targets_file, name=>$symbol_list_name, cancer_annotation_db => $cancer_annotation_db);
      $cnview_cmd->execute();
    }
  }

  #Set the output values that will be needed for downstream steps
  $self->gene_amp_file($cnview_dir . "cnv.All_genes.amp.tsv");
  $self->gene_del_file($cnview_dir . "cnv.All_genes.del.tsv");
  $self->gene_ampdel_file($cnview_dir . "cnv.All_genes.ampdel.tsv");

  return 1;
}

1;

