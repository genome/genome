package Genome::Model::ClinSeq::Command::CufflinksDifferentialExpression;
use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::RnaSeqUtil qw(:all);

class Genome::Model::ClinSeq::Command::CufflinksDifferentialExpression {
    is => 'Command::V2',
    has_input => [
        case_build => {
            is => 'Genome::Model::Build::RnaSeq',
            doc => 'RnaSeq build of interest. This should be the case sample (e.g. tumor, metastasis, drug resistant, etc.)',
        },
        control_build => { 
            is => 'Genome::Model::Build::RnaSeq',
            doc => 'RnaSeq build to compare against. This should be the control sample (i.e. normal, primary, drug sensitive, etc.)',
        },
        cancer_annotation_db => {
            is => 'Genome::Db::Tgi::CancerAnnotation',
            doc => 'cancer annotation db',
        },
        outdir => { 
            is => 'FilesystemPath',
            doc => 'Directory where output files will be written', 
        },
    ],
    has_output => [
        coding_hq_up_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
        coding_hq_down_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
        coding_hq_de_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
    ],
    doc => 'perform simple differential expression comparison between two samples using FPKM values from Cufflinks',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq cufflinks-differential-expression --outdir=/tmp/ --case-build=129767889 --control-build=129767952

EOS
}

sub help_detail {
    return <<EOS

Perform naive pairwise differential expression between two samples (e.g. tumor vs. normal)

Report a complete list of differential expression values for all genes in a single file

Also create over-expressed and under-expressed only files

Differential expression is based on gene/transcript expression estimates derived from the isoforms FPKM file of Cufflinks

EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__;

  #Check output dir
  unless(-e $self->outdir && -d $self->outdir) {
    push @errors,UR::Object::Tag->create(
      type => 'error',
               properties => ['outdir'],
               desc => 'outdir does not exist or is not a directory',
      );
  }

  #Make sure both RNA-seq build are using the same annotation
  unless ($self->case_build->model->annotation_build->id == $self->control_build->model->annotation_build->id){
    push @errors,UR::Object::Tag->create(
      type => 'error',
               properties => ['case_build','control_build'],
               desc => 'case and control build were not generated with the same reference annotation build!',
      );
  }

  #Make sure both RNA-seq build are using the same reference build
  unless ($self->case_build->reference_sequence_build->id == $self->control_build->reference_sequence_build->id){
    push @errors,UR::Object::Tag->create(
      type => 'error',
               properties => ['case_build','control_build'],
               desc => 'case and control build were not generated with the same reference sequence build!',
      );
  }

  return @errors;
}

sub execute {
  my $self = shift;
  my $case_build = $self->case_build;
  my $control_build = $self->control_build;
  my $cancer_annotation_db = $self->cancer_annotation_db;
  my $outdir = $self->outdir;

  #Set some human readable case and control labels
  my $case_label = "case";
  my $control_label = "control";
  $case_label = $case_build->subject->common_name if ($case_build->subject->common_name);
  $case_label =~ s/ /\_/g;
  $control_label = $control_build->subject->common_name if ($control_build->subject->common_name);
  $control_label =~ s/ /\_/g;
  if ($case_label eq $control_label){
    $case_label = "case";
    $control_label = "control";
  }

  #Set up directories for output
  $outdir .= "/" unless ($outdir =~ /\/$/);
  my $genes_outdir = $outdir . "genes/";
  mkdir($genes_outdir);
  my $transcripts_outdir = $outdir . "transcripts/";
  mkdir($transcripts_outdir);

  #Find the cufflinks fpkm files for both builds
  my $case_build_dir = $case_build->data_directory;
  my $case_fpkm_file = $case_build_dir . "/expression/isoforms.fpkm_tracking";
  my $case_status_file = $case_build_dir . "/expression/genes.fpkm_tracking";

  my $control_build_dir = $control_build->data_directory;
  my $control_fpkm_file = $control_build_dir . "/expression/isoforms.fpkm_tracking";
  my $control_status_file = $control_build_dir . "/expression/genes.fpkm_tracking";

  unless (-e $case_fpkm_file && -e $control_fpkm_file) {
    $self->error_message("Could not find neccesary case/control fpkm files:\n$case_fpkm_file\n$control_fpkm_file\n");
    die $self->error_message;
  }

  #Build a map of ensembl transcript ids to gene ids and gene names from the gene annotation object associated with the rna-seq builds
  my $reference_build = $self->case_build->reference_sequence_build;
  my $reference_build_id = $reference_build->id;
  my $reference_build_name = $reference_build->name;
  $self->status_message("Processing RNA-seq data that was aligned to: $reference_build_name");

  my $annotation_build = $self->case_build->model->annotation_build;
  my $annotation_build_name = $annotation_build->name;
  my $annotation_data_dir = $annotation_build->data_directory;
  my $transcript_info_path = $annotation_data_dir . "/annotation_data/rna_annotation/$reference_build_id-transcript_info.tsv";
  my $gtf_path = $annotation_build->annotation_file('gtf',$reference_build_id);
  $self->status_message("Getting transcript to gene and gene name mappings from annotation build: $annotation_build_name");
  unless (defined($gtf_path)) {
    $self->error_message("'There is no annotation GTF file defined for annotation_reference_transcripts build: ". $annotation_build->__display_name__);
    die $self->error_message;
  }
  unless (-e $transcript_info_path) {
    $self->error_message("'There is no transcript info file for annotation_reference_transcripts build: ". $annotation_build->__display_name__);
    die $self->error_message;
  }
  $self->status_message("\t$transcript_info_path");

  my $ensembl_map = &loadEnsemblMap('-gtf_path'=>$gtf_path, '-transcript_info_path'=>$transcript_info_path);

  #Get Entrez and Ensembl data for gene name mappings
  $self->status_message("Load entrez and ensembl gene data");
  my $entrez_ensembl_data = &loadEntrezEnsemblData(-cancer_db => $cancer_annotation_db);

  #Parse the isoform fpkm files, create cleaner transcript level versions of these files and store them in the output dir
  $self->status_message("Parse transcript FPKM files from Cufflinks");
  my $fpkm;
  my $case_isoforms_file_sorted = "$transcripts_outdir"."case.transcripts.fpkm.namesort.tsv";
  $fpkm = &parseFpkmFile('-infile'=>$case_fpkm_file, '-outfile'=>$case_isoforms_file_sorted, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-ensembl_map'=>$ensembl_map, '-verbose'=>0);
  my $control_isoforms_file_sorted = "$transcripts_outdir"."control.transcripts.fpkm.namesort.tsv";
  $fpkm = &parseFpkmFile('-infile'=>$control_fpkm_file, '-outfile'=>$control_isoforms_file_sorted, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-ensembl_map'=>$ensembl_map, '-verbose'=>0);

  #Create gene-level files where the transcript level values are merged
  $self->status_message("Assemble gene estimates by merging Cufflinks transcript results");
  my $case_isoforms_merged_file_sorted = "$genes_outdir"."case.genes.fpkm.namesort.tsv";
  $fpkm = &mergeIsoformsFile('-infile'=>$case_fpkm_file, '-status_file'=>$case_status_file, '-outfile'=>$case_isoforms_merged_file_sorted, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-ensembl_map'=>$ensembl_map, '-verbose'=>0);
  my $control_isoforms_merged_file_sorted = "$genes_outdir"."control.genes.fpkm.namesort.tsv";
  $fpkm = &mergeIsoformsFile('-infile'=>$control_fpkm_file, '-status_file'=>$control_status_file, '-outfile'=>$control_isoforms_merged_file_sorted, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-ensembl_map'=>$ensembl_map, '-verbose'=>0);

  #Create a file containing basic gene ids, names, etc. along with FPKM values for case and control
  $self->status_message("Create transcript DE file");
  $fpkm = ();
  my $transcript_de_file = $self->create_de('-outdir'=>$transcripts_outdir, '-type'=>'transcript', '-case_file'=>$case_isoforms_file_sorted, '-control_file'=>$control_isoforms_file_sorted, '-fpkm'=>$fpkm);
  $self->status_message("Create gene DE file");
  $fpkm = ();
  my $gene_de_file = $self->create_de('-outdir'=>$genes_outdir, '-type'=>'gene', '-case_file'=>$case_isoforms_merged_file_sorted, '-control_file'=>$control_isoforms_merged_file_sorted, '-fpkm'=>$fpkm);
 
  #Determine path to R script to process the DE files
  my $r_de_script = __FILE__ . '.R';
  unless (-e $r_de_script){
    $self->error_message("Could not find companion R script");
    die $self->error_message;
  }

  #Feed this file into an R script that performs the actual differential expression analysis:

  #genes
  $self->status_message("\nPerforming gene-level DE analysis in R"); 
  my $gene_r_stderr = $genes_outdir . "CufflinksDifferentialExpression.pm.R.stderr";
  my $gene_r_stdout = $genes_outdir . "CufflinksDifferentialExpression.pm.R.stdout";
  my $r_cmd_gene = "$r_de_script $genes_outdir $gene_de_file 'gene' '$case_label' '$control_label' 2>$gene_r_stderr 1>$gene_r_stdout";
  $self->status_message($r_cmd_gene);
  Genome::Sys->shellcmd(cmd => $r_cmd_gene);

  #transcripts
  $self->status_message("\nPerforming transcript-level DE analysis in R"); 
  my $transcript_r_stderr = $transcripts_outdir . "CufflinksDifferentialExpression.pm.R.stderr";
  my $transcript_r_stdout = $transcripts_outdir . "CufflinksDifferentialExpression.pm.R.stdout";
  my $r_cmd_transcript = "$r_de_script $transcripts_outdir $transcript_de_file 'transcript' '$case_label' '$control_label' 2>$transcript_r_stderr 1>$transcript_r_stdout";
  $self->status_message($r_cmd_transcript);
  Genome::Sys->shellcmd(cmd => $r_cmd_transcript);

  #Set outputs
  my $coding_hq_up_file = $genes_outdir . "case_vs_control.coding.hq.up.tsv";
  die $self->error_message("Trying to set a file as output but the file does not exist: $coding_hq_up_file") unless (-e $coding_hq_up_file);
  $self->coding_hq_up_file($coding_hq_up_file);

  my $coding_hq_down_file = $genes_outdir . "case_vs_control.coding.hq.down.tsv";
  die $self->error_message("Trying to set a file as output but the file does not exist: $coding_hq_down_file") unless (-e $coding_hq_down_file);
  $self->coding_hq_down_file($coding_hq_down_file);

  my $coding_hq_de_file = $genes_outdir . "case_vs_control.coding.hq.de.tsv";
  die $self->error_message("Trying to set a file as output but the file does not exist: $coding_hq_de_file") unless (-e $coding_hq_de_file);
  $self->coding_hq_de_file($coding_hq_de_file);

  return 1;
}


sub create_de{
  my $self = shift;
  my %args = @_;
  my $type = $args{'-type'};
  my $case_file = $args{'-case_file'};
  my $control_file = $args{'-control_file'};
  my $fpkm = $args{'-fpkm'};
  my $outdir = $args{'-outdir'};

  my $header = 1;
  my $header_line;
  open (CASE, "$case_file") || die "\n\nCould not open case file\n\n";
  my %columns1;
  while(<CASE>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      my $p = 0;
      foreach my $column (@line){
        $columns1{$column}{pos} = $p;
        $p++;
      }
      $header_line = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]";
      $header = 0;
      next;
    }
    my $id = $line[0];
    $fpkm->{$id}->{info} = "$line[1]\t$line[2]\t$line[3]\t$line[4]";
    $fpkm->{$id}->{case_fpkm} = $line[$columns1{'FPKM'}{pos}];
    $fpkm->{$id}->{case_fpkm_conf_lo} = $line[$columns1{'FPKM_conf_lo'}{pos}];
    $fpkm->{$id}->{case_fpkm_conf_hi} = $line[$columns1{'FPKM_conf_hi'}{pos}];
    $fpkm->{$id}->{case_fpkm_status} = $line[$columns1{'FPKM_status'}{pos}];
  }
  close (CASE);

  $header = 1;
  open (CTRL, "$control_file") || die "\n\nCould not open control file\n\n";
  my %columns2;
  while(<CTRL>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      my $p = 0;
      foreach my $column (@line){
        $columns2{$column}{pos} = $p;
        $p++;
      }
      $header = 0;
      next;
    }
    my $id = $line[0];
    $fpkm->{$id}->{info} = "$line[1]\t$line[2]\t$line[3]\t$line[4]";
    $fpkm->{$id}->{control_fpkm} = $line[$columns2{'FPKM'}{pos}];
    $fpkm->{$id}->{control_fpkm_conf_lo} = $line[$columns2{'FPKM_conf_lo'}{pos}];
    $fpkm->{$id}->{control_fpkm_conf_hi} = $line[$columns2{'FPKM_conf_hi'}{pos}];
    $fpkm->{$id}->{control_fpkm_status} = $line[$columns2{'FPKM_status'}{pos}];
  }
  close (CASE);

  my $de_file = $outdir."$type".".de.input.tsv";

  open (DE, ">$de_file") || die "\n\nCould not open outfile: $de_file\n\n";
  print DE "$header_line\tcase_fpkm\tcase_fpkm_conf_hi\tcase_fpkm_conf_lo\tcase_fpkm_status\tcontrol_fpkm\tcontrol_fpkm_conf_hi\tcontrol_fpkm_conf_lo\tcontrol_fpkm_status\n";
  foreach my $id (sort keys %{$fpkm}){
    my $info = $fpkm->{$id}->{info};
    my ($case_fpkm, $case_fpkm_conf_hi, $case_fpkm_conf_lo, $case_fpkm_status);
    my ($control_fpkm, $control_fpkm_conf_hi, $control_fpkm_conf_lo, $control_fpkm_status);
    if (defined($fpkm->{$id}->{case_fpkm})){
      $case_fpkm = $fpkm->{$id}->{case_fpkm};
      $case_fpkm_conf_hi = $fpkm->{$id}->{case_fpkm_conf_hi};
      $case_fpkm_conf_lo = $fpkm->{$id}->{case_fpkm_conf_lo};
      $case_fpkm_status = $fpkm->{$id}->{case_fpkm_status};
    }else{
      $case_fpkm = "na";
      $case_fpkm_conf_hi = "na";
      $case_fpkm_conf_lo = "na";
      $case_fpkm_status = "na";
    }
    if (defined($fpkm->{$id}->{control_fpkm})){
      $control_fpkm = $fpkm->{$id}->{control_fpkm};
      $control_fpkm_conf_hi = $fpkm->{$id}->{control_fpkm_conf_hi};
      $control_fpkm_conf_lo = $fpkm->{$id}->{control_fpkm_conf_lo};
      $control_fpkm_status = $fpkm->{$id}->{control_fpkm_status};
    }else{
      $control_fpkm = "na";
      $control_fpkm_conf_hi = "na";
      $control_fpkm_conf_lo = "na";
      $control_fpkm_status = "na";
    }
    print DE "$id\t$info\t$case_fpkm\t$case_fpkm_conf_hi\t$case_fpkm_conf_lo\t$case_fpkm_status\t$control_fpkm\t$control_fpkm_conf_hi\t$control_fpkm_conf_lo\t$control_fpkm_status\n";
  }
  close (DE);

  return($de_file);
}


1;
