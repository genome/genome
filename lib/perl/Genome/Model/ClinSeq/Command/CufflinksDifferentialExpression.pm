package Genome::Model::ClinSeq::Command::CufflinksDifferentialExpression;
use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::RnaSeqUtil qw(:all);

use Data::Dumper;

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
        outdir => { 
            is => 'FilesystemPath',
            doc => 'Directory where output files will be written', 
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
  my $outdir = $self->outdir;
    
  #Find the cufflinks fpkm files for both builds
  my $case_build_dir = $case_build->data_directory;
  my $case_fpkm_file = $case_build_dir . "/expression/isoforms.fpkm_tracking";
  my $control_build_dir = $control_build->data_directory;
  my $control_fpkm_file = $control_build_dir . "/expression/isoforms.fpkm_tracking";
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
  $self->status_message("Getting transcript to gene and gene name mappings from annotation build: $annotation_build_name");

  my $gtf_path = $annotation_build->annotation_file('gtf',$reference_build_id);
  unless (defined($gtf_path)) {
    $self->error_message("'There is no annotation GTF file defined for annotation_reference_transcripts build: ". $annotation_build->__display_name__);
    die $self->error_message;
  }

  my %ensembl_map;
  open (GTF, "$gtf_path") || die "\n\nCould not open GTF file: $gtf_path";
  while(<GTF>){
    chomp($_);
    my @line = split("\t", $_);
    my @anno = split(";", $line[8]);
    my ($gene_name, $gene_id, $transcript_id);
    if ($anno[0] =~ /gene_name\s+\"(.*)\"/){
      $gene_name = $1;
    }
    if ($anno[1] =~ /gene_id\s+\"(.*)\"/){
      $gene_id = $1;
    }
    if ($anno[2] =~ /transcript_id\s+\"(.*)\"/){
      $transcript_id = $1;
    }
    unless ($gene_name && $gene_id && $transcript_id){
      $self->error_message("Could not parse gene_name, gene_id, transcript_id from GTF in line:\n$_\n");
      die $self->error_message;
    }
    $ensembl_map{$transcript_id}{ensg_id} = $gene_id;
    $ensembl_map{$transcript_id}{ensg_name} = $gene_name;
  }
  close(GTF);

  #Get Entrez and Ensembl data for gene name mappings
  my $entrez_ensembl_data = &loadEntrezEnsemblData();

  #Parse the isoform fpkm files, create cleaner transcript level versions of these files and store them in the output dir
  my $fpkm;
  my $case_isoforms_file_sorted = "$outdir"."case.transcripts.fpkm.namesort.tsv";
  $fpkm = &parseFpkmFile('-infile'=>$case_fpkm_file, '-outfile'=>$case_isoforms_file_sorted, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
  my $control_isoforms_file_sorted = "$outdir"."control.transcripts.fpkm.namesort.tsv";
  $fpkm = &parseFpkmFile('-infile'=>$control_fpkm_file, '-outfile'=>$control_isoforms_file_sorted, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);

  #Create gene-level files where the transcript level values are merged
  my $case_isoforms_merged_file_sorted = "$outdir"."case.genes.fpkm.namesort.tsv";
  $fpkm = &mergeIsoformsFile('-infile'=>$case_fpkm_file, '-outfile'=>$case_isoforms_merged_file_sorted, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-ensembl_map'=>\%ensembl_map, '-verbose'=>0);
  my $control_isoforms_merged_file_sorted = "$outdir"."control.genes.fpkm.namesort.tsv";
  $fpkm = &mergeIsoformsFile('-infile'=>$control_fpkm_file, '-outfile'=>$control_isoforms_merged_file_sorted, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-ensembl_map'=>\%ensembl_map, '-verbose'=>0);

  #Create a file containing basic gene ids, names, etc. along with FPKM values for case and control
  $fpkm = ();
  $self->create_de('-type'=>'transcript', '-case_file'=>$case_isoforms_file_sorted, '-control_file'=>$control_isoforms_file_sorted, '-fpkm'=>$fpkm);
  $fpkm = ();
  $self->create_de('-type'=>'gene', '-case_file'=>$case_isoforms_merged_file_sorted, '-control_file'=>$control_isoforms_merged_file_sorted, '-fpkm'=>$fpkm);
 
  #Determine path to R script to process the DE files
  my $r_de_script = __FILE__ . '.R';
  unless (-e $r_de_script){
    $self->error_message("Could not find companion R script");
    die $self->error_message;
  }

  #Feed this file into an R script that performs the actual differential expression analysis:
  #- produce de file for all genes with fold change and de status (de hq, de lq, no change)
  #- determine fold change as log2 difference?
  #- hq will be those diffs where the conf intervals do not overlap and if available, both have status of 'OK'
  #- produce .hq and .lq differential expression files for both gains, losses, and gains+losses
  #- 'de' will be those that exceed fold-change cutoff
  #- files:
  #- case_vs_control.tsv
  #- case_vs_control.hq.de.tsv, case_vs_control.hq.up.tsv, case_vs_control.hq.down.tsv
  #- case_vs_control.lq.de.tsv, case_vs_control.lq.up.tsv, case_vs_control.lq.down.tsv
  #- Create a plot showing the distribution of FPKM values for both samples
  #- Create a plot of FPKM case vs. FPKM control.  Color plot with hq de genes/transcripts



  #Perform basic some checking on the results files


  return 1;
}

sub create_de{
  my $self = shift;
  my %args = @_;
  my $type = $args{'-type'};
  my $case_file = $args{'-case_file'};
  my $control_file = $args{'-control_file'};
  my $fpkm = $args{'-fpkm'};

  my $header = 1;
  my $header_line;
  open (CASE, "$case_file") || die "\n\nCould not open case file\n\n";
  while(<CASE>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      $header_line = "$line[0]\t$line[1]\t$line[2]\t$line[3]";
      $header = 0;
      next;
    }
    my $id = $line[0];
    $fpkm->{$id}->{info} = "$line[1]\t$line[2]\t$line[3]";
    $fpkm->{$id}->{case_fpkm} = $line[6];
    $fpkm->{$id}->{case_fpkm_conf_lo} = $line[7];
    $fpkm->{$id}->{case_fpkm_conf_hi} = $line[8];
    $fpkm->{$id}->{case_fpkm_status} = $line[9];
  }
  close (CASE);

  $header = 1;
  open (CTRL, "$control_file") || die "\n\nCould not open control file\n\n";
  while(<CTRL>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      $header = 0;
      next;
    }
    my $id = $line[0];
    $fpkm->{$id}->{info} = "$line[1]\t$line[2]\t$line[3]";
    $fpkm->{$id}->{control_fpkm} = $line[6];
    $fpkm->{$id}->{control_fpkm_conf_lo} = $line[7];
    $fpkm->{$id}->{control_fpkm_conf_hi} = $line[8];
    $fpkm->{$id}->{control_fpkm_status} = $line[9];
  }
  close (CASE);

  my $de_file = $self->outdir."$type".".de.input.tsv";

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

  return;
}


1;
