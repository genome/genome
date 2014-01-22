package Genome::Model::ClinSeq::Command::CufflinksExpressionAbsolute;
#Written by Malachi Griffith

use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::RnaSeqUtil qw(:all);

class Genome::Model::ClinSeq::Command::CufflinksExpressionAbsolute {
    is => 'Command::V2',
    has_input => [
        build => {
            is => 'Genome::Model::Build::RnaSeq',
            shell_args_position => 1,
            doc => 'RnaSeq build to analyze',
        },
        cancer_annotation_db => {
            is => 'Genome::Db',
        },
        outdir => { 
            is => 'FilesystemPath',
            doc => 'Directory where output files will be written', 
        },
        percent_cutoff => {
            is => 'Number',
            is_optional => 1,
            default => 1,
            doc => 'The top N% of genes will be printed to a filtered file based on this cutoff (e.g. 1 for top 1%)',

        },
        verbose => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Flag for verbose output',
        },
    ],
    has_output => [
        tumor_fpkm_file => {
          is => 'FilesystemPath',
          is_optional => 1,
        },
        tumor_fpkm_topnpercent_file => {
          is => 'FilesystemPath',
          is_optional =>1,
        },
    ],
    doc => 'perform simple differential expression comparison between two samples using FPKM values from Cufflinks',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq outlier-genes-absolute --outdir=/tmp/ --percent-cutoff=1 129396808

EOS
}

sub help_detail {
    return <<EOS

This script takes a Cufflinks build dir as input and processes the genes.fpkm_tracking AND isoforms.fpkm_tracking files into more useful forms

Create formatted files:

1.) For genes.fpkm_tracking. Key on tracking_id AND locus, fix the gene name and create a mapped_gene_name field, create a new file, sort on this key and write out: 

tracking_id, mapped_gene_name,gene_id, length, coverage, FPKM, FPKM_conf_lo, FPKM_conf_hi, FPKM_status

2.) For isoforms.fpkm_tracking. Key on tracking_id only (but make sure these are distinct!), summarize to the gene level by adding all transcripts together based on the 'gene_id' column, fix the gene name and create a mapped_gene_name field, create a new file, sort on this key and write out:

tracking_id, gene_id, length, coverage, FPKM, FPKM_conf_lo, FPKM_conf_hi, FPKM_status, number of transcripts combined

3.) Create an output file sorted on FPKM - otherwise identical to the main file

4.) Create a filtered file sorted on FPKM for the top N% of genes

5.) Create some general figures showing the distribution of FPKM values, etc.

6.) Import a list of genes of interest (cancer genes + kinases + genes with neoplastic agents, etc.) - GOI list

7.) For each gene in the GOI list, show the expression level of this gene relative to all genes in the distribution

8.) Create files like the complete sorted and top N% sorted files, but only for the GOI genes

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
               desc => 'outdir does not exist or is not a directory: ' . $self->outdir,
      );
  }
  return @errors;
}


sub execute {
  my $self = shift;
  my $rnaseq_build = $self->build;
  my $cancer_annotation_db = $self->cancer_annotation_db;
  my $working_dir = $self->outdir;
  $working_dir .= "/" unless ($working_dir =~ /\/$/);

  my $cufflinks_dir = $rnaseq_build->data_directory . "/expression/";

  #Check the cufflinks dir
  unless (-e $cufflinks_dir){
    $self->error_message("Could not find cufflinks expression dir: $cufflinks_dir");
    die $self->error_message;
  }

  #Paths to input genes and isoforms files
  my $genes_infile = "$cufflinks_dir"."genes.fpkm_tracking";
  my $isoforms_infile = "$cufflinks_dir"."isoforms.fpkm_tracking";

  #Make sure expected data files exist
  unless (-e $genes_infile && -e $isoforms_infile){
    $self->error_message("File not found: $genes_infile | $isoforms_infile");
    die $self->error_message;
  }

  #Set paths to output files
  my $genes_file_sorted = "$working_dir"."genes.fpkm.namesort.tsv";
  my $isoforms_file_sorted = "$working_dir"."isoforms.fpkm.namesort.tsv";
  my $isoforms_merge_file_sorted = "$working_dir"."isoforms.merged.fpkm.namesort.tsv";
  my $goi_file = "$working_dir"."genes_of_interest.txt";

  #Get Entrez and Ensembl data for gene name mappings
  my $entrez_ensembl_data = $self->loadEntrezEnsemblData('-cancer_db' => $cancer_annotation_db);

  #Build a map of ensembl transcript ids to gene ids and gene names from the gene annotation object associated with the rna-seq builds
  my $reference_build = $rnaseq_build->reference_sequence_build;
  my $reference_build_id = $reference_build->id;
  my $reference_build_name = $reference_build->name;
  $self->status_message("Processing RNA-seq data that was aligned to: $reference_build_name");

  my $annotation_build = $rnaseq_build->model->annotation_build;
  my $annotation_build_name = $annotation_build->name;
  my $annotation_data_dir = $annotation_build->data_directory;
  my $transcript_info_path = $annotation_data_dir . "/annotation_data/rna_annotation/$reference_build_id-transcript_info.tsv";
  my $gtf_path = $annotation_build->annotation_file('gtf',$reference_build_id);
  $self->status_message("Getting transcript to gene and gene name mappings from annotation build: $annotation_build_name");
  unless (defined($gtf_path)) {
    $self->error_message("'There is no annotation GTF file defined for annotation_reference_transcripts build: ". $annotation_build->__display_name__);
    die $self->error_message;
  }
  $self->status_message("\t$gtf_path");
  unless (-e $transcript_info_path) {
    $self->error_message("'There is no transcript info file for annotation_reference_transcripts build: ". $annotation_build->__display_name__);
    die $self->error_message;
  }
  $self->status_message("\t$transcript_info_path");
  my $ensembl_map = $self->loadEnsemblMap('-gtf_path'=>$gtf_path, '-transcript_info_path'=>$transcript_info_path);

  #Import a set of gene symbol lists 
  #- these files must be gene symbols in the first column, .txt extension, tab-delimited if multiple columns, one symbol per field, no header
  #- fix gene names as they are being imported
  my $gene_symbol_lists_dir = $cancer_annotation_db->data_directory . "/GeneSymbolLists/";
  $gene_symbol_lists_dir = $self->checkDir('-dir'=>$gene_symbol_lists_dir, '-clear'=>"no");
  my @symbol_list_names = qw ( GenesOfInterest_MG );
  my $gene_symbol_lists = $self->importGeneSymbolLists('-gene_symbol_lists_dir'=>$gene_symbol_lists_dir, '-symbol_list_names'=>\@symbol_list_names, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
  my $goi_ref = $gene_symbol_lists->{'GenesOfInterest_MG'}->{'symbols'};
  open (GOI, ">$goi_file") || die "\n\nCould not open genes of interest file\n\n";
  foreach my $goi (sort keys %{$goi_ref}){
    print GOI "$goi\n";
  }
  close(GOI);

  #Parse the genes.fpkm_tracking and isoforms.fpkm_tracking files
  my $fpkm = $self->parseFpkmFile('-infile'=>$genes_infile, '-outfile'=>$genes_file_sorted, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-ensembl_map'=>$ensembl_map, '-verbose'=>$self->verbose);
  $fpkm = $self->parseFpkmFile('-infile'=>$isoforms_infile, '-outfile'=>$isoforms_file_sorted, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-ensembl_map'=>$ensembl_map, '-verbose'=>$self->verbose);

  #Merge the isoforms.fpkm_tracking file to the gene level
  my $merged_fpkm = $self->mergeIsoformsFile('-infile'=>$isoforms_infile, '-status_file'=>$genes_infile, '-outfile'=>$isoforms_merge_file_sorted, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-ensembl_map'=>$ensembl_map, '-verbose'=>$self->verbose);

  #Make the genes subdir
  my $images_sub_dir = $working_dir . "images/";
  mkdir($images_sub_dir);

  #Run the R code that produces filtered and sorted text files as well as plots for the whole dataset as well as individual genes of interest
  my $outlier_genes_absolute_r_script = __FILE__ . '.R';
  my $r_cmd_stdout = "$working_dir"."outlierGenesAbsolute.R.stdout";
  my $r_cmd_stderr = "$working_dir"."outlierGenesAbsolute.R.stderr";

  my $r_cmd = "$outlier_genes_absolute_r_script $working_dir $images_sub_dir " . $self->percent_cutoff . " 1>$r_cmd_stdout 2>$r_cmd_stderr";
  $self->status_message($r_cmd);
  Genome::Sys->shellcmd(cmd => $r_cmd);

  #Reorganize the output files into sub-directories
  my $genes_sub_dir = $self->createNewDir('-path'=>$working_dir, '-new_dir_name'=>"genes", '-force'=>"yes");
  my $genes_stats_dir = $self->createNewDir('-path'=>$genes_sub_dir, '-new_dir_name'=>"summary", '-force'=>"yes");
  my $isoforms_sub_dir = $self->createNewDir('-path'=>$working_dir, '-new_dir_name'=>"isoforms", '-force'=>"yes");
  my $isoforms_stats_dir = $self->createNewDir('-path'=>$isoforms_sub_dir, '-new_dir_name'=>"summary", '-force'=>"yes");
  my $isoforms_merged_sub_dir = $self->createNewDir('-path'=>$working_dir, '-new_dir_name'=>"isoforms_merged", '-force'=>"yes");
  my $isoforms_merged_stats_dir = $self->createNewDir('-path'=>$isoforms_merged_sub_dir, '-new_dir_name'=>"summary", '-force'=>"yes");

  my $mv_cmd1 = "mv $working_dir/genes.fpkm*.tsv $genes_sub_dir";
  Genome::Sys->shellcmd(cmd => $mv_cmd1);
  my $mv_cmd2 = "mv $working_dir/isoforms.fpkm*.tsv $isoforms_sub_dir";
  Genome::Sys->shellcmd(cmd => $mv_cmd2);
  my $mv_cmd3 = "mv $working_dir/isoforms.merged.fpkm*.tsv $isoforms_merged_sub_dir";
  Genome::Sys->shellcmd(cmd => $mv_cmd3);

  #Create some basic summary statistics for each expression file generated
  my $genes_data_file = "$genes_sub_dir"."genes.fpkm.namesort.tsv";
  my $genes_stat_file = "$genes_stats_dir"."Stats.tsv";
  my $isoforms_data_file = "$isoforms_sub_dir"."isoforms.fpkm.namesort.tsv";
  my $isoforms_stat_file = "$isoforms_stats_dir"."Stats.tsv";
  my $isoforms_merge_data_file = "$isoforms_merged_sub_dir"."isoforms.merged.fpkm.namesort.tsv";
  my $isoforms_merge_stat_file = "$isoforms_merged_stats_dir"."Stats.tsv";

  #Each of these files above has an 'FPKM' column.  All basic Cufflinks stats will be based on that
  &calculateCufflinksStats('-infile'=>$genes_data_file, '-outfile'=>$genes_stat_file);
  &calculateCufflinksStats('-infile'=>$isoforms_data_file, '-outfile'=>$isoforms_stat_file);
  &calculateCufflinksStats('-infile'=>$isoforms_merge_data_file, '-outfile'=>$isoforms_merge_stat_file);

  #Set outputs
  my $tumor_fpkm_file = "$isoforms_merged_sub_dir"."isoforms.merged.fpkm.expsort.tsv";
  unless (-e $tumor_fpkm_file){
    die $self->error_message("Trying to set a file as output but the file does not exist: $tumor_fpkm_file");
  }
  $self->tumor_fpkm_file($tumor_fpkm_file);

  my $tumor_fpkm_topnpercent_file = "$isoforms_merged_sub_dir"."isoforms.merged.fpkm.expsort.top" . $self->percent_cutoff ."percent.tsv";
  unless (-e $tumor_fpkm_topnpercent_file){
    die $self->error_message("Trying to set a file as output but the file does not exist: $tumor_fpkm_topnpercent_file");
  }
  $self->tumor_fpkm_topnpercent_file($tumor_fpkm_topnpercent_file);

  return 1;
}

1;
