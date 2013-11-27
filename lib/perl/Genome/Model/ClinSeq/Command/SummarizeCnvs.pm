package Genome::Model::ClinSeq::Command::SummarizeCnvs;

#Written by Malachi Griffith

use strict;
use warnings;
use Genome;
use Data::Dumper;
use Genome::Model::ClinSeq::Util qw(:all);

class Genome::Model::ClinSeq::Command::SummarizeCnvs {
    is => 'Command::V2',
    has_input => [
        build => { 
              is => 'Genome::Model::Build::SomaticVariation',
              is_many => 0,
              shell_args_position => 1,
              require_user_verify => 0,
              doc => 'somatic variation build to summarize CNVs from',
        },
        cnv_hmm_file => {
            is => 'FilesystemPath',
            doc => 'cnv hmm file from clin-seq generate-clonality-plots',
        },
        gene_amp_file => {
            is => 'FilesystemPath',
            doc => 'gene amplification file from clin-seq run-cn-view',
        },
        gene_del_file => {
            is => 'FilesystemPath',
            doc => 'gene deletion file from clin-seq run-cn-view',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
    ],
    doc => 'summarize the CNVs of clinseq build',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq summarize-cnvs --outdir=/tmp/  --cnv-hmm-file=? --gene-amp-file=? --gene-del-file=? 119390903 

genome model clin-seq summarize-cnvs --outdir=/tmp/  --cnv-hmm-file=? --gene-amp-file=? --gene-del-file=? id=119390903

genome model clin-seq summarize-cnvs --outdir=/tmp/  --cnv-hmm-file=? --gene-amp-file=? --gene-del-file=? model.id=2882504846

genome model clin-seq summarize-cnvs --outdir=/tmp/  --cnv-hmm-file=? --gene-amp-file=? --gene-del-file=? "model.name='My Somatic model name'"

genome model clin-seq summarize-cnvs --outdir=/tmp/  --cnv-hmm-file=/gscmnt/gc7001/info/model_data/2887519760/build126680687/AML103/clonality/cnaseq.cnvhmm  --gene-amp-file=/gscmnt/gc7001/info/model_data/2887519760/build126680687/AML103/cnv/cnv.AllGenes_Ensembl58.amp.tsv  --gene-del-file=/gscmnt/gc7001/info/model_data/2887519760/build126680687/AML103/cnv/cnv.AllGenes_Ensembl58.del.tsv  119390903

EOS
}

sub help_detail {
    return <<EOS
Summarize copy number variants using a somatic-variation build and files from a clinseq build 

(put more content here)
EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless (-e $self->outdir && -d $self->outdir) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['outdir'],
	                                          desc => "Outdir: " . $self->outdir . " not found or not a directory",
                                          );
  }
  unless (-e $self->cnv_hmm_file) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['cnv_hmm_file'],
	                                          desc => "cnv hmm file: " . $self->cnv_hmm_file . " not found",
                                          );
  }
  unless (-e $self->gene_amp_file) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['gene_amp_file'],
	                                          desc => "gene amplification file: " . $self->gene_amp_file . " not found",
                                          );
  }
  unless (-e $self->gene_del_file) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['gene_del_file'],
	                                          desc => "gene deletion file: " . $self->gene_del_file . " not found",
                                          );
  }

  return @errors;
}

sub execute {
  my $self = shift;
  my $build = $self->build;
  my $outdir = $self->outdir;

  unless ($outdir =~ /\/$/){
    $outdir .= "/";
  }

  #Determine the reference alignment build version
  my $reference_build_ncbi_n = "";
  my $reference_build = $build->reference_sequence_build;
  my $reference_build_name = $reference_build->name;
  if ($reference_build_name =~ /GRCh37\-lite\-build37/){
    $reference_build_ncbi_n = "37";
  }else{
    die $self->error_message("Reference sequence build name not recognized by clin-seq summarize-cnvs");
  }

  #Create a Stats.tsv and Summarize the number of CNV amp and del windows
  #Question Answer  Data_Type Analysis_Type Statistic_Type  Extra_Description
  my $stats_file = $outdir . "Stats.tsv";
  open (STATS, ">$stats_file") || die "\n\nCould not open stats file: $stats_file\n\n";
  print STATS "Question\tAnswer\tData_Type\tAnalysis_Type\tStatistic_Type\tExtra_Description\n";
  my $wgs_som_build_dir = $build->data_directory;

  #Overall strategy
  #Get a copy of the cnvs.hq (10-kb copy number window values) and cna-seq (cnv-hmm files) and summarize
  my $cnv_hq = $wgs_som_build_dir . "/variants/cnvs.hq";
  my $cnv_hq_new = $outdir . "cnvs.hq";
  if (-e $cnv_hq and not -e $cnv_hq_new){
    Genome::Sys->copy_file($cnv_hq, $cnv_hq_new);
  }

  #Similarly, copy the .png file generated by the pipeline if it is available
  #/gscmnt/gc8002/info/model_data/2882504846/build119390903/variants/cnv/bam-to-cna-v1-d41d8cd98f00b204e9800998ecf8427e/cnvs.hq.png
  my $cnv_png_search = $wgs_som_build_dir . "/variants/cnv/bam-to-cna*/cnvs.hq.png";
  my $cnv_png = `ls $cnv_png_search 2>/dev/null`;
  chomp($cnv_png);
  my $cnv_png_new = $outdir . "cnvs.hq.png";
  my $cp_cmd2 = "cp $cnv_png $outdir";
  if (-e $cnv_png and not -e $cnv_png_new){
    Genome::Sys->copy_file($cnv_png, $cnv_png_new);
  }

  #Gather CNV window stats:
  #Number of windows
  #Number windows with CNV diff >1, >2, >5, <-0.5, <-0.75, <-1.0
  #Number of windows with tumor/normal coverage >100 and >1000
  my ($window_count, $cum_cov, $amp_050, $amp_1, $amp_2, $amp_5, $amp_10, $del_025, $del_050, $del_075, $del_1, $del_15, $cov_100x, $cov_250x, $cov_1000x, $cov_2500x, $cov_5000x, $cov_10000x, $cov_25000x, $cov_50000x) = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
  if (-e $cnv_hq){
    open (CNV, "$cnv_hq") || die "\n\nCould not open cnv hq file: $cnv_hq\n\n";
    while(<CNV>){
      chomp($_);
      if ($_ =~ /^\#|^CHR/){
        next();
      }
      my @line = split("\t", $_);
      my $chr = $line[0];
      my $pos = $line[1];
      my $tumor = $line[2];
      my $normal = $line[3];
      my $cov = $tumor+$normal;
      my $diff = $line[4];
      $window_count++;
      $cum_cov+=$cov;
      if ($diff > 0.5){$amp_050++;}
      if ($diff > 1){$amp_1++;}
      if ($diff > 2){$amp_2++;}
      if ($diff > 5){$amp_5++;}
      if ($diff > 10){$amp_10++;}
      if ($diff < -0.25){$del_025++;}
      if ($diff < -0.5){$del_050++;}
      if ($diff < -0.75){$del_075++;}
      if ($diff < -1){$del_1++;}
      if ($diff < -1.5){$del_15++;}
      if ($cov > 100){$cov_100x++;}
      if ($cov > 250){$cov_250x++;}
      if ($cov > 1000){$cov_1000x++;}
      if ($cov > 2500){$cov_2500x++;}
      if ($cov > 5000){$cov_5000x++;}
      if ($cov > 10000){$cov_10000x++;}
      if ($cov > 25000){$cov_25000x++;}
      if ($cov > 50000){$cov_50000x++;}
     }
    close(CNV);
  }
  my $avg_cov = sprintf("%.2f", $cum_cov/$window_count);
  my $cov_100x_p = sprintf("%.2f", ($cov_100x/$window_count)*100);
  my $cov_250x_p = sprintf("%.2f", ($cov_250x/$window_count)*100);
  my $cov_1000x_p = sprintf("%.2f", ($cov_1000x/$window_count)*100);
  my $cov_2500x_p = sprintf("%.2f", ($cov_2500x/$window_count)*100);
  my $cov_5000x_p = sprintf("%.2f", ($cov_5000x/$window_count)*100);
  my $cov_10000x_p = sprintf("%.2f", ($cov_10000x/$window_count)*100);
  my $cov_25000x_p = sprintf("%.2f", ($cov_25000x/$window_count)*100);
  my $cov_50000x_p = sprintf("%.2f", ($cov_50000x/$window_count)*100);

  print STATS "Total CNV window count\t$window_count\twgs\tcnv_bamtocna\tcount\tNumber of windows from CNV analysis\n";
  print STATS "CNV amplified windows > 0.5\t$amp_050\twgs\tcnv_bamtocna\tcount\tTotal CNV window counts with tumor-normal diff > 0.5\n";
  print STATS "CNV amplified windows > 1\t$amp_1\twgs\tcnv_bamtocna\tcount\tTotal CNV window counts with tumor-normal diff > 1\n";
  print STATS "CNV amplified windows > 2\t$amp_2\twgs\tcnv_bamtocna\tcount\tTotal CNV window counts with tumor-normal diff > 2\n";
  print STATS "CNV amplified windows > 5\t$amp_5\twgs\tcnv_bamtocna\tcount\tTotal CNV window counts with tumor-normal diff > 5\n";
  print STATS "CNV amplified windows > 10\t$amp_10\twgs\tcnv_bamtocna\tcount\tTotal CNV window counts with tumor-normal diff > 10\n";
  print STATS "CNV deleted windows < -0.25\t$del_025\twgs\tcnv_bamtocna\tcount\tTotal CNV window counts with tumor-normal diff < -0.25\n";
  print STATS "CNV deleted windows < -0.50\t$del_050\twgs\tcnv_bamtocna\tcount\tTotal CNV window counts with tumor-normal diff > -0.50\n";
  print STATS "CNV deleted windows < -0.75\t$del_075\twgs\tcnv_bamtocna\tcount\tTotal CNV window counts with tumor-normal diff > -0.75\n";
  print STATS "CNV deleted windows < -1.0\t$del_1\twgs\tcnv_bamtocna\tcount\tTotal CNV window counts with tumor-normal diff > -1.0\n";
  print STATS "CNV deleted windows < -1.5\t$del_15\twgs\tcnv_bamtocna\tcount\tTotal CNV window counts with tumor-normal diff > -1.5\n";
  print STATS "Average coverage of CNV windows\t$avg_cov\twgs\tcnv_bamtocna\tmean\tCumulative coverage of tumor+normal divided by window count\n";
  print STATS "CNV windows with coverage > 100x\t$cov_100x_p\twgs\tcnv_bamtocna\tpercent\tCNV windows with tumor+normal coverage > 100 reads\n";
  print STATS "CNV windows with coverage > 250x\t$cov_250x_p\twgs\tcnv_bamtocna\tpercent\tCNV windows with tumor+normal coverage > 250 reads\n";
  print STATS "CNV windows with coverage > 1000x\t$cov_1000x_p\twgs\tcnv_bamtocna\tpercent\tCNV windows with tumor+normal coverage > 1000 reads\n";
  print STATS "CNV windows with coverage > 2500x\t$cov_2500x_p\twgs\tcnv_bamtocna\tpercent\tCNV windows with tumor+normal coverage > 2500 reads\n";
  print STATS "CNV windows with coverage > 5000x\t$cov_5000x_p\twgs\tcnv_bamtocna\tpercent\tCNV windows with tumor+normal coverage > 5000 reads\n";
  print STATS "CNV windows with coverage > 10000x\t$cov_10000x_p\twgs\tcnv_bamtocna\tpercent\tCNV windows with tumor+normal coverage > 10000 reads\n";
  print STATS "CNV windows with coverage > 25000x\t$cov_25000x_p\twgs\tcnv_bamtocna\tpercent\tCNV windows with tumor+normal coverage > 25000 reads\n";
  print STATS "CNV windows with coverage > 50000x\t$cov_50000x_p\twgs\tcnv_bamtocna\tpercent\tCNV windows with tumor+normal coverage > 50000 reads\n";

  #Summarise the CNV AMP and DEL genes from the CNView analyses stored in the ClinSeq results
  my $cnv_amp = $self->gene_amp_file;
  my $cnv_del = $self->gene_del_file;
  my $amp_count = -1;
  my $del_count = -1;
  if (-e $cnv_amp && -e $cnv_del){
    open (AMP, "$cnv_amp") || die "\n\nCould not open amp file: $cnv_amp\n\n";
    while(<AMP>){
      $amp_count++;
    }
    close(DEL);
    open (DEL, "$cnv_del") || die "\n\nCould not open del file: $cnv_del\n\n";
    while(<DEL>){
      $del_count++;
    }
    close(DEL);
    print STATS "CNV amplified genes\t$amp_count\twgs\tcnv_cnview\tcount\tNumber of CNV tumor vs. normal amplified genes according to CNView analysis\n";
    print STATS "CNV deleted genes\t$del_count\twgs\tcnv_cnview\tcount\tNumber of CNV tumor vs. normal deleted genes according to CNView analysis\n";
  }
  #Summarize the number of CNV amp and del segments from the hmm-segs file
  #Unfortunately, for now this is only available from the ClinSeq build itself... it would be better to get all this from somatic variation probably...
  my $cnv_hmm = $self->cnv_hmm_file;
  if (-e $cnv_hmm){
    #Gather some basic stats from the cna-seg analysis
    my $amp_seg_count = 0;
    my $del_seg_count = 0;

    open (CNV_HMM, "$cnv_hmm") || die "\n\nCould not open CNV HMM file: $cnv_hmm\n\n";
    #Kind of a nasty format to these files.  Search for data entries like the following:
    #CN1 = Tumor?  and  CN2 = Normal? 
    #CHR	START	END	SIZE	nMarkers	CN1	Adjusted_CN1	CN2	Adjusted_CN2	LLR_Somatic	Status
    #4	7190000	7540000	350000	36	2	1.57	1	1.46	10.25	Gain
    while(<CNV_HMM>){
      chomp($_);
      next if ($_ =~ /^\#/);
      my @line = split("\t", $_);
      next unless (scalar @line == 11);

      #print "@line\n";
      my $chr = $line[0];
      my $start = $line[1];
      my $end = $line[2];
      my $size = $line[3];
      my $nmarkers = $line[4];
      my $cn1 = $line[5];
      my $cn1_adjusted = $line[6];
      my $cn2 = $line[7];
      my $cn2_adjusted = $line[8];
      my $llr_somatic = $line[9];
      my $status = $line[10];
      if ($status eq 'Gain'){
        $amp_seg_count++;
      }
      if ($status eq 'Loss'){
        $del_seg_count++;
      }
    }
    close (CNV_HMM);
    print STATS "CNV amplified segments\t$amp_seg_count\twgs\tcnv_cnaseq\tcount\tNumber of CNV tumor vs. normal amplified segments according to CNV hmm cna-seq analysis\n";
    print STATS "CNV deleted segments\t$del_seg_count\twgs\tcnv_cnaseq\tcount\tNumber of CNV tumor vs. normal deleted segments according to CNV hmm cna-seg analysis\n";
  }
  close (STATS);

  #Perform single BAM CNV analysis if it has not already been run.  Otherwise just copy over the resulting PDF
  my $single_bam_cnv_dir = "$outdir"."single_bam_cnv/";
  mkdir ($single_bam_cnv_dir);
  my $output_pdf_name = "CNV_SingleBAMs_TumorAndNormal.pdf";
  my $output_pdf_path = $single_bam_cnv_dir . $output_pdf_name;

  #First test to see if a single BAM CNV plot.pdf is already created (done by more recent versions of the somatic variation pipeline)
  my $test_path = $wgs_som_build_dir . "/variants/cnv/plot-cnv*/cnv_graph.pdf";
  my $pdf_path = `ls $test_path  2>/dev/null`;
  chomp($pdf_path);
  if (-e $pdf_path){
    unless (-e $output_pdf_path){
      Genome::Sys->copy_file($pdf_path, $output_pdf_path);
    }
  }else{
    #The .pdf file was not found.  Presumably this is an older somatic variation build that did not include this step. Generate it now
    my $cn_stdout = "$single_bam_cnv_dir"."CNV_SingleBAMs_TumorAndNormal.stdout";
    my $cn_stderr = "$single_bam_cnv_dir"."CNV_SingleBAMs_TumorAndNormal.stderr";
    my $normal_bam = $build->normal_bam;
    my $tumor_bam = $build->tumor_bam;
    my $single_bam_cnv_plot_cmd = Genome::Model::Tools::CopyNumber::PlotSegmentsFromBamsWorkflow->create(normal_bam=>$normal_bam, tumor_bam=>$tumor_bam, output_directory=>$single_bam_cnv_dir, genome_build=>$reference_build_ncbi_n, output_pdf=>$output_pdf_name);
    $single_bam_cnv_plot_cmd->execute();

    unless (-e $output_pdf_path){
      $single_bam_cnv_plot_cmd .= " 1>$cn_stdout 2>$cn_stderr";
      Genome::Sys->shellcmd(cmd => $single_bam_cnv_plot_cmd);
    }
  }

  $self->status_message("\n\n");

  return 1;
}

1;


