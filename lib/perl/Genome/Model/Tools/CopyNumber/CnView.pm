#Written by Malachi Griffith, Scott Smith, Obi Griffith

package Genome::Model::Tools::CopyNumber::CnView;
use strict;
use warnings;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use File::Basename;
use Genome::Model::ClinSeq::Util qw(:all);

class Genome::Model::Tools::CopyNumber::CnView{
  is=>'Command::V2',
  has=>[
    annotation_build  => { 
      is => 'Genome::Model::Build::ImportedAnnotation', 
      doc => 'Supply an annotation build id (e.g., 124434505 for NCBI-human.ensembl/67_37l_v2)' },
    cancer_annotation_db => {
      is => 'Genome::Db',
      example_values => ['tgi/cancer-annotation/human/build37-20130401.1'],
      doc => 'cancer-specific annotation extenions',
    },
    output_dir        => { 
      is => 'Text', 
      doc => 'The output directory' },
    segments_file     => { 
      is => 'Text', 
      doc => 'Path to cnaseq.cnvhmm file' },
  ],
  has_optional_input  => [
    cnv_file          => { 
      is => 'Text', 
      doc => 'Path to cnvs.hq file' },
    somatic_build       => { 
      is => 'Genome::Model::Build::SomaticVariation', 
      doc => 'Instead of supplying a "cnvs.hq" file, you can specify a somatic variation build ID and the file will be found automatically' },
    gene_targets_file => { 
      is => 'Text', 
      doc => 'List of gene symbols. These will be highlighted if they meet the minimum CNV cutoffs.  If not specified all genes will be used by default.' },
    name              => { 
      is => 'Text', 
      doc => 'Human readable name describing the gene symbols to be highlighted. Use "All" if no gene_targets_file is specified' },
    cnv_results_file  => { 
      is => 'Text', 
      doc => 'If you have already generated a CNV results file for a list of genes, you can supply it and skip to the graph generation step' },
    ideogram_file     => { 
      is => 'Text', 
      doc => 'Chromosome ideogram coordinates for the reference build.  Will be automatically selected based on the reference build specified' },
    chr               => { 
      is => 'Text', 
      doc => 'If you wish to limit analysis to a single chromosome, use --chr (e.g. --chr=17). Otherwise all will be processed' },
    window_size         => { 
      is => 'Number', 
      doc => 'For determining overlap, you can define a window size' },
    chr_start         => { 
      is => 'Number', 
      doc => 'If you specify a single chromosome, you may also specify start position using --chr-start' },
    chr_end           => { 
      is => 'Number', 
      doc => 'If you specify a single chromosome, you may also specify end position using --chr-end' },
  ],
  has_optional_param => [
    verbose => { 
      is => 'Boolean', 
      doc => 'More verbose output' },
    image_type => { 
      is => 'Text', 
      valid_values => ['jpeg', 'png', 'bmp', 'tiff', 'pdf', 'none'], default_value => 'jpeg', 
      doc => 'Type of image files created.' },
  ],
  doc=>'Generate summaries, visualizations, and gene-level results for CNV output'
};

sub help_synopsis {
  return <<EOS
gmt copy-number cn-view --annotation-build=124434505 --cancer_annotation_db='tgi/cancer-annotation/human/build37-20130401.1' --cnv-file=/gscmnt/gc13001/info/model_data/2888915570/build129973671/variants/cnvs.hq --segments-file=/gscmnt/gc2013/info/model_data/2889110844/build130030495/PNC6/clonality/cnaseq.cnvhmm --output-dir=/tmp/ --gene-targets-file=/gscmnt/sata132/techd/mgriffit/reference_annotations/GeneSymbolLists/CancerGeneCensusPlus_Sanger.txt --name='Cancer Genes'  --chr=1  --verbose
EOS
}

sub help_detail {
  return <<EOS
Takes standard CNV analysis output, together with CNAseg output to produce chromosome-by-chromosome and gene level colored copy number plots and summary files.

EOS
}

sub execute{
  my $self = shift;

  #Required
  my $annotation_build = $self->annotation_build;
  my $cancer_annotation_db = $self->cancer_annotation_db;

  my $cnv_file = $self->cnv_file;
  my $segments_file = $self->segments_file;
  my $somatic_build = $self->somatic_build;
  my $output_dir = $self->output_dir;
  my $name = $self->name;

  #Optional
  my $gene_targets_file = $self->gene_targets_file;
  my $cnv_results_file = $self->cnv_results_file;
  my $ideogram_file = $self->ideogram_file;
  my $chr = $self->chr;
  my $chr_start = $self->chr_start;
  my $chr_end = $self->chr_end;
  my $image_type = $self->image_type;
  my $window_size = $self->window_size;

  #Set amplification / deletion cutoffs that will be used to created filtered files of results for convenience in downstream analyses
  #Try two fold gain and 1/2 loss to start
  my $amp_cutoff = 2;
  my $del_cutoff = -0.5;

  #Check/format all input parameters and options
  #Based on the reference build specified, define paths to transcript/gene id mapping files
  my ($gtf_file, $subdir, $outfile_prefix, $outfile_cnvhmm, $name_f);
  my $extra_flank = 100000;

  ####################################################################################################################################
  #Check inputs
  my $inputs_are_good = eval { 
    unless ($cnv_file || $somatic_build){
      $self->usage_message($self->help_usage_complete_text);
      $self->error_message("Must supply either cnv_file or somatic_build (Somatic Variation)");
      return;
    }

    #Check annotation build and clinseq annotations dir
    my $clinseq_annotations_dir = $cancer_annotation_db->data_directory;
    my $annotation_data_dir=$annotation_build->data_directory;
    my $reference_sequence_build=$annotation_build->reference_sequence;
    my $default_build37 = Genome::Model::Build->get(106942997);
    my $default_build36 = Genome::Model::Build->get(101947881);
    if ($reference_sequence_build->is_compatible_with($default_build37)){
      $ideogram_file ||= $clinseq_annotations_dir . "/hg19/ideogram/ChrBandIdeogram.tsv";
    }elsif ($reference_sequence_build->is_compatible_with($default_build36)){
      $ideogram_file ||= $clinseq_annotations_dir . "/hg18/ideogram/ChrBandIdeogram.tsv";
    }else {
      $self->error_message("Specified reference build resolved from annotation build is not compatible with default build36 or build37");
      return;
    }

    #Check GTF file
    my $reference_id = $reference_sequence_build->id;
    $gtf_file = "$annotation_data_dir/annotation_data/rna_annotation/$reference_id-all_sequences.gtf";
    unless (-e $gtf_file && -e $ideogram_file){
      $self->error_message("One or more of the following annotation files is missing:\ngtf_file = $gtf_file\nideogram_file = $ideogram_file");
      return;
    }

    #Check gene targets file if defined
    if ($self->gene_targets_file){    
      unless (-e $gene_targets_file){
        $self->error_message("Gene targets file not found: $gene_targets_file");
        return;
      }
    }

    #Check for segments file
    unless (-e $segments_file){  
      $self->error_message("Segments file missing: $segments_file");
      return;
    }

    #Check input CNV data
    if ($cnv_file){
      unless (-e $cnv_file){
        $self->error_message("CNV file missing: $cnv_file");
        return;
      }
    }elsif($somatic_build){
      my $somatic_directory = $somatic_build->data_directory;
      my $path = "$somatic_directory/variants/"."cnvs.hq";
      unless(-e $path){
        $self->error_message("Could not find a copy number results file in the directory: $somatic_directory");
        return;
      }
      $cnv_file = $path;
    }

    #Check output dir
    unless (-e $output_dir && -d $output_dir){
      $self->status_message("Creating dir: $output_dir") if $self->verbose;
      mkdir($output_dir);
    }
    unless ($output_dir =~ /\/$/){
      $output_dir .= "/";
    }
    
    #Get/Set the name of the comparison
    if ($name){
      $name_f = $name;
      $name_f =~ s/ //g;
    }else{
      $self->error_message("Must specify a valid --name parameter");
      return;
    }
    
    #Create a subdirectory within the working directory
    $subdir = "$output_dir"."CNView_"."$name_f/";
    unless (-e $subdir){
      mkdir($subdir);
    }
   
    #Check the chromosome specified by the user for format (should be N).   If not specified, set $chr to 'ALL' and set $chr_start and $chr_end to 0
    chomp($chr) if $chr;
    chomp($chr_start) if $chr_start;
    chomp($chr_end) if $chr_end;
    if ($chr_start || $chr_end){
      if (($chr_start =~ /^\d+$/) && ($chr_end =~ /^\d+$/)){
        unless ($chr_start <= $chr_end){
          $self->error_message("--chr_start ($chr_start) must be less than --chr_end ($chr_end)");
          return;
        }
        $chr_start -= $extra_flank;
        $chr_end += $extra_flank;
        if ($chr_start < 1){
          $chr_start = 1;
        }
      }else{
        $self->error_message("--chr_start ($chr_start) and --chr_end ($chr_end) must both be integers");
        return;
      }
      unless($chr){
        $self->error_message("If you are going to define --chr_start and --chr_end, you must also define a chromosome with --chr\n\n");
        return;
      }
    }else{
      $chr_start = 0;
      $chr_end = 0;
    }

    if ($chr){
      if ($chr =~ /all/i){
        $chr = "ALL";
      }
    }else{
      $chr = "ALL";
      $chr_start = 0;
      $chr_end = 0;
    }

    #Name the output files
    $outfile_cnvhmm = "$subdir"."cnaseq.cnvhmm.tsv";
    if ($chr eq "ALL"){
      $outfile_prefix = "$subdir"."CNView_"."$name_f";
    }elsif ($chr_start && $chr_end){
      $outfile_prefix = "$subdir"."CNView_"."$name_f"."_chr"."$chr"."_"."$chr_start-$chr_end";
    }else{
      $outfile_prefix = "$subdir"."CNView_"."$name_f"."_chr"."$chr";
    }

    #If specified check the cnv results file
    if ($cnv_results_file){
      unless (-e $cnv_results_file){
        $self->error_message("CNV results file could not be found: $cnv_results_file");
        return;
      }
    }
    #Check for valid image type.  If not specified, default to jpeg
    chomp($image_type);
    if ($image_type){
      unless ($image_type =~ /^jpeg$|^png$|^bmp$|^tiff$|^pdf$|^none$/){
        $self->error_message("image_type: ($image_type) not recognized.  Must be one of: 'jpeg', 'png', 'bmp', 'tiff', 'pdf', 'none'");
        return;
      }
    }else{
      $image_type = "jpeg";
    }
    return 1;
  };

  #die if the inputs are invalid
  unless ($inputs_are_good) {
    die $self->error_message("error processing inputs!");
  }

  #Done checking inputs
  ####################################################################################################################################

  #Load ensembl/entrez data for fixing gene names
  my $entrez_ensembl_data = $self->loadEntrezEnsemblData(-cancer_db => $cancer_annotation_db);

  #Import ideogram data
  my $ideo_data = $self->importIdeogramData('-ideogram_file'=>$ideogram_file);

  #If the user is supplying a pre-computed CNV file, the following steps will be skipped
  my ($g_map, $t_map, $cnvs, $segments, $targets);
  my @genes_results_files;
  my @trans_results_files;

  if ($cnv_results_file){
    $self->status_message("Using user supplied CNV results file: $cnv_results_file") if $self->verbose;
  }else{
    $self->status_message("Generating CNV results file") if $self->verbose;

    #Load the gene targets of interest
    if ($self->gene_targets_file){
      $targets = $self->loadTargetGenes('-gene_targets_file'=>$gene_targets_file, '-entrez_ensembl_data'=>$entrez_ensembl_data);
    }

    #Load the gene to transcript name mappings.
    #Key on unique gene ID - reference to a hash containing all transcript IDs associated with that gene
    my $r = $self->loadGeneTranscriptMap('-gtf_file'=>$gtf_file, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-targets'=>$targets, '-ideo_data'=>$ideo_data, '-target_chr'=>$chr, '-target_chr_start'=>$chr_start, '-target_chr_end'=>$chr_end);
    $g_map = $r->{genes};
    $t_map = $r->{trans};

    #Load the CNV window coordinates and copy number estimates
    $cnvs = $self->loadCnvData('-cnv_file'=>$cnv_file, '-target_chr'=>$chr, '-target_chr_start'=>$chr_start, '-target_chr_end'=>$chr_end, '-window_size'=>$window_size);
    
    #Load the CNV segments and copy number estimates
    $segments = $self->loadSegmentData('-segments_file'=>$segments_file, '-outfile_cnvhmm'=>$outfile_cnvhmm);

    #Get CNV status and CNV value for each gene
    $self->status_message("Determining CNV status for genes") if $self->verbose;
    $self->calculateMeanCnvDiff('-f_map'=>$g_map, '-f_name'=>'genes', '-cnvs'=>$cnvs, '-segments'=>$segments, '-output_dir'=>$subdir);

    #Get CNV status and CNV value for each gene
    $self->status_message("Determining CNV status for transcripts") if $self->verbose;
    $self->calculateMeanCnvDiff('-f_map'=>$t_map, '-f_name'=>'transcripts', '-cnvs'=>$cnvs, '-segments'=>$segments, '-output_dir'=>$subdir);
    
    #Print gene results to a series of output files
    @genes_results_files = @{$self->printCnvResultFile('-outfile_prefix'=>$outfile_prefix."_genes", '-f_map'=>$g_map, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-amp_cutoff'=>$amp_cutoff, '-del_cutoff'=>$del_cutoff)};

    #Print transcript results to a series of output files
    @trans_results_files = @{$self->printCnvResultFile('-outfile_prefix'=>$outfile_prefix."_transcripts", '-f_map'=>$t_map, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-amp_cutoff'=>$amp_cutoff, '-del_cutoff'=>$del_cutoff)};

  }

  #Only the genes CNV results file will be fed into the R script for generating plots
  $cnv_results_file = $genes_results_files[0];


  #Assume that the CNview.R script resides in the same dir as this script and obtain the path automatically

  #Execute the R script that generates the CNV plots
  my $chr_name=$chr; #R script still expects chr number to be prefixed with "chr"
  unless ($chr_name eq 'ALL'){$chr_name="chr".$chr;}
  my $rscript = __FILE__.".R";
  my $r_cmd = "$rscript '$name' $cnv_file $cnv_results_file $outfile_cnvhmm $ideogram_file $subdir $chr_name $chr_start $chr_end $image_type";
  my $r_cmd_stdout = "$subdir"."CNView.R.stdout";
  my $r_cmd_stderr = "$subdir"."CNView.R.stderr";
  if ($self->verbose){
    $self->status_message("Executing R code:\n$r_cmd") if $self->verbose;
  }else{
    $r_cmd .= " 1>$r_cmd_stdout 2>$r_cmd_stderr";
  }

  Genome::Sys->shellcmd(cmd => $r_cmd);

  $self->status_message("Results written to:\n$subdir") if $self->verbose;

  return(1);
}


#####################################################################################################################################################################
#Load the gene targets of interest
#####################################################################################################################################################################
sub loadTargetGenes{
  my $self = shift;
  my %args = @_;
  my $gene_targets_file = $args{'-gene_targets_file'}; 
  my $entrez_ensembl_data = $args{'-entrez_ensembl_data'};

  my %targets_gene_name;
  my %targets_mapped_gene_name;
  open (TARGET ,"$gene_targets_file") || die "\n\nCould not open gene targets file\n\n";
  while(<TARGET>){
    chomp($_);
    if ($_ =~ /(\S+)/){
      my $gene_name = $1;
      #Skip "n/a"
      if ($gene_name eq "n/a"){
        next();
      }else{
        my $mapped_gene_name = $self->fixGeneName('-gene'=>$gene_name, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);      
        $targets_gene_name{$gene_name} = 1;
        $targets_mapped_gene_name{$mapped_gene_name} = 1;
      }
    }
  }
  close(TARGET);
  my %targets;
  $targets{gene_name} = \%targets_gene_name;
  $targets{mapped_gene_name} = \%targets_mapped_gene_name;
  return(\%targets);
}


#####################################################################################################################################################################
#Load a genes object including gene to transcript mappings, also load a transcript object, key on unique gid or tid                                                 #
#####################################################################################################################################################################
sub loadGeneTranscriptMap{
  my $self = shift;
  my %args = @_;
  my $gtf_file = $args{'-gtf_file'};
  my $entrez_ensembl_data = $args{'-entrez_ensembl_data'};
  my $targets = $args{'-targets'};
  my $ideo_data = $args{'-ideo_data'};
  my $target_chr = $args{'-target_chr'};
  my $target_chr_start = $args{'-target_chr_start'};
  my $target_chr_end = $args{'-target_chr_end'};
  
  my $targets_gene_name = $targets->{gene_name};
  my $targets_mapped_gene_name = $targets->{mapped_gene_name};
 
  my %genes;
  my %trans;

  open (GTF, "$gtf_file") || die "\n\nCould not open input file: $gtf_file\n\n";
  while(<GTF>){
    chomp($_);
    my @line = split("\t", $_);
    next unless ($line[2] eq "exon");
    my $chr = $line[0];
    my $chr_start = $line[3];
    my $chr_end = $line[4];
    my $id_string = $line[8];
    my @ids = split(";", $id_string);
    my $gene_name_string = $ids[0];
    my $gene_id_string = $ids[1];
    my $transcript_id_string = $ids[2];
    my $gene_name;
    my $gid;
    my $tid;
    if ($gene_name_string =~ /gene\_name\s+\"(.*)\"/){
      $gene_name = $1;
    }else{
      $self->error_message("Could not resolve gene name from GTF file and string:\n\t$gtf_file\n\t$id_string");
      exit 1;
    }
    if ($gene_id_string =~ /gene\_id\s+\"(.*)\"/){
      $gid = $1;
    }else{
      $self->error_message("Could not resolve gene id from GTF file and string:\n\t$gtf_file\n\t$id_string");
      exit 1;
    }
    if ($transcript_id_string =~ /transcript\_id\s+\"(.*)\"/){
      $tid = $1;
    }else{
      $self->error_message("Could not resolve transcript id from GTF file and string:\n\t$gtf_file\n\t$id_string");
      exit 1;
    }

    #Get a 'fixed' version of the gene name
    my $mapped_gene_name = $self->fixGeneName('-gene'=>$gene_name, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);

    #If the user defined a gene target list, limit import of genes/transcripts at this point to symbol matches from that list
    if ($self->gene_targets_file){
      next unless ($targets_gene_name->{$gene_name} || $targets_mapped_gene_name->{$gene_name});
    }

    #If the user defined a target_chr or target_chr and target_chr_start-target_chr_end - apply that limit now
    if ($target_chr_start && $target_chr_end){
      next unless (($chr eq $target_chr) && (($chr_start >= $target_chr_start && $chr_start <= $target_chr_end) || ($chr_end >= $target_chr_start && $chr_end <= $target_chr_end) || ($chr_start <= $target_chr_start && $chr_end >= $target_chr_end)));
    }elsif($target_chr && $target_chr ne "ALL"){
      next unless ($chr eq $target_chr);
    }

    #Store Gene level data
    if ($genes{$gid}){ 
      $genes{$gid}{start} = $chr_start if ($chr_start < $genes{$gid}{start});
      $genes{$gid}{end} = $chr_end if ($chr_end > $genes{$gid}{end});
    }else{
      $genes{$gid}{gid} = $gid;
      $genes{$gid}{gene_name} = $gene_name;
      $genes{$gid}{mapped_gene_name} = $mapped_gene_name;
      $genes{$gid}{chr} = $chr;
      $genes{$gid}{start} = $chr_start;
      $genes{$gid}{end} = $chr_end;
    }
    
    #Store Transcript level data
    if ($trans{$tid}){
      $trans{$tid}{start} = $chr_start if ($chr_start < $trans{$tid}{start});
      $trans{$tid}{end} = $chr_end if ($chr_end > $trans{$tid}{end});
    }else{
      $trans{$tid}{gene_name} = $gene_name;
      $trans{$tid}{mapped_gene_name} = $mapped_gene_name;
      $trans{$tid}{gid} = $gid;
      $trans{$tid}{chr} = $chr;
      $trans{$tid}{start} = $chr_start;
      $trans{$tid}{end} = $chr_end;
    }
  }
  close(GTF);

  my $genes_count = keys %genes;
  my $trans_count = keys %trans;
  $self->status_message("Imported $genes_count genes from $gtf_file") if $self->verbose;
  $self->status_message("Imported $trans_count transcripts from $gtf_file") if $self->verbose;

  #Add cytoband annotation for every transcript and gene
  $self->status_message("Annotating genes with cytoband names") if $self->verbose;
  foreach my $gid (keys %genes){
    my $chr = $genes{$gid}{chr};
    my $chr_start = $genes{$gid}{start};
    my $chr_end = $genes{$gid}{end};
    my $cytoband_string = &getCytoband('-ideo_data'=>$ideo_data, '-chr'=>"chr".$chr, '-chr_start'=>$chr_start, '-chr_end'=>$chr_end);
    $genes{$gid}{cytoband} = $cytoband_string;
  }
  $self->status_message("Annotating transcripts with cytoband names") if $self->verbose;
  foreach my $tid (keys %trans){
    my $chr = $trans{$tid}{chr};
    my $chr_start = $trans{$tid}{start};
    my $chr_end = $trans{$tid}{end};
    my $cytoband_string = &getCytoband('-ideo_data'=>$ideo_data, '-chr'=>"chr".$chr, '-chr_start'=>$chr_start, '-chr_end'=>$chr_end);
    $trans{$tid}{cytoband} = $cytoband_string;
  }
  
  my %result;
  $result{genes} = \%genes;
  $result{trans} = \%trans;
  
  return(\%result);
}


#####################################################################################################################################################################
#Load the CNV window coordinates and copy number estimates
#####################################################################################################################################################################
sub loadCnvData{
  my $self = shift;
  my %args = @_;
  my $cnv_file = $args{'-cnv_file'};
  my $target_chr = $args{'-target_chr'};
  my $target_chr_start = $args{'-target_chr_start'};
  my $target_chr_end = $args{'-target_chr_end'};
  my $window_size = $args{'-window_size'};

  my %cnvs;
  open(CNV, "$cnv_file") || die "\n\nCould not open input file: $cnv_file\n\n";
  my $c = 0;
  my $p1 = 0;
  my $p2 = 0;
  while(<CNV>){
    chomp($_);
    if ($_ =~ /^\#|^CHR/){
      next();
    }
    $c++;
    my @line = split("\t", $_);
    my $chr = $line[0];
    my $chr_start = $line[1];
    if ($c == 1){$p1 = $chr_start;}
    if ($c == 2){$p2 = $chr_start;}
    $cnvs{$chr}{$chr_start}{diff} = $line[4];
  }
  close(CNV);
  unless(defined($window_size)){
    $window_size = $p2 - $p1;
  }
  if ($self->verbose){
    $self->status_message("Detected a CNV window size of $window_size bp.  Using this for overlap calculations") if $self->verbose;
  }

  #Now that we know the size, reorganize this object
  my %cnvs2;
  foreach my $chr (keys %cnvs){
    my %chr_cnvs = %{$cnvs{$chr}};
    foreach my $chr_start (sort {$a <=> $b} keys %chr_cnvs){
      my $chr_end = $chr_start+$window_size;
      my $coord = "$chr:$chr_start-$chr_end";

      #If the user defined a target_chr or target_chr and target_chr_start-target_chr_end - apply that limit now
      if ($target_chr_start && $target_chr_end){
        next unless (($chr eq $target_chr) && (($chr_start >= $target_chr_start && $chr_start <= $target_chr_end) || ($chr_end >= $target_chr_start && $chr_end <= $target_chr_end) || ($chr_start <= $target_chr_start && $chr_end >= $target_chr_end)));
      }elsif($target_chr && $target_chr ne "ALL"){
        next unless ($chr eq $target_chr);
      }
      $cnvs2{$coord}{chr} = $chr;
      $cnvs2{$coord}{start} = $chr_start;
      $cnvs2{$coord}{end} = $chr_end;
      $cnvs2{$coord}{diff} = $cnvs{$chr}{$chr_start}{diff};
    }
  }

  return(\%cnvs2);
}


#####################################################################################################################################################################
#Load the CNV segments coordinates and copy number estimates
#####################################################################################################################################################################
sub loadSegmentData{
  my $self = shift;
  my %args = @_;
  my $segments_file = $args{'-segments_file'};
  my $outfile_cnvhmm = $args{'-outfile_cnvhmm'};
  my $target_chr = $args{'-target_chr'};
  my $target_chr_start = $args{'-target_chr_start'};
  my $target_chr_end = $args{'-target_chr_end'};

  #TODO: If the user limited analysis to a single chromosome or piece of a chromosome, only import segments that qualify

  my %segments;
  open(SEGMENTS, "$segments_file") || die "\n\nCould not open segments input file: $segments_file\n\n";
  open(SEGMENTS_CLEAN, ">$outfile_cnvhmm") || die "\n\nCould not open input file: $outfile_cnvhmm\n\n";
  print SEGMENTS_CLEAN "CHR\tSTART\tEND\tSIZE\tnMarkers\tCN1\tAdjusted_CN1\tCN2\tAdjusted_CN2\tLLR_Somatic\tStatus\n";
  while(<SEGMENTS>){
    chomp($_);
    if ($_=~/^\#CHR/ || $_=~/^---/ || $_=~/^Iter\=/ || $_=~/^purity/ || $_=~/CNA predicted/){
      next(); #skip all header lines
    }
    print SEGMENTS_CLEAN "$_\n"; #Create new segments file for use with R script
    my @line = split("\t", $_);
    my $chr = $line[0];
    my $chr_start = $line[1];
    my $chr_end = $line[2];
    my $coord = "$chr:$chr_start-$chr_end";
    my $cn1 = $line[5];
    my $cn1_adj = $line[6];
    my $cn2 = $line[7];
    my $cn2_adj = $line[8];
    my $llr_somatic = $line[9];
    my $status = $line[10];
    my $cn_adj_diff = $cn1_adj - $cn2_adj;


    #If the user defined a target_chr or target_chr and target_chr_start-target_chr_end - apply that limit now
    if ($target_chr_start && $target_chr_end){
      next unless (($chr eq $target_chr) && (($chr_start >= $target_chr_start && $chr_start <= $target_chr_end) || ($chr_end >= $target_chr_start && $chr_end <= $target_chr_end) || ($chr_start <= $target_chr_start && $chr_end >= $target_chr_end)));
    }elsif($target_chr && $target_chr ne "ALL"){
      next unless ($chr eq $target_chr);
    }

    $segments{$coord}{chr} = $chr;
    $segments{$coord}{start} = $chr_start;
    $segments{$coord}{end} = $chr_end;
    $segments{$coord}{cn1} = $cn1;
    $segments{$coord}{cn1_adj} = $cn1_adj;
    $segments{$coord}{cn2} = $cn2;
    $segments{$coord}{cn2_adj} = $cn2_adj;
    $segments{$coord}{llr_somatic} = $llr_somatic;
    $segments{$coord}{status} = $status;
    $segments{$coord}{cn_adj_diff} = $cn_adj_diff;
  }
  close(SEGMENTS);
  close(SEGMENTS_CLEAN);
  return(\%segments);
}


#####################################################################################################################################################################
#For each gene target, get the corresponding transcripts
#Merge these coordinates and get the grand outer coordinates for the gene of interest
#Then get the mean CNV Diff for each gene
#####################################################################################################################################################################
sub calculateMeanCnvDiff{
  my $self = shift;
  my %args = @_;
  my $f_map = $args{'-f_map'};
  my $f_name = $args{'-f_name'};
  my $cnvs = $args{'-cnvs'};
  my $segments = $args{'-segments'};
  my $output_dir = $args{'-output_dir'};

  #Dump bed files for cnvs, segments and features
  my $features_bed_file = $output_dir . "$f_name" . ".bed";
  my $cnvs_bed_file = $output_dir . "cnvs" . ".bed";
  my $segments_bed_file = $output_dir . "segments" . ".bed";
  my $f_count = $self->writeBed('-features'=>$f_map, '-file'=>$features_bed_file);
  my $cnv_count = $self->writeBed('-features'=>$cnvs, '-file'=>$cnvs_bed_file);
  my $segment_count = $self->writeBed('-features'=>$segments, '-file'=>$segments_bed_file);

  #Note, when running bedtools, file B is loaded into memory
  #bedtools intersect command might look something like this:
  #my $bed_cmd3 = "$bedtools_bin_dir"."intersectBed -a $temp_known_acceptors -b $temp_obs_junctions -f 1.0 -s -wa -wb > $result_file"

  #Use intersectBed to get the overlap between the features (genes/transcripts) and cnv windows
  #Parse the overlaps and store for lookup.  Key on feature $coord and store an array of overlapping window/segment $coords;
  my $features_vs_cnv_windows_file = $output_dir . "$f_name" . "_vs_cnv_windows.bed";
  if ($f_count && $cnv_count){
    my $bed_cmd_cnvs = "intersectBed -a $features_bed_file -b $cnvs_bed_file -wa -wb > $features_vs_cnv_windows_file";
    Genome::Sys->shellcmd(cmd => $bed_cmd_cnvs);
  }else{
    Genome::Sys->shelcmd(cmd => "touch $features_vs_cnv_windows_file");
  }
  my $feature_window_overlaps = $self->parseIntersectBed('-file'=>$features_vs_cnv_windows_file);

  #Use intersectBed to get the overlap between the features (genes/transcripts) and cnv segments
  my $features_vs_segments_file = $output_dir . "$f_name" . "_vs_segments.bed";
  if ($f_count && $segment_count){
    my $bed_cmd_segments = "intersectBed -a $features_bed_file -b $segments_bed_file -wa -wb > $features_vs_segments_file";
    Genome::Sys->shellcmd(cmd => $bed_cmd_segments);
  }else{
    Genome::Sys->shellcmd(cmd => "touch $features_vs_segments_file");
  }
  my $feature_segment_overlaps = $self->parseIntersectBed('-file'=>$features_vs_segments_file);

  #Process each 'feature' (gene or transcript) and calculate the CNV status and average CNV value
  my $cnv_diffs_found = 0;
  my $cnv_statuses_found = 0;
  foreach my $fid (sort keys %{$f_map}){
    $f_map->{$fid}->{mean_diff} = 0;
    my $f_chr = $f_map->{$fid}->{chr};
    my $f_start = $f_map->{$fid}->{start};
    my $f_end = $f_map->{$fid}->{end};
    my $f_coord = "$f_chr:$f_start-$f_end";

    my @f_cnv_overlaps;
    @f_cnv_overlaps = @{$feature_window_overlaps->{$f_coord}} if ($feature_window_overlaps->{$f_coord});

    #Calculate the average copy number for the gene of interest
    my $overlaps = scalar(@f_cnv_overlaps);
    my @diffs;
    foreach my $coord (@f_cnv_overlaps){
      push (@diffs, $cnvs->{$coord}->{diff});
    }
    my $sum = 0;
    my $mean_diff = 0;
    foreach my $diff (@diffs){$sum+=$diff;}
    if ($overlaps > 0){
      $mean_diff = $sum/$overlaps;
      $cnv_diffs_found++;
    }
    $mean_diff = sprintf("%.10f", $mean_diff) if ($mean_diff);
    $f_map->{$fid}->{mean_diff} = $mean_diff;

    #Determine the CN status of the gene of interest according to HMM segmentation analysis
    #First, identify the copy number segments that overlap this gene
    $f_map->{$fid}->{cnseg_status} = "NA";
    my @f_segment_overlaps;
    @f_segment_overlaps = @{$feature_segment_overlaps->{$f_coord}} if ($feature_segment_overlaps->{$f_coord});
    
    my $segoverlaps = scalar(@f_segment_overlaps);
    my @statuses;
    foreach my $coord (@f_segment_overlaps){
      push (@statuses, $segments->{$coord}->{status});
    }
    if ($segoverlaps > 0){
      my $gain_sum = 0;
      my $loss_sum = 0;
      foreach my $status (@statuses){
        if ($status eq 'Gain'){$gain_sum++};
        if ($status eq 'Loss'){$loss_sum++};
      }
      if (($gain_sum > 0) && ($gain_sum > $loss_sum)){
        $f_map->{$fid}->{cnseg_status} = "Gain";
        $cnv_statuses_found++;
      }elsif (($loss_sum > 0) && ($loss_sum > $gain_sum)){
        $f_map->{$fid}->{cnseg_status} = "Loss";
        $cnv_statuses_found++;
      }
    }
  }

  #Clean up all the temp files
  Genome::Sys->shellcmd(cmd => "rm -f $features_bed_file $cnvs_bed_file $segments_bed_file $features_vs_cnv_windows_file $features_vs_segments_file");

  $self->status_message("\tFound non-zero CNV diff for $cnv_diffs_found features") if $self->verbose;
  $self->status_message("\tFound cnvhmm gain or loss status for $cnv_statuses_found features") if $self->verbose;

  return();
}


#####################################################################################################################################################################
#Print result to an output file
#####################################################################################################################################################################
sub printCnvResultFile{
  my $self = shift;
  my %args = @_;
  my $outfile_prefix = $args{'-outfile_prefix'};
  my $f_map = $args{'-f_map'};
  my $entrez_ensembl_data = $args{'-entrez_ensembl_data'};
  my $amp_cutoff = $args{'-amp_cutoff'};
  my $del_cutoff = $args{'-del_cutoff'};

  my $outfile = $outfile_prefix . ".tsv";
  my $outfile_amp = $outfile_prefix . ".amp.tsv";
  my $outfile_del = $outfile_prefix . ".del.tsv";
  my $outfile_ampdel = $outfile_prefix . ".ampdel.tsv";
  my @results_files = ($outfile, $outfile_amp, $outfile_del, $outfile_ampdel);

  #Write four files, one with all results, one with amplifications passing a cutoff, one with deletions passing a cutoff, and one with both amplifications and deletions passing cutoffs
  open (OUT, ">$outfile") || die "\n\nCould not open outfile: $outfile for writting\n\n";
  open (OUT_AMP, ">$outfile_amp") || die "\n\nCould not open outfile: $outfile_amp for writting\n\n";
  open (OUT_DEL, ">$outfile_del") || die "\n\nCould not open outfile: $outfile_del for writting\n\n";
  open (OUT_AMPDEL, ">$outfile_ampdel") || die "\n\nCould not open outfile: $outfile_ampdel for writting\n\n";

  #Print headers
  my $header = "fid\tgene_id\tgene_name\tmapped_gene_name\tchr\tstart\tend\tcytoband\tmean_cnv_diff\tcnvhmm_status\n";
  print OUT "$header";
  print OUT_AMP "$header";
  print OUT_DEL "$header";
  print OUT_AMPDEL "$header";

  foreach my $fid (sort {abs($f_map->{$b}->{mean_diff}) <=> abs($f_map->{$a}->{mean_diff})} keys %{$f_map}){
    my $gid = $f_map->{$fid}->{gid};
    my $gene_name = $f_map->{$fid}->{gene_name};
    my $mapped_gene_name = $f_map->{$fid}->{mapped_gene_name};
    my $f_chr = $f_map->{$fid}->{chr};
    my $f_start = $f_map->{$fid}->{start};
    my $f_end = $f_map->{$fid}->{end};
    my $f_cytoband = $f_map->{$fid}->{cytoband};
    my $mean_diff = $f_map->{$fid}->{mean_diff};
    my $cnseg_status = $f_map->{$fid}->{cnseg_status};
    my $string = "$fid\t$gid\t$gene_name\t$mapped_gene_name\t$f_chr\t$f_start\t$f_end\t$f_cytoband\t$mean_diff\t$cnseg_status\n";

    print OUT "$string";
    if ($mean_diff >= $amp_cutoff || $cnseg_status eq "Gain"){
      print OUT_AMP "$string";
      print OUT_AMPDEL "$string";
    }
    if ($mean_diff <= $del_cutoff || $cnseg_status eq "Loss"){
      print OUT_DEL "$string";
      print OUT_AMPDEL "$string";
    }
  }

  close(OUT);
  close(OUT_AMP);
  close(OUT_DEL);
  close(OUT_AMPDEL);

  return(\@results_files);
}


#####################################################################################################################################
#For any features hash with a unique key and 'chr', 'chr_start', and 'chr_end' values write a bed file of the specified name        #
#####################################################################################################################################
sub writeBed{
  my $self = shift;
  my %args = @_;
  my $features = $args{'-features'};
  my $file = $args{'-file'};

  $self->status_message("Writing BED file: $file");
  open (BED, ">$file") || die "\n\nCould not open output BED file: $file\n\n";
  my $f_count = 0;
  foreach my $fid (sort keys %{$features}){
    my $chr = $features->{$fid}->{chr};
    my $chr_start = $features->{$fid}->{start};
    my $chr_end = $features->{$fid}->{end};
    print BED "$chr\t$chr_start\t$chr_end\n";
    $f_count++;
  }
  close (BED);

  return $f_count;
}


#####################################################################################################################################
#Parse overlaps between features (genes/transcripts) and cnv windows or cnv segments
#####################################################################################################################################
sub parseIntersectBed{
  my $self = shift;
  my %args = @_;
  my $file = $args{'-file'};
  my %overlaps;

  open (INTERSECT, "$file") || die "\n\nCould not open BED intersect file: $file\n\n";
  while(<INTERSECT>){
    chomp($_);
    my @line = split("\t", $_);
    my $coord1 = "$line[0]:$line[1]-$line[2]";
    my $coord2 = "$line[3]:$line[4]-$line[5]";
    if ($overlaps{$coord1}){
      push(@{$overlaps{$coord1}}, $coord2);
    }else{
      my @tmp;
      push(@tmp, $coord2);
      $overlaps{$coord1} = \@tmp;
    }
  }
  close(INTERSECT);

  return(\%overlaps);
}

