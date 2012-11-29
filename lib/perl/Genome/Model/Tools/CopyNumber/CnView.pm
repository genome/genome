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
      doc => 'List of gene symbols. These will be highlighted if they meet the minimum CNV cutoffs.  The cancer gene census will be used by default.' },
    name              => { 
      is => 'Text', 
      doc => 'Name for list of gene symbols to be highlighted ("Cancer Genes" will be used by default).' },
    cnv_results_file  => { 
      is => 'Text', 
      doc => 'If you have already generated a CNV results file for a list of genes, you can supply it and skip to the graph generation step' },
    ideogram_file     => { 
      is => 'Text', 
      doc => 'Chromosome ideogram coordinates for the reference build.  Will be automatically selected based on the reference build specified' },
    chr               => { 
      is => 'Text', 
      doc => 'If you wish to limit analysis to a single chromosome, use --chr (e.g. --chr=17). Otherwise all will be processed' },
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
gmt copy-number cn-view --annotation-build=124434505 --cnv-file=/gscmnt/gc13001/info/model_data/2888915570/build129973671/variants/cnvs.hq --segments-file=/gscmnt/gc2013/info/model_data/2889110844/build130030495/PNC6/clonality/cnaseq.cnvhmm --output-dir=/tmp/ --gene-targets-file=/gscmnt/sata132/techd/mgriffit/reference_annotations/GeneSymbolLists/CancerGeneCensusPlus_Sanger.txt --name='CancerGeneCensusPlus_Sanger' 
EOS
}

sub help_detail {
  return <<EOS
Takes standard CNV analysis output, together with CNAseg output to produce chromosome-by-chromosome and gene level colored copy number plots and summary files.

EOS
}

sub execute{
  my $self = shift;

  #WTF does this mean? Put a comment here...
$DB::single = 1;

  #Required
  my $annotation_build = $self->annotation_build;
  my $cnv_file = $self->cnv_file;
  my $segments_file = $self->segments_file;
  my $somatic_build = $self->somatic_build;
  my $output_dir = $self->output_dir;

  #Optional
  my $gene_targets_file = $self->gene_targets_file;
  my $name = $self->name;
  my $cnv_results_file = $self->cnv_results_file;
  my $ideogram_file = $self->ideogram_file;
  my $verbose = $self->verbose;
  my $chr = $self->chr;
  my $chr_start = $self->chr_start;
  my $chr_end = $self->chr_end;
  my $image_type = $self->image_type;

  #Set amplification / deletion cutoffs that will be used to created filtered files of results for convenience in downstream analyses
  #Try two fold gain and 1/2 loss to start
  my $amp_cutoff = 2;
  my $del_cutoff = -0.5;

  #Check/format all input parameters and options
  #Based on the reference build specified, define paths to transcript/gene id mapping files
  my ($transcript_info_file, $subdir, $outfile, $outfile_amp, $outfile_del, $outfile_ampdel, $outfile_cnvhmm, $name_f);
  my $extra_flank = 100000;

  ####################################################################################################################################
  #Check inputs
  my $inputs_are_good = eval { 
    unless ($cnv_file || $somatic_build){
      $self->usage_message($self->help_usage_complete_text);
      $self->error_message("Must supply either cnv_file or somatic_build (Somatic Variation)");
      return;
    }
    my $clinseq_annotations_dir="/gscmnt/sata132/techd/mgriffit/reference_annotations/";
    my $annotation_data_dir=$annotation_build->data_directory;
    
    my $reference_sequence_build=$annotation_build->reference_sequence;
    my $default_build37 = Genome::Model::Build->get(102671028);
    my $default_build36 = Genome::Model::Build->get(101947881);
    if ($reference_sequence_build->is_compatible_with($default_build37)){
      $ideogram_file ||= $clinseq_annotations_dir . "hg19/ideogram/ChrBandIdeogram.tsv";
    }elsif ($reference_sequence_build->is_compatible_with($default_build36)){
      $ideogram_file ||= $clinseq_annotations_dir . "hg18/ideogram/ChrBandIdeogram.tsv";
    }else {
      $self->error_message("Specified reference build resolved from annotation build is not compatible with default build36 or build37");
      return;
    }

    $transcript_info_file = $annotation_data_dir . "/annotation_data/transcripts.csv";
    unless (-e $transcript_info_file && -e $ideogram_file){
      $self->error_message("One or more of the following annotation files is missing:\n$transcript_info_file\n$ideogram_file");
      return;
    }

    # TODO: switch that name and version (date) for a value in the processing profile if this changes
    my $gene_symbol_lists_dir = Genome::Sys->dbpath("tgi-gene-symbol-lists","2012-11-13");

    #Set the default gene targets file if it wasnt specified by the user
    unless ($gene_targets_file){
      $gene_targets_file = $gene_symbol_lists_dir . "/CancerGeneCensusPlus_Sanger.txt"
    }
    unless (-e $gene_targets_file){
      $self->error_message("Gene targets file not found: $gene_targets_file");
      return;
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

    unless (-e $output_dir && -d $output_dir){
      $self->status_message("Creating dir: $output_dir");
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
      $name = "Cancer Genes";
      $name_f = $name;
      $name_f =~ s/ //g;
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
    if ($chr eq "ALL"){
      $outfile = "$subdir"."CNView_"."$name_f".".tsv";
      $outfile_amp = "$subdir"."CNView_"."$name_f".".amp.tsv";
      $outfile_del = "$subdir"."CNView_"."$name_f".".del.tsv";
      $outfile_ampdel = "$subdir"."CNView_"."$name_f".".ampdel.tsv";
      $outfile_cnvhmm = "$subdir"."cnaseq.cnvhmm.tsv";
    }elsif ($chr_start && $chr_end){
      $outfile = "$subdir"."CNView_"."$name_f"."_chr"."$chr"."_"."$chr_start-$chr_end".".tsv";
      $outfile_amp = "$subdir"."CNView_"."$name_f"."_chr"."$chr"."_"."$chr_start-$chr_end".".amp.tsv";
      $outfile_del = "$subdir"."CNView_"."$name_f"."_chr"."$chr"."_"."$chr_start-$chr_end".".del.tsv";
      $outfile_ampdel = "$subdir"."CNView_"."$name_f"."_chr"."$chr"."_"."$chr_start-$chr_end".".ampdel.tsv";
      $outfile_cnvhmm = "$subdir"."cnaseq.cnvhmm.tsv";
    }else{
      $outfile = "$subdir"."CNView_"."$name_f"."_chr"."$chr".".tsv";
      $outfile_amp = "$subdir"."CNView_"."$name_f"."_chr"."$chr".".amp.tsv";
      $outfile_del = "$subdir"."CNView_"."$name_f"."_chr"."$chr".".del.tsv";
      $outfile_ampdel = "$subdir"."CNView_"."$name_f"."_chr"."$chr".".ampdel.tsv";
      $outfile_cnvhmm = "$subdir"."cnaseq.cnvhmm.tsv";
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

  unless ($inputs_are_good) {
    $self->error_message("error processing inputs!");
    return;
  }

  #Done checking inputs
  ####################################################################################################################################

  my $entrez_ensembl_data = &loadEntrezEnsemblData();

  #If the user is supplying a pre-computed CNV file, the following steps will be skipped
  my ($gt_map, $t_coords, $window_size, $cnvs, $segments, $targets);
  if ($cnv_results_file){
    if($verbose){ print BLUE, "\n\nUsing user supplied CNV results file: $cnv_results_file", RESET; }
  }else{
    if($verbose){ print BLUE, "\n\nGenerating CNV results file: $cnv_results_file", RESET; }

    #Load the gene to transcript name mappings.
    #Key on unique gene ID - reference to a hash containing all transcript IDs associated with that gene
    $gt_map = $self->loadGtMap('-transcript_info_file'=>$transcript_info_file);

    #Load the transcripts outer coords
    $t_coords = $self->loadTranscriptCoords('-transcript_info_file'=>$transcript_info_file);

    #Load the CNV window coordinates and copy number estimates
    $cnvs = $self->loadCnvData('-cnv_file'=>$cnv_file, '-window_size'=>\$window_size);

    #Load the CNV segments and copy number estimates
    $segments = $self->loadSegmentData('-segments_file'=>$segments_file, '-outfile_cnvhmm'=>$outfile_cnvhmm);

    #Load the gene targets of interest
    $targets = $self->loadTargetGenes('-gene_targets_file'=>$gene_targets_file);

    #For each gene target, get the corresponding transcripts
    #Merge these coordinates and get the grand outer coordinates for the gene of interest
    #Then get the mean CNV Diff for each gene
    $self->calculateMeanCnvDiff('-targets'=>$targets, '-gt_map'=>$gt_map, '-t_coords'=>$t_coords, '-cnvs'=>$cnvs, '-window_size'=>$window_size, '-segments'=>$segments);

    #Print result to an output file
    $self->printCnvResultFile('-outfile'=>$outfile, '-outfile_amp'=>$outfile_amp, '-outfile_del'=>$outfile_del, '-outfile_ampdel'=>$outfile_ampdel, '-targets'=>$targets, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-amp_cutoff'=>$amp_cutoff, '-del_cutoff'=>$del_cutoff, '-chr'=>$chr, '-chr_start'=>$chr_start, '-chr_end'=>$chr_end, '-cnv_results_file'=>\$cnv_results_file);
  }

  #Assume that the CNview.R script resides in the same dir as this script and obtain the path automatically

  #Execute the R script that generates the CNV plots
  my $chr_name=$chr; #R script still expects chr number to be prefixed with "chr"
  unless ($chr_name eq 'ALL'){$chr_name="chr".$chr;}
  my $rscript = __FILE__.".R";
  my $r_cmd = "$rscript '$name' $cnv_file $cnv_results_file $outfile_cnvhmm $ideogram_file $subdir $chr_name $chr_start $chr_end $image_type";
  my $r_cmd_stdout = "$subdir"."CNView.R.stdout";
  my $r_cmd_stderr = "$subdir"."CNView.R.stderr";
  if ($verbose){
    print BLUE, "\n\nExecuting R code:\n$r_cmd\n", RESET;
  }else{
    $r_cmd .= " 1>$r_cmd_stdout 2>$r_cmd_stderr";
  }
  Genome::Sys->shellcmd(cmd => $r_cmd);

  #Annotate all gene entries in the output file with a cytoband label using the coordinates in the ideogram data file, and the coordinates of the gene

  #Parse import the coordinates of the ideogram file using a subroutine
  my $ideo_data = &importIdeogramData('-ideogram_file'=>$ideogram_file);

  #Modify the output files ($cnv_results_file) above to annotate each gene with the corresponding cytoband(s) spanned by the gene
  #Load file, get: (chr, start, end), determine cytobands that overlap, update output file, replace old file with new
  #e.g. of cytoband: 11q13.2  OR  11q13.2 - 11q13.3
  my @results_files = ($outfile, $outfile_amp, $outfile_del, $outfile_ampdel);
  foreach my $results_file (@results_files){
    if (-e $results_file){
      #Parse input file and determine cytobands for each gene
      my %data;
      my $new_cnv_results_file = "$results_file".".tmp";
      my $header_line = '';
      my $header = 1;
      my $l = 0;
      open (CNV, "$results_file") || die "\n\nCould not open CNV results file: $results_file\n\n";
      my %columns;
      while(<CNV>){
        $l++;
        chomp($_);
        my @line = split("\t", $_);
        if ($header){
          $header = 0;
          $header_line = $_;
          my $p = 0;
          foreach my $col (@line){
            $columns{$col}{position} = $p;
            $p++;
          }
          next();
        }
        $data{$l}{line} = $_;
        my $chr = $line[$columns{'Chr'}{position}];
        my $chr_start = $line[$columns{'Start'}{position}];
        my $chr_end = $line[$columns{'End'}{position}];
        my $cytoband_string = &getCytoband('-ideo_data'=>$ideo_data, '-chr'=>"chr".$chr, '-chr_start'=>$chr_start, '-chr_end'=>$chr_end);
        $data{$l}{cytoband} = $cytoband_string;
      }
      close (CNV);

      #Print out new file containing the cytoband info
      open(CNV_NEW, ">$new_cnv_results_file") || die "\n\nCould not open new CNV results file: $new_cnv_results_file\n\n";
      print CNV_NEW "$header_line\tCytoband\n";
      foreach my $l (sort {$a <=> $b} keys %data){
        print CNV_NEW "$data{$l}{line}\t$data{$l}{cytoband}\n";
      }
      close(CNV_NEW);

      #Overwrite the old file with the new file
      my $mv_cmd = "mv $new_cnv_results_file $results_file";
      Genome::Sys->shellcmd(cmd => $mv_cmd);

    }else{
      print YELLOW, "\n\nCNV results file not found - not possible to annotate with Cytoband info\n\n", RESET;
    }
  }

  if ($verbose){ print BLUE, "\n\nResults written to:\n$subdir\n\n", RESET; }

  return(1);
}


#TODO: Flip the logic of this whole script around so that the primary data object is the ensembl genes from the annotation object
#TODO: If the user supplies a gene symbol list, filter this list down so that only matching symbols are allowed
#TODO: If the user does not supply such as list, process all genes.
#TODO: The final output file should be the same except it should contain the Ensembl gene ID as well as the symbol


#####################################################################################################################################################################
#Load the gene to transcript name mappings.
#Key on unique gene ID - reference to a hash containing all transcript IDs associated with that gene
#####################################################################################################################################################################
sub loadGtMap{
  my $self = shift;
  my %args = @_;
  my $transcript_info_file = $args{'-transcript_info_file'};

  my %gt_map;
  open (GT, "$transcript_info_file") || die "\n\nCould not open input file: $transcript_info_file\n\n";
  while(<GT>){
    chomp($_);
    my @line = split(",", $_);
    my $tid = $line[4];
    my $gene = $line[11];
    unless ($gene){
      next();
    }
    if ($gt_map{$gene}){ 
      my $trans_ref = $gt_map{$gene}{trans};
      $trans_ref->{$tid}=1;
    }else{
      my %trans;
      $trans{$tid}=1;
      $gt_map{$gene}{trans} = \%trans;
    }
  }
  close(GT);
  return(\%gt_map);
}


#####################################################################################################################################################################
#Load the transcripts outer coords
#####################################################################################################################################################################
sub loadTranscriptCoords{
  my $self = shift;
  my %args = @_;
  my $transcript_info_file = $args{'-transcript_info_file'};
  my %t_coords;
  
  open (TC, "$transcript_info_file") || die "\n\nCould not open input file: $transcript_info_file\n\n";
  while(<TC>){
    chomp($_);
    my @line = split(",", $_);
    my $tid = $line[4];
    my $chr_input = $line[7];

    $t_coords{$tid}{chr} = $chr_input;
    $t_coords{$tid}{start} = $line[2];
    $t_coords{$tid}{end} = $line[3];
  }
  close(TC);
  return(\%t_coords);
}


#####################################################################################################################################################################
#Load the CNV window coordinates and copy number estimates
#####################################################################################################################################################################
sub loadCnvData{
  my $self = shift;
  my %args = @_;
  my $cnv_file = $args{'-cnv_file'};
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
    my $chr_input = $line[0];
    my $start = $line[1];
    if ($c == 1){$p1 = $start;}
    if ($c == 2){$p2 = $start;}
    $cnvs{$chr_input}{$start}{tumor} = $line[2];
    $cnvs{$chr_input}{$start}{normal} = $line[3];
    $cnvs{$chr_input}{$start}{diff} = $line[4];
  }
  close(CNV);
  ${$window_size} = $p2 - $p1;
  if ($self->verbose){
    $self->status_message("Detected a CNV window size of $window_size bp.  Using this for overlap calculations");
  }
  return(\%cnvs);
}


#####################################################################################################################################################################
#Load the CNV segments coordinates and copy number estimates
#####################################################################################################################################################################
sub loadSegmentData{
  my $self = shift;
  my %args = @_;
  my $segments_file = $args{'-segments_file'};
  my $outfile_cnvhmm = $args{'-outfile_cnvhmm'};

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
    my $chr_input = "$line[0]";
    my $start = $line[1];
    my $end = $line[2];
    my $cn1 = $line[5];
    my $cn1_adj = $line[6];
    my $cn2 = $line[7];
    my $cn2_adj = $line[8];
    my $llr_somatic = $line[9];
    my $status = $line[10];
    my $cn_adj_diff = $cn1_adj - $cn2_adj;
    $segments{$chr_input}{$start}{$end}{cn1} = $cn1;
    $segments{$chr_input}{$start}{$end}{cn1_adj} = $cn1_adj;
    $segments{$chr_input}{$start}{$end}{cn2} = $cn2;
    $segments{$chr_input}{$start}{$end}{cn2_adj} = $cn2_adj;
    $segments{$chr_input}{$start}{$end}{llr_somatic} = $llr_somatic;
    $segments{$chr_input}{$start}{$end}{status} = $status;
    $segments{$chr_input}{$start}{$end}{cn_adj_diff} = $cn_adj_diff;
  }
  close(SEGMENTS);
  close(SEGMENTS_CLEAN);
  return(\%segments);
}


#####################################################################################################################################################################
#Load the gene targets of interest
#####################################################################################################################################################################
sub loadTargetGenes{
  my $self = shift;
  my %args = @_;
  my $gene_targets_file = $args{'-gene_targets_file'}; 

  my %targets;
  open (TARGET ,"$gene_targets_file") || die "\n\nCould not open gene targets file\n\n";
  while(<TARGET>){
    chomp($_);
    if ($_ =~ /(\S+)/){
      #Skip "n/a"
      if ($1 eq "n/a"){
        next();
      }else{
        $targets{$1}{found} = 0;
      }
    }
  }
  close(TARGET);
  return(\%targets);
}


#####################################################################################################################################################################
#For each gene target, get the corresponding transcripts
#Merge these coordinates and get the grand outer coordinates for the gene of interest
#Then get the mean CNV Diff for each gene
#####################################################################################################################################################################
sub calculateMeanCnvDiff{
  my $self = shift;
  my %args = @_;
  my $targets = $args{'-targets'};
  my $gt_map = $args{'-gt_map'};
  my $t_coords = $args{'-t_coords'};
  my $cnvs = $args{'-cnvs'};
  my $window_size = $args{'-window_size'};
  my $segments = $args{'-segments'};

  #TODO: This subroutine should be going through the ensembl gene IDs (*not* the target symbols)

  foreach my $target (sort keys %{$targets}){
    if ($self->verbose){
      print BLUE, "\n\n$target", RESET;
    }
    $targets->{$target}->{mean_diff} = 0;
    if ($gt_map->{$target}){
      #NOTE: A target is only found if the NAME matches exactly.  This should be improved...
      #Get the transcript list of this gene
      my $trans_ref = $gt_map->{$target}->{trans};
      my $tcount = keys %{$trans_ref};
      if ($self->verbose){
        print BLUE, "\n\tFound $tcount transcripts", RESET;
      }

      #Merge the transcript coordinates.  Make sure all coordinates are from the same chromosome
      my %chr_list;
      my $grand_chr = '';
      my $grand_start = 1000000000000000000000000000000000000;
      my $grand_end = 0;

      #Check chromosomes of transcripts.  If there is more than one, chose the one with the most support
      foreach my $tid (keys %{$trans_ref}){
        my $chr = $t_coords->{$tid}->{chr};
        #Ignore transcripts from the haplotype chromosomes - Add to this list as they come up...
        if ($chr =~ /chr4_ctg9_hap1|chr6_apd_hap1|chr6_mann_hap4|chr6_qbl_hap6|chr6_ssto_hap7|chr6_mcf_hap5|chr6_cox_hap2|chr6_dbb_hap3|chr17_ctg5_hap1/){
          next();        
        }
        $chr_list{$chr}{count}++;
      }
      my $i = 0;
      foreach my $chr (sort {$chr_list{$b}{count} <=> $chr_list{$a}{count}} keys %chr_list){
        $i++;
        if ($i == 1){
          $grand_chr = $chr;
        }
      }

      #Now using the chr identified above, get the grand coords
      foreach my $tid (keys %{$trans_ref}){
        my $chr = $t_coords->{$tid}->{chr};
        unless ($chr eq $grand_chr){
          next();
        }
        my $start = $t_coords->{$tid}->{start};
        my $end = $t_coords->{$tid}->{end};
        if ($start < $grand_start){$grand_start = $start;}
        if ($end > $grand_end){$grand_end = $end;}
      }
      my $chr_count = keys %chr_list;
      my $grand_size = $grand_end - $grand_start;
      if ($chr_count == 1){
        if ($self->verbose){
          $self->status_message("\t$grand_chr:$grand_start-$grand_end ($grand_size)");
        }
      }else{
        my @chrs = keys %chr_list;
        if ($self->verbose){
          $self->warning_message("\tDid not find a single chromosome for this gene!  (@chrs) - chose the one with the most support (i.e. most transcripts)");
          $self->status_message("\t$grand_chr:$grand_start-$grand_end ($grand_size)");
        }
      }
      #Make sure CNV value exists for this chromosome
      unless ($cnvs->{$grand_chr}){
        if ($self->verbose){
          $self->warning_message("\tNo CNV data defined for this chromosome: $grand_chr");
        }
        next();
      }

      #Now identify the copy number windows that overlap this gene
      my %chr_cnvs = %{$cnvs->{$grand_chr}};
      my $window_count = keys %chr_cnvs;
      if ($self->verbose){
        $self->status_message("\tFound $window_count windows for this chromosome");
      }

      #Store values for later
      $targets->{$target}->{found} = 1;
      $targets->{$target}->{chr} = $grand_chr;
      $targets->{$target}->{start}= $grand_start;
      $targets->{$target}->{end} = $grand_end;
      $targets->{$target}->{size} = $grand_size;

      #Calculate the average copy number for the gene of interest
      my $overlaps = 0;
      my @diffs;
      foreach my $pos (sort {$a <=> $b} keys %chr_cnvs){
        my $cnv_start = $pos;
        my $cnv_end = $pos+$window_size;
        if (($cnv_start >= $grand_start && $cnv_start <= $grand_end) || ($cnv_end >= $grand_start && $cnv_end <= $grand_end) || ($cnv_start <= $grand_start && $cnv_end >= $grand_end)){
          $overlaps++;
          push (@diffs, $chr_cnvs{$pos}{diff});
        }
      }
      my $sum = 0;
      my $mean_diff = 0;
      foreach my $diff (@diffs){$sum+=$diff;}
      if ($overlaps > 0){
        $mean_diff = $sum/$overlaps;
      }
      if ($self->verbose){
        $self->status_message("\tFound $overlaps overlapping CNV windows with mean diff: $mean_diff");
      }
      $targets->{$target}->{mean_diff} = $mean_diff;

      #Determine the CN status of the gene of interest according to HMM segmentation analysis
      #First, identify the copy number segments that overlap this gene
      $targets->{$target}->{cnseg_status} = "NA";
      if ($segments->{$grand_chr}){
        my %chr_segments = %{$segments->{$grand_chr}};
        my $segment_count = keys %chr_cnvs;
        if ($self->verbose){
          $self->status_message("\n\tFound $segment_count segments for this chromosome");
        }
        #Determine the CN status for the gene of interest
        my $segoverlaps = 0;
        my @statuses;
        foreach my $seg_start (sort {$a <=> $b} keys %chr_segments){
          foreach my $seg_end (sort {$a <=> $b} keys %{$chr_segments{$seg_start}}){
            if (($seg_start >= $grand_start && $seg_start <= $grand_end) || ($seg_end >= $grand_start && $seg_end <= $grand_end) || ($seg_start <= $grand_start && $seg_end >= $grand_end)){
            $segoverlaps++;
            push (@statuses, $chr_segments{$seg_start}{$seg_end}{status});
            }
          }
        }
        if ($segoverlaps > 0){
          my $gain_sum = 0;
          my $loss_sum = 0;
          foreach my $status (@statuses){
            if ($status eq 'Gain'){$gain_sum++};
            if ($status eq 'Loss'){$loss_sum++};
          }
          if ($gain_sum>0 && $gain_sum>$loss_sum){
            $targets->{$target}->{cnseg_status} = "Gain";
          }
          if ($loss_sum>0 && $loss_sum>$gain_sum){
            $targets->{$target}->{cnseg_status} = "Loss";
          }
        }
        if ($self->verbose){
          $self->status_message("\tFound $segoverlaps overlapping segments with status: $targets->{$target}->{cnseg_status}");
        }
      }
    }else{
      if ($self->verbose){
        $self->status_message("\tCould not find any transcripts for this gene: $target");
      }
    }
  }
  return();
}


#####################################################################################################################################################################
#Print result to an output file
#####################################################################################################################################################################
sub printCnvResultFile{
  my $self = shift;
  my %args = @_;
  my $outfile = $args{'-outfile'};
  my $outfile_amp = $args{'-outfile_amp'};
  my $outfile_del = $args{'-outfile_del'};
  my $outfile_ampdel = $args{'-outfile_ampdel'};
  my $targets = $args{'-targets'};
  my $entrez_ensembl_data = $args{'-entrez_ensembl_data'};
  my $amp_cutoff = $args{'-amp_cutoff'};
  my $del_cutoff = $args{'-del_cutoff'};
  my $chr = $args{'-chr'};
  my $chr_start = $args{'-chr_start'};
  my $chr_end = $args{'-chr_end'};
  my $cnv_results_file = $args{'-cnv_results_file'};

  #Write four files, one with all results, one with amplifications passing a cutoff, one with deletions passing a cutoff, and one with both amplifications and deletions passing cutoffs
  open (OUT, ">$outfile") || die "\n\nCould not open outfile: $outfile for writting\n\n";
  open (OUT_AMP, ">$outfile_amp") || die "\n\nCould not open outfile: $outfile_amp for writting\n\n";
  open (OUT_DEL, ">$outfile_del") || die "\n\nCould not open outfile: $outfile_del for writting\n\n";
  open (OUT_AMPDEL, ">$outfile_ampdel") || die "\n\nCould not open outfile: $outfile_ampdel for writting\n\n";
  print OUT "Symbol\tmapped_gene_name\tChr\tStart\tEnd\tMean CNV Diff\tCNVhmm Status\n";
  print OUT_AMP "Symbol\tmapped_gene_name\tChr\tStart\tEnd\tMean CNV Diff\tCNVhmm Status\n";
  print OUT_DEL "Symbol\tmapped_gene_name\tChr\tStart\tEnd\tMean CNV Diff\tCNVhmm Status\n";
  print OUT_AMPDEL "Symbol\tmapped_gene_name\tChr\tStart\tEnd\tMean CNV Diff\tCNVhmm Status\n";

  foreach my $target (sort {abs($targets->{$b}->{mean_diff}) <=> abs($targets->{$a}->{mean_diff})} keys %{$targets}){
    if ($targets->{$target}->{found}){
      my $fixed_gene_name = &fixGeneName('-gene'=>$target, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>$self->verbose);
      my $target_chr = $targets->{$target}->{chr};
      my $target_chr_start = $targets->{$target}->{start};
      my $target_chr_end = $targets->{$target}->{end};
      my $mean_diff = $targets->{$target}->{mean_diff};
      my $cnseg_status = $targets->{$target}->{cnseg_status};
      my $string = "$target\t$fixed_gene_name\t$target_chr\t$target_chr_start\t$target_chr_end\t$mean_diff\t$cnseg_status\n";
      #Name the output file
      if ($chr eq "ALL"){
        print OUT "$string";
        if ($mean_diff >= $amp_cutoff || $cnseg_status eq "Gain"){
          print OUT_AMP "$string";
          print OUT_AMPDEL "$string";
        }
        if ($mean_diff <= $del_cutoff || $cnseg_status eq "Loss"){
          print OUT_DEL "$string";
          print OUT_AMPDEL "$string";
        }
      }elsif (($chr =~ /chr/i) && $chr_start && $chr_end){
        unless (($target_chr eq $chr) && (($target_chr_start >= $chr_start && $target_chr_start <= $chr_end) || ($target_chr_end >= $chr_start && $target_chr_end <= $chr_end) || ($target_chr_start <= $chr_start && $target_chr_end >= $chr_end))){
          next();
        }
        print OUT "$string";
        if ($mean_diff >= $amp_cutoff || $cnseg_status eq "Gain"){
          print OUT_AMP "$string";
          print OUT_AMPDEL "$string";
        }
        if ($mean_diff <= $del_cutoff || $cnseg_status eq "Loss"){
          print OUT_DEL "$string";
          print OUT_AMPDEL "$string";
        }

      }else{
        unless ($target_chr eq $chr){
          next();
        }
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
    }
  }
  ${$cnv_results_file} = $outfile;

  close(OUT);
  close(OUT_AMP);
  close(OUT_DEL);
  close(OUT_AMPDEL);

  return();
}


