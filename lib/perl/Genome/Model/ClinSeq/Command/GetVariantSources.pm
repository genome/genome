package Genome::Model::ClinSeq::Command::GetVariantSources;

#Written by Obi Griffith

use strict;
use warnings;
use Genome;
use Data::Dumper;

class Genome::Model::ClinSeq::Command::GetVariantSources {
    is => 'Command::V2',
    has_input => [
        builds => { 
              is => 'Genome::Model::Build::SomaticVariation',
              is_many => 1,
              shell_args_position => 1,
              require_user_verify => 0,
              doc => 'somatic variation build(s) to get variant sources from',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
    ],
    has_output => [
        indel_variant_sources_file => {
              is => 'FilesystemPath',
              is_optional =>1,
        },
        snv_variant_sources_file => {
              is => 'FilesystemPath',
              is_optional =>1,
        },
    ],
    doc => 'summarize the sources of variants (i.e., which variant callers) for a somatic variation build',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq get-variant-sources --outdir=/tmp/  128884819

genome model clin-seq get-variant-sources --outdir=/tmp/  id=128884819

genome model clin-seq get-variant-sources --outdir=/tmp/  model.id=2888329352

genome model clin-seq get-variant-sources --outdir=/tmp/  "model.name='H_JG-300000-1206887.somatic_variation-1'"

genome model clin-seq get-variant-sources --outdir=/tmp/  'id in [128884819,128884852]'

EOS
}

sub help_detail {
    return <<EOS
Summarize source of variants (i.e., snv/indel caller) for one or more somatic variation builds 

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
  return @errors;
}

#Global variables - #figure out a way to do this with a method
my %indel_caller;
my %snv_caller;

sub execute {
  my $self = shift;
  my @builds = $self->builds;
  my $outdir = $self->outdir;

  unless ($outdir =~ /\/$/){
    $outdir .= "/";
  }

  my $somatic_build_count = scalar(@builds);
  for my $somatic_build (@builds) {
    #If there is more than one somatic variation build supplied... create sub-directories for each
    my $build_outdir;
    if ($somatic_build_count > 1){
      $build_outdir = $outdir . $somatic_build->id . "/";
      mkdir ($build_outdir);
    }else{
      $build_outdir = $outdir;
    }
    my $build_dir = $somatic_build->data_directory;

    #Set files for output
    my $indel_outfile = $build_outdir . "indel_sources.tsv";
    my $snv_outfile = $build_outdir . "snv_sources.tsv";

    #Set a list of variant files to consider
    my %indel_files;
    $indel_files{1}{file} = $build_dir . "/effects/indels.hq.novel.tier1.v2.bed";
    $indel_files{1}{tier} = "tier1";
    $indel_files{2}{file} = $build_dir . "/effects/indels.hq.novel.tier2.v2.bed";
    $indel_files{2}{tier} = "tier2";
    $indel_files{3}{file} = $build_dir . "/effects/indels.hq.novel.tier3.v2.bed";
    $indel_files{3}{tier} = "tier3";
    $indel_files{4}{file} = $build_dir . "/effects/indels.hq.novel.tier4.v2.bed";
    $indel_files{4}{tier} = "tier4";

    my %snv_files;
    $snv_files{1}{file} = $build_dir . "/effects/snvs.hq.novel.tier1.v2.bed";
    $snv_files{1}{tier} = "tier1";
    $snv_files{2}{file} = $build_dir . "/effects/snvs.hq.novel.tier2.v2.bed";
    $snv_files{2}{tier} = "tier2";
    $snv_files{3}{file} = $build_dir . "/effects/snvs.hq.novel.tier3.v2.bed";
    $snv_files{3}{tier} = "tier3";
    $snv_files{4}{file} = $build_dir . "/effects/snvs.hq.novel.tier4.v2.bed";
    $snv_files{4}{tier} = "tier4";

    #Locate the final indel/snv results files and load into memory
    #For indels, use ~/effects/indels.hq.novel.tier1.v2.bed ?  (Or the annotated file?)
    #For SNVs, use ~/effects/snvs.hq.novel.tier1.v2.bed
    my $indel_results_file = $build_outdir . "indels.hq.novel.tier1-3.v2.bed";
    open (INDELS_OUT, ">$indel_results_file") || die $self->error_message("Could not open output file: $indel_results_file");
    my %indels;
    my $l=0;
    foreach my $c (sort {$a <=> $b} keys %indel_files){
      my $file = $indel_files{$c}{file};
      my $tier = $indel_files{$c}{tier};
      open (INDELS, $file) or die "can't open $file\n";
      while(<INDELS>){
        $l++;
        chomp;
        my $line=$_;
        my @data=split("\t",$_);
        my $variant_string="$data[0]".":"."$data[1]"."-"."$data[2]"." ($data[3])";
        my $coord_string="$data[0]".":"."$data[1]"."-"."$data[2]";
        $indels{$l}{line}=$line;
        $indels{$l}{variant_string}=$variant_string;
        $indels{$l}{coord_string}=$coord_string;
        $indels{$l}{tier} = $tier;
        print INDELS_OUT "$line\n";
      }
      close(INDELS);
    }
    close (INDELS_OUT);

    my $snv_results_file = $build_outdir . "snvs.hq.novel.tier1-3.v2.bed";     
    open (SNVS_OUT, ">$snv_results_file") || die $self->error_message("Could not open output file: $snv_results_file");
    my %snvs;
    $l=0;
    foreach my $c (sort {$a <=> $b} keys %snv_files){
      my $file = $snv_files{$c}{file};
      my $tier = $snv_files{$c}{tier};
      open (SNVS, $file) or die "can't open $file\n";
      while(<SNVS>){
        $l++;
        chomp;
        my $line=$_;
        my @data=split("\t",$_);
        my $variant_string="$data[0]".":"."$data[1]"."-"."$data[2]"." ($data[3])";
        my $coord_string="$data[0]".":"."$data[1]"."-"."$data[2]";
        $snvs{$l}{line}=$line;
        $snvs{$l}{variant_string}=$variant_string;
        $snvs{$l}{coord_string}=$coord_string;
        $snvs{$l}{tier} = $tier;
        print SNVS_OUT "$line\n";
      }
      close(SNVS);
    }
    close (SNVS_OUT);

    #Sort the BED files using joinx
    my $indel_results_file_sorted = $indel_results_file . ".sort";
    my $joinx_indel_sort_cmd = Genome::Model::Tools::Joinx::Sort->create(output_file=>$indel_results_file_sorted, input_files=>[$indel_results_file]);
    $joinx_indel_sort_cmd->execute();
    unlink $indel_results_file;
    Genome::Sys->move_file($indel_results_file_sorted, $indel_results_file);
    
    my $snv_results_file_sorted = $snv_results_file . ".sort";
    my $joinx_snv_sort_cmd = Genome::Model::Tools::Joinx::Sort->create(output_file=>$snv_results_file_sorted, input_files=>[$snv_results_file]);
    $joinx_snv_sort_cmd->execute();
    unlink $snv_results_file;
    Genome::Sys->move_file($snv_results_file_sorted, $snv_results_file);

    my $indel_count = keys %indels;
    $self->status_message("Stored $indel_count indels");
    my $snv_count = keys %snvs;
    $self->status_message("Stored $snv_count indels");

    #Locate the individual indel/snv files for each caller to use in joinx intersect
    #This should be replaced by a method which somehow determines the appropriate files automatically
    #Depending on whether the somatic variation build is for exome or wgs data, the paths will differ - this should also be determined automatically

    my ($indel_strelka_results_file, $indel_gatk_results_file, $indel_pindel_results_file, $indel_varscan_results_file);

    #Create a list of possible indels file paths
    my @strelka_indel_paths = ("$build_dir/variants/indel/strelka-0.4.6.2-bc1213eb5850cc1810af3214c95cf30d/indels.hq.bed",
                               "$build_dir/variants/indel/strelka-0.4.6.2-c5009c08801c3ffa834ecb28d4293d27/indels.hq.bed",
                               "$build_dir/variants/indel/strelka-0.4.6.2-14acc00d0b01975892118ec71cfc3506/indels.hq.bed",
                               "$build_dir/variants/indel/strelka-0.4.6.2-673995e8237c2c733def86d8d9b3d5a6/indels.hq.bed");
    $indel_strelka_results_file = $self->checkResultFile('-paths'=>\@strelka_indel_paths, '-caller'=>"strelka");

    my @gatk_indel_paths = ("$build_dir/variants/indel/gatk-somatic-indel-5336-d41d8cd98f00b204e9800998ecf8427e/false-indel-v1-05fbf69c10534fd630b99e44ddf73c7f/indels.hq.bed");
    $indel_gatk_results_file = $self->checkResultFile('-paths'=>\@gatk_indel_paths, '-caller'=>"gatk");
    
    my @pindel_indel_paths = ("$build_dir/variants/indel/pindel-0.5-d41d8cd98f00b204e9800998ecf8427e/pindel-somatic-calls-v1-d41d8cd98f00b204e9800998ecf8427e/pindel-vaf-filter-v1-34c9479830c83a54e5d4f73f71e9c660/pindel-read-support-v1-d41d8cd98f00b204e9800998ecf8427e/indels.hq.bed");
    $indel_pindel_results_file = $self->checkResultFile('-paths'=>\@pindel_indel_paths, '-caller'=>"pindel");

    my @varscan_indel_paths = ("$build_dir/variants/indel/varscan-somatic-2.2.6-d41d8cd98f00b204e9800998ecf8427e/varscan-high-confidence-indel-v1-d41d8cd98f00b204e9800998ecf8427e/false-indel-v1-05fbf69c10534fd630b99e44ddf73c7f/indels.hq.bed");
    $indel_varscan_results_file = $self->checkResultFile('-paths'=>\@varscan_indel_paths, '-caller'=>"varscan");

    my ($snv_strelka_results_file, $snv_sniper_results_file, $snv_varscan_results_file);

    #Create a list of possible snv file paths
    my @strelka_snv_paths = ("$build_dir/variants/snv/strelka-0.4.6.2-bc1213eb5850cc1810af3214c95cf30d/snvs.hq.bed",
                             "$build_dir/variants/snv/strelka-0.4.6.2-c5009c08801c3ffa834ecb28d4293d27/snvs.hq.bed",
                             "$build_dir/variants/snv/strelka-0.4.6.2-673995e8237c2c733def86d8d9b3d5a6/snvs.hq.bed",
                             "$build_dir/variants/snv/strelka-0.4.6.2-14acc00d0b01975892118ec71cfc3506/snvs.hq.bed");
    $snv_strelka_results_file = $self->checkResultFile('-paths'=>\@strelka_snv_paths, '-caller'=>"strelka");

    my @sniper_snv_paths = ("$build_dir/variants/snv/sniper-1.0.2-74a151fc61a7a2171177397f4c4f3633/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/somatic-score-mapping-quality-v1-39b60f48b6f8c9e63436a5424305e9fd/snvs.hq.bed");
    $snv_sniper_results_file = $self->checkResultFile('-paths'=>\@sniper_snv_paths, '-caller'=>"sniper");

    my @varscan_snv_paths = ("$build_dir/variants/snv/varscan-somatic-2.2.6-d41d8cd98f00b204e9800998ecf8427e/varscan-high-confidence-v1-d41d8cd98f00b204e9800998ecf8427e/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/snvs.hq.bed");
    $snv_varscan_results_file = $self->checkResultFile('-paths'=>\@varscan_snv_paths, '-caller'=>"varscan");


    #Use 'joinx intersect' to determine which indels in the merged/union file are found in each individual caller's results file
    #gmt joinx intersect a.bed b.bed [--output-file=n.bed] --exact-pos --exact-allele
    my $params_string = "--exact-pos --exact-allele";
    my $indel_strelka_outfile = $build_outdir . "indel_strelka.bed";
    my $indel_gatk_outfile = $build_outdir . "indel_gatk.bed";
    my $indel_pindel_outfile = $build_outdir . "indel_pindel.bed";
    my $indel_varscan_outfile = $build_outdir . "indel_varscan.bed";
    my $snv_strelka_outfile = $build_outdir . "snv_strelka.bed";
    my $snv_sniper_outfile = $build_outdir . "snv_sniper.bed";
    my $snv_varscan_outfile = $build_outdir . "snv_varscan.bed";

    $self->determineCaller("indel","strelka", $indel_results_file, $indel_strelka_results_file, $indel_strelka_outfile);
    $self->determineCaller("indel","gatk", $indel_results_file, $indel_gatk_results_file, $indel_gatk_outfile);
    $self->determineCaller("indel","pindel", $indel_results_file, $indel_pindel_results_file, $indel_pindel_outfile);
    $self->determineCaller("indel","varscan", $indel_results_file, $indel_varscan_results_file, $indel_varscan_outfile);

    $self->determineCaller("snv","strelka", $snv_results_file, $snv_strelka_results_file, $snv_strelka_outfile);
    $self->determineCaller("snv","sniper", $snv_results_file, $snv_sniper_results_file, $snv_sniper_outfile);
    $self->determineCaller("snv","varscan", $snv_results_file, $snv_varscan_results_file, $snv_varscan_outfile);

    #Print out a new file containing the extra source columns
    open (INDEL_OUT, ">$indel_outfile") || die "\n\nCould not open $indel_outfile\n\n";
    print INDEL_OUT "coord\tchr\tstart\tend\tvariant\tscore1\tscore2\tcallers\tstrelka\tgatk\tpindel\tvarscan\ttier\n";

    foreach my $indel (sort {$indels{$a}->{coord_string} cmp $indels{$b}->{coord_string}} keys %indels){
      my @callers = sort keys %{$indel_caller{$indels{$indel}{variant_string}}};
      my $strelka=0; my $gatk=0; my $pindel=0; my $varscan=0;
      foreach my $caller (@callers){
        if ($caller eq 'strelka'){$strelka=1;}
        if ($caller eq 'gatk'){$gatk=1;}
        if ($caller eq 'pindel'){$pindel=1;}
        if ($caller eq 'varscan'){$varscan=1;}
      }
      print INDEL_OUT "$indels{$indel}{coord_string}\t$indels{$indel}{line}\t",join(",",@callers),"\t$strelka\t$gatk\t$pindel\t$varscan\t$indels{$indel}{tier}\n";
    }
    close(INDEL_OUT);

    open (SNV_OUT, ">$snv_outfile") || die "\n\nCould not open $snv_outfile\n\n";
    print SNV_OUT "coord\tchr\tstart\tend\tvariant\tscore1\tscore2\tcallers\tstrelka\tsniper\tvarscan\ttier\n";
    foreach my $snv (sort {$snvs{$a}->{coord_string} cmp $snvs{$b}->{coord_string}} keys %snvs){
      my @callers = sort keys %{$snv_caller{$snvs{$snv}{variant_string}}};
      my $strelka=0; my $sniper=0; my $varscan=0;
      foreach my $caller (@callers){
        if ($caller eq 'strelka'){$strelka=1;}
        if ($caller eq 'sniper'){$sniper=1;}
        if ($caller eq 'varscan'){$varscan=1;}
      }
      print SNV_OUT "$snvs{$snv}{coord_string}\t$snvs{$snv}{line}\t",join(",",@callers),"\t$strelka\t$sniper\t$varscan\t$snvs{$snv}{tier}\n";
    }
    close(INDEL_OUT);

    #Cleanup temp files
    unlink $indel_results_file;
    unlink $snv_results_file;

    #Set output files as output to this step
    die $self->error_message("Trying to set a file as output but the file does not exist: $indel_outfile") unless (-e $indel_outfile);
    $self->indel_variant_sources_file($indel_outfile);
    die $self->error_message("Trying to set a file as output but the file does not exist: $snv_outfile") unless (-e $snv_outfile);
    $self->snv_variant_sources_file($snv_outfile);
  }

  $self->status_message("\n\n");

  return 1;
}

sub noteCaller{
  my $self = shift;
  my $intersect_file=shift;
  my $caller=shift;
  my $variant_type=shift;
  #Go through output of 'joinx intersect' and identify which lines in the main/original file were found intersecting with the other file
  open (INTERSECT, "$intersect_file") || die "\n\ncan't open $intersect_file\n";

  while (<INTERSECT>){
    chomp;
    my @data=split("\t",$_);
    my $variant_string="$data[0]".":"."$data[1]"."-"."$data[2]"." ($data[3])";
    if ($variant_type eq "indel"){
      $indel_caller{$variant_string}{$caller}++;
    }
    if ($variant_type eq "snv"){
      $snv_caller{$variant_string}{$caller}++;
    }
  }
  #when finished, delete intermediate result file from joinx - otherwise this can sometimes interfere with future runs of the tool
  unlink $intersect_file;
}

sub checkResultFile{
  my $self = shift;
  my %args = @_;
  my @paths = @{$args{'-paths'}};
  my $caller = $args{'-caller'};
  my $result_file = '';
  foreach my $path (@paths){
    if (-e $path){
      $result_file=$path;
    }
  }
  unless (-e $result_file){
    my $path_list = join("\n", @paths);
    $self->error_message("$caller result not found in the following list of paths\n\n$path_list\n");
    return undef;
  }
  return($result_file);
}

sub determineCaller {
    my ($self, $variant_type, $caller_name, $results_file, $caller_file, $outfile) = @_;

    if(defined $caller_file) {
        $self->status_message("Looking for overlapping $variant_type results between:\n$results_file\n$caller_file\n\n");

        my $cmd = Genome::Model::Tools::Joinx::Intersect->create(exact_pos=>1, exact_allele=>1, output_file=>$outfile, input_file_a=>$results_file, input_file_b=>$caller_file);
        $cmd->execute();
        
        #Go through original indels and note all files from different callers where that indel was called
        $self->noteCaller($outfile, $caller_name, $variant_type);
    }
}


1;


