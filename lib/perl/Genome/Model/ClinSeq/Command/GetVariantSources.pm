package Genome::Model::ClinSeq::Command::GetVariantSources;

#Written by Obi Griffith

use strict;
use warnings;
use Genome;
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use Genome::Model::ClinSeq::Util qw(:all);

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
	                                          desc => RED . "Outdir: " . $self->outdir . " not found or not a directory" . RESET,
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

    #Locate the final indel/snv results files and load into memory
    #For indels, use ~/effects/indels.hq.novel.tier1.v2.bed ?  (Or the annotated file?)
    #For SNVs, use ~/effects/snvs.hq.novel.tier1.v2.bed
    my $indel_results_file = $build_dir . "/effects/indels.hq.novel.tier1.v2.bed";
    my $snv_results_file = $build_dir . "/effects/snvs.hq.novel.tier1.v2.bed";     

    my %indels;
    open (INDELS, $indel_results_file) or die "can't open $indel_results_file\n";
    my $l=0;
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
    }

    my %snvs;
    open (SNVS, $snv_results_file) or die "can't open $snv_results_file\n";
    $l=0;
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
    }

    #Locate the individual indel/snv files for each caller to use in joinx intersect
    #This should be replaced by a method which somehow determines the appropriate files automatically
    #Depending on whether the somatic variation build is for exome or wgs data, the paths will differ - this should also be determined automatically

    my ($indel_strelka_results_file, $indel_gatk_results_file, $indel_pindel_results_file, $indel_varscan_results_file);
    #my $indel_strelka_results_file_exome = $build_dir . "/variants/indel/strelka-0.4.6.2-bc1213eb5850cc1810af3214c95cf30d/indels.hq.bed";
    my $indel_strelka_results_file_exome = $build_dir . "/variants/indel/strelka-0.4.6.2-c5009c08801c3ffa834ecb28d4293d27/indels.hq.bed";
    my $indel_strelka_results_file_wgs = $build_dir . "/variants/indel/strelka-0.4.6.2-673995e8237c2c733def86d8d9b3d5a6/indels.hq.bed";

    $indel_strelka_results_file=&checkResultFile($indel_strelka_results_file_wgs, $indel_strelka_results_file_exome, "strelka");
    $indel_gatk_results_file = $build_dir . "/variants/indel/gatk-somatic-indel-5336-d41d8cd98f00b204e9800998ecf8427e/false-indel-v1-05fbf69c10534fd630b99e44ddf73c7f/indels.hq.bed";
    $indel_pindel_results_file = $build_dir . "/variants/indel/pindel-0.5-d41d8cd98f00b204e9800998ecf8427e/pindel-somatic-calls-v1-d41d8cd98f00b204e9800998ecf8427e/pindel-vaf-filter-v1-34c9479830c83a54e5d4f73f71e9c660/pindel-read-support-v1-d41d8cd98f00b204e9800998ecf8427e/indels.hq.bed";
    $indel_varscan_results_file = $build_dir . "/variants/indel/varscan-somatic-2.2.6-d41d8cd98f00b204e9800998ecf8427e/varscan-high-confidence-indel-v1-d41d8cd98f00b204e9800998ecf8427e/false-indel-v1-05fbf69c10534fd630b99e44ddf73c7f/indels.hq.bed";


    my ($snv_strelka_results_file, $snv_sniper_results_file, $snv_varscan_results_file);
    #my $snv_strelka_results_file_exome = $build_dir . "/variants/snv/strelka-0.4.6.2-bc1213eb5850cc1810af3214c95cf30d/snvs.hq.bed";
    my $snv_strelka_results_file_exome = $build_dir . "/variants/snv/strelka-0.4.6.2-c5009c08801c3ffa834ecb28d4293d27/snvs.hq.bed";
    my $snv_strelka_results_file_wgs = $build_dir . "/variants/snv/strelka-0.4.6.2-673995e8237c2c733def86d8d9b3d5a6/snvs.hq.bed";
    $snv_strelka_results_file=&checkResultFile($snv_strelka_results_file_wgs, $snv_strelka_results_file_exome, "strelka");
    $snv_sniper_results_file = $build_dir . "/variants/snv/sniper-1.0.2-74a151fc61a7a2171177397f4c4f3633/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/somatic-score-mapping-quality-v1-39b60f48b6f8c9e63436a5424305e9fd/snvs.hq.bed";
    $snv_varscan_results_file = $build_dir . "/variants/snv/varscan-somatic-2.2.6-d41d8cd98f00b204e9800998ecf8427e/varscan-high-confidence-v1-d41d8cd98f00b204e9800998ecf8427e/false-positive-v1-05fbf69c10534fd630b99e44ddf73c7f/snvs.hq.bed";

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

    print "Looking for overlapping indel results between:\n$indel_results_file\n$indel_strelka_results_file\n\n";
    print "Looking for overlapping indel results between:\n$indel_results_file\n$indel_gatk_results_file\n\n";
    print "Looking for overlapping indel results between:\n$indel_results_file\n$indel_pindel_results_file\n\n";
    print "Looking for overlapping indel results between:\n$indel_results_file\n$indel_varscan_results_file\n\n";
    print "Looking for overlapping snv results between:\n$snv_results_file\n$snv_strelka_results_file\n\n";
    print "Looking for overlapping snv results between:\n$snv_results_file\n$snv_sniper_results_file\n\n";
    print "Looking for overlapping snv results between:\n$snv_results_file\n$snv_varscan_results_file\n\n";

    my $joinx_indel_strelka_cmd = "gmt joinx intersect $indel_results_file $indel_strelka_results_file $params_string --output-file $indel_strelka_outfile";
    my $joinx_indel_gatk_cmd = "gmt joinx intersect $indel_results_file $indel_gatk_results_file $params_string --output-file $indel_gatk_outfile";
    my $joinx_indel_pindel_cmd = "gmt joinx intersect $indel_results_file $indel_pindel_results_file $params_string --output-file $indel_pindel_outfile";
    my $joinx_indel_varscan_cmd = "gmt joinx intersect $indel_results_file $indel_varscan_results_file $params_string --output-file $indel_varscan_outfile";
    my $joinx_snv_strelka_cmd = "gmt joinx intersect $snv_results_file $snv_strelka_results_file $params_string --output-file $snv_strelka_outfile";
    my $joinx_snv_sniper_cmd = "gmt joinx intersect $snv_results_file $snv_sniper_results_file $params_string --output-file $snv_sniper_outfile";
    my $joinx_snv_varscan_cmd = "gmt joinx intersect $snv_results_file $snv_varscan_results_file $params_string --output-file $snv_varscan_outfile";

    Genome::Sys->shellcmd(cmd => $joinx_indel_strelka_cmd);
    Genome::Sys->shellcmd(cmd => $joinx_indel_gatk_cmd);
    Genome::Sys->shellcmd(cmd => $joinx_indel_pindel_cmd);
    Genome::Sys->shellcmd(cmd => $joinx_indel_varscan_cmd);
    Genome::Sys->shellcmd(cmd => $joinx_snv_strelka_cmd);
    Genome::Sys->shellcmd(cmd => $joinx_snv_sniper_cmd);
    Genome::Sys->shellcmd(cmd => $joinx_snv_varscan_cmd);

    #Go through original indels and note all files from different callers where that indel was called
    &noteCaller($indel_strelka_outfile, "strelka", "indel");
    &noteCaller($indel_gatk_outfile, "gatk", "indel");
    &noteCaller($indel_pindel_outfile, "pindel", "indel");
    &noteCaller($indel_varscan_outfile, "varscan", "indel");
    &noteCaller($snv_strelka_outfile, "strelka", "snv");
    &noteCaller($snv_sniper_outfile, "sniper", "snv");
    &noteCaller($snv_varscan_outfile, "varscan", "snv");

    #Print out a new file containing the extra source columns
    open (INDEL_OUT, ">$indel_outfile") || die "\n\nCould not open $indel_outfile\n\n";
    print INDEL_OUT "chr\tstart\tend\tvariant\tscore1\tscore2\tcallers\tstrelka\tgatk\tpindel\tvarscan\tcoord_string\n";

    foreach my $indel (sort keys %indels){
      my @callers = sort keys %{$indel_caller{$indels{$indel}{variant_string}}};
      my $strelka=0; my $gatk=0; my $pindel=0; my $varscan=0;
      foreach my $caller (@callers){
        if ($caller eq 'strelka'){$strelka=1;}
        if ($caller eq 'gatk'){$gatk=1;}
        if ($caller eq 'pindel'){$pindel=1;}
        if ($caller eq 'varscan'){$varscan=1;}
      }
      print INDEL_OUT "$indels{$indel}{line}\t",join(",",@callers),"\t$strelka\t$gatk\t$pindel\t$varscan","\t$indels{$indel}{coord_string}","\n";
    }
    close(INDEL_OUT);

    open (SNV_OUT, ">$snv_outfile") || die "\n\nCould not open $snv_outfile\n\n";
    print SNV_OUT "chr\tstart\tend\tvariant\tscore1\tscore2\tcallers\tstrelka\tsniper\tvarscan\tcoord_string\n";
    foreach my $snv (sort keys %snvs){
      my @callers = sort keys %{$snv_caller{$snvs{$snv}{variant_string}}};
      my $strelka=0; my $sniper=0; my $varscan=0;
      foreach my $caller (@callers){
        if ($caller eq 'strelka'){$strelka=1;}
        if ($caller eq 'sniper'){$sniper=1;}
        if ($caller eq 'varscan'){$varscan=1;}
      }
      print SNV_OUT "$snvs{$snv}{line}\t",join(",",@callers),"\t$strelka\t$sniper\t$varscan","\t$snvs{$snv}{coord_string}","\n";
    }
    close(INDEL_OUT);
  }
  $self->status_message("\n\n");

  return 1;
}

sub noteCaller{
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
  my $rm_command = "rm $intersect_file";
  Genome::Sys->shellcmd(cmd => $rm_command);
}

sub checkResultFile{
  my $wgs_file = shift;
  my $exome_file = shift;
  my $caller = shift;
  my $result_file;
  if (-e $wgs_file){
    $result_file=$wgs_file;
  }elsif (-e $exome_file){
    $result_file=$exome_file;
  }else{
    print "$caller result not found\n";
    exit;
  }
return($result_file);
}

1;


