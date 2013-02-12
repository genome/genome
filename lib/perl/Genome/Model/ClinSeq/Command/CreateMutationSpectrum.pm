package Genome::Model::ClinSeq::Command::CreateMutationSpectrum;

use strict;
use warnings;
use Genome;
use Data::Dumper;

class Genome::Model::ClinSeq::Command::CreateMutationSpectrum {
    is => 'Command::V2',
    has_input => [
        build => { 
            is => 'Genome::Model::Build::SomaticVariation',
            shell_args_position => 1,
            require_user_verify => 0,
            doc => 'somatic variation build(s) to create mutation spectrum results from',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
        datatype => {
              is => 'Text',
              doc => 'wgs or exome datatype',
        },
        max_snvs => {
              is => 'Number',
              doc => 'Limit analysis to the first N SNVs',
              is_optional => 1,
        },
    ],
    doc => 'analyze the mutation spectrum of wgs or exome variants',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq create-mutation-spectrum --outdir=/tmp/create_mutation_spectrum/ --datatype=wgs 129855269

EOS
}

sub help_detail {
    return <<EOS
Create mutation spectrum results for a single somatic variation build.

Should work for either wgs or exome data but different sets of variants will be used

EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless ($self->datatype =~ /wgs|exome/i) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['datatype'],
	                                          desc => "datatype: " . $self->datatype . " must be wgs or exome",
                                          );
  }
  return @errors;
}


sub execute {
  my $self = shift;
  $self->status_message("Performing mutation spectrum analysis");

  my $build = $self->build;
  my $outdir = $self->outdir;
  my $datatype = $self->datatype;

  #Make directories for results
  my $sub_outdir;
  $outdir .= "/" unless ($outdir =~ /\/$/);
  unless (-e $outdir && -d $outdir){
    mkdir($outdir);
  }
  if ($datatype =~ /wgs/i){
    $sub_outdir = $outdir . "wgs/";
  }
  if ($datatype =~ /exome/i){
    $sub_outdir = $outdir . "exome/";
  }
  my $sub_outdir2 = $sub_outdir . "mutation_spectrum_sequence_context/";
  my $sub_outdir3 = $sub_outdir . "summarize_mutation_spectrum/";
  my $sub_outdir4 = $sub_outdir . "mutation_rate/";
  mkdir $sub_outdir unless (-e $sub_outdir && -d $sub_outdir);
  mkdir $sub_outdir2 unless (-e $sub_outdir2 && -d $sub_outdir2);
  mkdir $sub_outdir3 unless (-e $sub_outdir3 && -d $sub_outdir3);
  mkdir $sub_outdir4 unless (-e $sub_outdir4 && -d $sub_outdir4);

  #Run various mutation-spectrum tools in various modes.  
  #TODO: accommodate new features pushed by Charles Lu

  #Get data directory, reference annotation etc.
  my $data_directory = $build->data_directory;
  my $reference_annotation_name = $build->model->annotation_build->name;
  my $reference_annotation_dir = $build->model->annotation_build->data_directory;
  my $reference_sequence_build = $build->tumor_model->reference_sequence_build;
  my $reference_fasta_path = $reference_sequence_build->full_consensus_path('fa');
  my $reference_sequence_dir = $reference_sequence_build->data_directory;
  
  #Find tier1-3 SNV bed files
  my $tier1_snvs = "$data_directory/effects/snvs.hq.novel.tier1.v2.bed";
  my $tier2_snvs = "$data_directory/effects/snvs.hq.novel.tier2.v2.bed";
  my $tier3_snvs = "$data_directory/effects/snvs.hq.novel.tier3.v2.bed";

  unless (-e $tier1_snvs && -e $tier2_snvs && $tier3_snvs){
    $self->error_message("Could not find snv files at expected path: $data_directory" . "/effects/");
    exit 1;
  }

  #Get a 'final' name for the sample
  my $final_name = $build->model->id;
  $final_name = $build->model->subject->name if ($build->model->subject->name);
  $final_name = $build->model->subject->patient->common_name if ($build->model->subject->patient->common_name);

  #1.) Get variants from somatic variation build (tier1-3 for wgs and tier1 for exome?)
  my $variant_file;
  if ($datatype =~ /wgs/i){
    $variant_file = $sub_outdir2 . "snvs.hq.novel.tier123.v2.bed";
    my $snv_cat_cmd = "cat $tier1_snvs $tier2_snvs $tier3_snvs > $variant_file";
    $self->status_message($snv_cat_cmd);
    Genome::Sys->shellcmd(cmd => $snv_cat_cmd, output_files=>["$variant_file"]);  
  }
  if ($datatype =~ /exome/i){
    $variant_file = $sub_outdir2 . "snvs.hq.novel.tier1.v2.bed";
    my $snv_cat_cmd = "cat $tier1_snvs > $variant_file";
    $self->status_message($snv_cat_cmd);
    Genome::Sys->shellcmd(cmd => $snv_cat_cmd, output_files=>["$variant_file"]);
  }

  #Reduce file length if max_snvs was supplied
  if ($self->max_snvs){
    my $tmp_file = $variant_file . ".tmp";
    my $cmd = "grep -v GL $variant_file > $tmp_file; mv $tmp_file $variant_file; head -n " . $self->max_snvs . " $variant_file > $tmp_file; mv $tmp_file $variant_file";
    $self->status_message($cmd);
    Genome::Sys->shellcmd(cmd => $cmd);
  }

  #2.) Run annotator on all Tier1-3 variants
  #gmt annotate transcript-variants --variant-bed-file='' --output-file='' --no-headers? --annotation-filter=top  --use-version=3 --reference-transcripts=''
  my $annotated_file = $variant_file . ".annotated";
  my $annotate_stdout = $sub_outdir2 . "gmt-annotate.stdout";
  my $annotate_stderr = $sub_outdir2 . "gmt-annotate.stderr";   
  my $annotate_cmd = "gmt annotate transcript-variants --variant-bed-file=$variant_file --output-file=$annotated_file --no-headers --annotation-filter=top  --use-version=3 --reference-transcripts=$reference_annotation_name --skip-if-output-present 1>$annotate_stdout 2>$annotate_stderr";
  $self->status_message($annotate_cmd);
  Genome::Sys->shellcmd(cmd => $annotate_cmd, output_files=>["$annotated_file"]);
  
  #5.) Generate mutation-spectrum-sequence-context result
  #gmt analysis mutation-spectrum-sequence-context
  my $variant_file2 = $sub_outdir2 . "variants.tsv";
  my $cut_cmd = "cut -f 1-5 $annotated_file > $variant_file2";
  $self->status_message($cut_cmd);
  Genome::Sys->shellcmd(cmd => $cut_cmd);

  my $mssc_file4plot = $sub_outdir2 . "data.tsv";
  my $mssc_outfile = $sub_outdir2 . $final_name ."_mutation-spectrum-sequence-context.pdf";
  my $mssc_stdout = $sub_outdir2 . "gmt-mutation-spectrum-sequence-context.stdout";
  my $mssc_stderr = $sub_outdir2 . "gmt-mutation-spectrum-sequence-context.stderr";
  my $mssc_cmd = "gmt analysis mutation-spectrum-sequence-context --output-file=$mssc_outfile --roi-file=$variant_file2 --file4plot=$mssc_file4plot --plot-title='Mutation Spectrum Sequence Context for $final_name' --ref-seq=$reference_fasta_path --window-size=10 1>$mssc_stdout 2>$mssc_stderr";
  $self->status_message($mssc_cmd);
  Genome::Sys->shellcmd(cmd => $mssc_cmd);

  #6.) Generate summarize-mutation-spectrum result
  #gmt analysis summarize-mutation-spectrum
  my $somatic_id = $build->model->id;
  my $mut_spec_file = $sub_outdir3 . "mutation_spectrum.tsv";
  my $sms_file = $sub_outdir3 . $final_name . "_summarize-mutation-spectrum.pdf";
  my $sms_stdout = $sub_outdir3 . "summarize-mutation-spectrum.stdout";
  my $sms_stderr = $sub_outdir3 . "summarize-mutation-spectrum.stderr";
  my $sms_cmd = "gmt analysis summarize-mutation-spectrum --exclude-gl-contigs --somatic-id=$somatic_id --mut-spec-file=$mut_spec_file --output-file=$sms_file 1>$sms_stdout 2>$sms_stderr";
  $self->status_message($sms_cmd);
  Genome::Sys->shellcmd(cmd => $sms_cmd);

  #7.) Generate mutation-rate result

  #Get tier bed files from annotation dir
  my $tier1_bed = $reference_annotation_dir . "/annotation_data/tiering_bed_files_v3/tier1.bed";
  my $tier2_bed = $reference_annotation_dir . "/annotation_data/tiering_bed_files_v3/tier2.bed";
  my $tier3_bed = $reference_annotation_dir . "/annotation_data/tiering_bed_files_v3/tier3.bed";

  unless (-e $tier1_bed && -e $tier2_bed && -e $tier3_bed){
    $self->error_message("Could not find tiering files at expected path: $reference_annotation_dir" . " /annotation_data/tiering_bed_files_v3/");
    exit 1;
  }

  my $mutation_rate_cmd;
  my $mutation_rate_outfile;
  my $mr_stderr = $sub_outdir4 . "mutation-rate.stderr";

  if ($datatype =~ /wgs/i){
    #WGS Tier1 only
    $mutation_rate_outfile = $sub_outdir4 . "mutation-rate-tier1.tsv";
    $mutation_rate_cmd = "gmt analysis mutation-rate --sample-name=$final_name --tier1-file=$tier1_snvs --tier1-space=$tier1_bed > $mutation_rate_outfile 2>>$mr_stderr";
    $self->status_message($mutation_rate_cmd);
    Genome::Sys->shellcmd(cmd => $mutation_rate_cmd);

    #WGS Tier1-3
    $mutation_rate_outfile = $sub_outdir4 . "mutation-rate-tier123.tsv";
    $mutation_rate_cmd = "gmt analysis mutation-rate --sample-name=$final_name --tier1-file=$tier1_snvs --tier1-space=$tier1_bed --tier2-file=$tier2_snvs --tier2-space=$tier2_bed --tier3-file=$tier3_snvs --tier3-space=$tier3_bed > $mutation_rate_outfile 2>>$mr_stderr";
    $self->status_message($mutation_rate_cmd);
    Genome::Sys->shellcmd(cmd => $mutation_rate_cmd);
  }

  if ($datatype =~ /exome/i){
    #Exome Tier1 only
    $mutation_rate_outfile = $sub_outdir4 . "mutation-rate-tier1.tsv";
    $mutation_rate_cmd = "gmt analysis mutation-rate --sample-name=$final_name --tier1-file=$tier1_snvs --tier1-space=$tier1_bed --tier2-space=$tier2_bed --tier3-space=$tier3_bed --coverage-factor=0.02 > $mutation_rate_outfile 2>>$mr_stderr";
    $self->status_message($mutation_rate_cmd);
    Genome::Sys->shellcmd(cmd => $mutation_rate_cmd);
  
  }

  my $readme_file = $sub_outdir4 . "readme.txt";
  open (RM, ">$readme_file") || die "\n\nCould not open $readme_file for writing\n\n";
  print RM "OUTPUT format is TSV (SampleName then Mutation counts then Mutation rates per megabase).\n";
  print RM "Depends on which tiers are specified but in general of the form:\n";
  print RM "SampleName Tier1_Count ... Tier3_Count NonTier1_Count Overall_Count Tier1_Rate ... Tier3_Rate NonTier1_Rate Overall_Rate\n";
  close(RM);

  return 1;
}


1;

