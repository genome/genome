package Genome::Model::ClinSeq::Command::CreateMutationSpectrum;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::CreateMutationSpectrum {
    is => ['Command::V2',
      'Genome::Model::ClinSeq::Util'],
    has_input => [
        somvar_build => { 
            is => 'Genome::Model::Build::SomaticVariation',
            shell_args_position => 1,
            require_user_verify => 0,
            doc => 'somatic variation build(s) to '.
                'create mutation spectrum results from',
        },
        clinseq_build => { 
            is => 'Genome::Model::Build::ClinSeq',
            doc => 'clinseq build to '.
                'create mutation spectrum results from',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
        datatype => {
              is => 'Text',
              doc => 'wgs or exome datatype',
        },
        limit_snvs => {
              is => 'Number',
              doc => 'Limit analysis to the first N SNVs.' .
                'This parameter is used only when "test" is true',
              is_optional => 1,
              default => 50,
        },
        min_quality_score => {
            is => 'Integer',
            doc => 'minimum mapping quality of reads to be considered',
            is_optional => 1,
            default => '20',
        },
        min_base_quality => {
            is => 'Integer',
            doc => 'minimum base quality of bases in reads to be considered',
            is_optional => 1,
            default => '30',
        },
        converge_snv_indel_report_result => {
            is => 'Boolean',
            doc => 'Link in workflow',
            is_optional => 1
        },
        test => {
            is => 'Boolean',
            doc => 'set for test-cases',
            is_optional => 1,
            default => 0,
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
Should work for either wgs or exome data but different sets of variants will be used.
The read-counts are obtained from the SNVINdel report directory.
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

sub get_final_name {
  my $self = shift;
  my $somvar_build = shift;
  my $somvar_build_id = $somvar_build->model->id;
  my $final_name = $somvar_build->model->subject->name if ($somvar_build->model->subject->name);
  $final_name = $somvar_build->model->subject->individual->common_name if ($somvar_build->model->subject->individual->common_name);
  return $final_name;
}

sub make_outdirs {
  my $self = shift;
  my $datatype = shift;
  my $outdir = $self->outdir;
  $outdir .= "/" unless ($outdir =~ /\/$/);
  Genome::Sys->create_directory($outdir);
  my $sub_outdir;
  if ($datatype =~ /wgs/i){
    $sub_outdir = $outdir . "wgs/";
  }
  if ($datatype =~ /exome/i){
    $sub_outdir = $outdir . "exome/";
  }
  my $sub_outdir2 = $sub_outdir . "mutation_spectrum_sequence_context/";
  my $sub_outdir3 = $sub_outdir . "summarize_mutation_spectrum/";
  my $sub_outdir4 = $sub_outdir . "mutation_rate/";
  Genome::Sys->create_directory($sub_outdir);
  Genome::Sys->create_directory($sub_outdir2);
  Genome::Sys->create_directory($sub_outdir3);
  Genome::Sys->create_directory($sub_outdir4);
  return ($sub_outdir, $sub_outdir2, $sub_outdir3, $sub_outdir4);
}

sub parse_variant_file {
  my $self = shift;
  my $clinseq_build = shift;
  my $variant_file = shift;
  my $data_type = $self->datatype;
  my $gender = $clinseq_build->subject->gender;
  my @headers = qw/chr start end ref_allele var_allele/;
  my @tumor_prefixes = $self->_get_si_report_tumor_prefix($clinseq_build);
  my @tiers = ("tier1", "tier2", "tier3", "tier1-3");
  my $variant_file_temp = $variant_file . "_" . $data_type;
  $variant_file_temp =~ s/\.tsv//;
  foreach my $tumor_prefix(@tumor_prefixes) {
    unless ($tumor_prefix =~ $data_type) {
      next;
    }
    my $writer_t1 = Genome::Utility::IO::SeparatedValueWriter->create(
      output => $variant_file_temp . "." . $tiers[0] . ".tsv",
      separator => "\t",
      headers => \@headers,
      print_headers => 0,
    );
    my $writer_t2 = Genome::Utility::IO::SeparatedValueWriter->create(
      output => $variant_file_temp . "." . $tiers[1] . ".tsv",
      separator => "\t",
      headers => \@headers,
      print_headers => 0,
    );
    my $writer_t3 = Genome::Utility::IO::SeparatedValueWriter->create(
      output => $variant_file_temp . "." . $tiers[2] . ".tsv",
      separator => "\t",
      headers => \@headers,
      print_headers => 0,
    );
    my $writer_t4 = Genome::Utility::IO::SeparatedValueWriter->create(
      output => $variant_file_temp . "." . $tiers[3] . ".tsv",
      separator => "\t",
      headers => \@headers,
      print_headers => 0,
    );
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
      separator => "\t",
      input => $variant_file,
    );
    my $out_data;
    while (my $data = $reader->next) {
      unless ($data->{type} =~ /SNP/) {
        next;
      }
      if ($data->{chromosome_name} =~ /Y|MT|GL/) {
        next;
      }
      unless ($data->{data_type} =~ /$data_type/) {
        next;
      }
      if ($gender ne "female" and $data->{chromosome_name} =~ /X/) {
        next;
      }
      $out_data->{chr} = $data->{chromosome_name};
      $out_data->{start} = $data->{start};
      $out_data->{end} = $data->{start};
      $out_data->{ref_allele} = $data->{reference};
      $out_data->{var_allele} = $data->{variant};
      if($data->{tier} eq "tier1") {
        $writer_t1->write_one($out_data);
      }
      if($data->{tier} eq "tier2") {
        $writer_t2->write_one($out_data);
      }
      if($data->{tier} eq "tier3") {
        $writer_t3->write_one($out_data);
      }
      $writer_t4->write_one($out_data);
    }
  }
  my @variant_files;
  foreach my $tier(@tiers) {
    push @variant_files, $variant_file_temp . "." . $tier . ".tsv";  
  }
  return @variant_files;
}

sub get_alltier_snvs {
  my $self = shift;
  my $clinseq_build = shift;
  my $sub_outdir = shift;
  my $outfile = $sub_outdir . "/alltiers.variants.filtered.clean.tsv";
  my $bq = $self->min_base_quality;
  my $mq = $self->min_quality_score;
  my $snv_indel_report_clean_file =
  $clinseq_build->snv_indel_report_clean_filtered_file(
    $bq, $mq);
  if(-e $snv_indel_report_clean_file) {
    Genome::Sys->copy_file($snv_indel_report_clean_file, $outfile);
  } else {
    die $self->error_message("Unable to find variant read-counts file
      for clinseq build " . $clinseq_build->id);
  }
  return $self->parse_variant_file($clinseq_build, $outfile);
}

sub reduce_file_length {
  my $self = shift;
  my @variant_files = @_;
  my $max_snvs = $self->limit_snvs;
  unless ($max_snvs > 0){
    die $self->error_message("Invalid max_snvs $max_snvs");
  }
  foreach my $variant_file(@variant_files) {
    my $tmp_file = $variant_file . ".tmp";
    my $cmd = "head -n " . $max_snvs . " $variant_file > $tmp_file";
    $self->debug_message($cmd);
    Genome::Sys->shellcmd(cmd => $cmd);
    $cmd = "mv $tmp_file $variant_file";
    $self->debug_message($cmd);
    Genome::Sys->shellcmd(cmd => $cmd);
  }
}

sub get_mutation_spectrum_sequence_context_result {
  my $self = shift;
  my $final_name = shift;
  my $tier1_snvs = shift;
  my $sub_outdir2 = shift;
  my $somvar_build = $self->somvar_build;
  my $reference_sequence_build = $somvar_build->tumor_model->reference_sequence_build;
  my $reference_fasta_path = $reference_sequence_build->full_consensus_path('fa');
  my $mssc_file4plot = $sub_outdir2 . $final_name . ".data.tsv";
  my $mssc_outfile = $sub_outdir2 . $final_name .".mutation-spectrum-sequence-context.pdf";
  my $mssc_proportiontest_outfile =  $sub_outdir2 . $final_name . ".prop.test";
  my $plot_title = "Mutation Spectrum Sequence Context for " . $final_name;
  my $mssc_cmd = Genome::Model::Tools::Analysis::MutationSpectrumSequenceContext->create(output_file=>$mssc_outfile, roi_file=>$tier1_snvs, file4plot=>$mssc_file4plot, plot_title=>$plot_title, ref_seq=>$reference_fasta_path, proportiontest=>$mssc_proportiontest_outfile, random_trials=>100, random_seed=>2013, window_size=>10);
  $mssc_cmd->execute();
}

sub generate_summarize_mutation_spectrum_result {
  my $self = shift;
  my $somvar_build = shift;
  my $final_name = shift;
  my $tier1to3snvs = shift;
  my $sub_outdir3 = shift;
  my $somatic_id = $somvar_build->model->id;
  my $mut_spec_file = $sub_outdir3 . "mutation_spectrum.tsv";
  my $sms_file = $sub_outdir3 . $final_name . "_summarize-mutation-spectrum.pdf";
  my $input_file = $sub_outdir3 . "/mutation_spectrum.input.tsv";
  open(my $SMS, ">$input_file") or
    die $self->error_message("Unable to open $input_file for writing");
  print $SMS $final_name, "\t", $tier1to3snvs;
  close $SMS;
  my $sms_cmd = Genome::Model::Tools::Analysis::SummarizeMutationSpectrum->create(exclude_gl_contigs=>1, somatic_id=>$somatic_id, mut_spec_file=>$mut_spec_file, output_file=>$sms_file, input_file=>$input_file);
  $sms_cmd->execute();
}

sub generate_mutation_rate_result {
  my $self = shift;
  my $somvar_build = shift;
  my $final_name = shift;
  my $tier1_snvs = shift;
  my $tier2_snvs = shift;
  my $tier3_snvs = shift;
  my $sub_outdir4 = shift;
  my $data_type = $self->datatype;
  my $reference_annotation_dir = $somvar_build->model->annotation_build->data_directory;
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
  if ($data_type =~ /wgs/i){
    #WGS Tier1 only
    $mutation_rate_outfile = $sub_outdir4 . "mutation-rate-tier1.tsv";
    $mutation_rate_cmd = Genome::Model::Tools::Analysis::MutationRate->create(sample_name=>$final_name, tier1_file=>$tier1_snvs, tier1_space=>$tier1_bed, outfile=>$mutation_rate_outfile);
    $mutation_rate_cmd->execute();
    #WGS Tier1-3
    $mutation_rate_outfile = $sub_outdir4 . "mutation-rate-tier123.tsv";
    $mutation_rate_cmd = Genome::Model::Tools::Analysis::MutationRate->create(sample_name=>$final_name, tier1_file=>$tier1_snvs, tier1_space=>$tier1_bed, tier2_file=>$tier2_snvs, tier2_space=>$tier2_bed, tier3_file=>$tier3_snvs, tier3_space=>$tier3_bed, outfile=>$mutation_rate_outfile);
    $mutation_rate_cmd->execute();
  }
  if ($data_type =~ /exome/i){
    #Exome Tier1 only
    $mutation_rate_outfile = $sub_outdir4 . "mutation-rate-tier1.tsv";
    $mutation_rate_cmd = Genome::Model::Tools::Analysis::MutationRate->create(sample_name=>$final_name, tier1_file=>$tier1_snvs, tier1_space=>$tier1_bed, tier2_space=>$tier2_bed, tier3_space=>$tier3_bed, coverage_factor=>0.02, outfile=>$mutation_rate_outfile);
    $mutation_rate_cmd->execute();
  }
}

sub create_readme_file {
  my $self = shift;
  my $sub_outdir4 = shift;
  my $readme_file = $sub_outdir4 . "readme.txt";
  open (RM, ">$readme_file") || die "\n\nCould not open $readme_file for writing\n\n";
  print RM "OUTPUT format is TSV (SampleName then Mutation counts then Mutation rates per megabase).\n";
  print RM "Depends on which tiers are specified but in general of the form:\n";
  print RM "SampleName Tier1_Count ... Tier3_Count NonTier1_Count Overall_Count Tier1_Rate ... Tier3_Rate NonTier1_Rate Overall_Rate\n";
  close(RM);
}

sub execute {
  my $self = shift;
  my $somvar_build = $self->somvar_build;
  my $clinseq_build = $self->clinseq_build;
  my $data_type = $self->datatype;
  $self->debug_message("Performing mutation spectrum analysis");
  my ($sub_outdir, $sub_outdir2, $sub_outdir3, $sub_outdir4) =
    $self->make_outdirs($data_type);
  my $final_name = $self->get_final_name($somvar_build);
  my ($tier1_snvs, $tier2_snvs, $tier3_snvs, $tier1to3_snvs) = $self->get_alltier_snvs($clinseq_build, $sub_outdir);
  if($self->test) {
    $self->reduce_file_length($tier1_snvs, $tier2_snvs, $tier3_snvs, $tier1to3_snvs)
  }
  $self->get_mutation_spectrum_sequence_context_result($final_name, $tier1_snvs, $sub_outdir2);
  $self->generate_summarize_mutation_spectrum_result($somvar_build, $final_name, $tier1to3_snvs, $sub_outdir3);
  $self->generate_mutation_rate_result($somvar_build, $final_name, $tier1_snvs, $tier2_snvs, $tier3_snvs, $sub_outdir4);  
  $self->create_readme_file($sub_outdir4);
  return 1;
}


1;

