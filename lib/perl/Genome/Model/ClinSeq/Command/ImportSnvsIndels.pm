package Genome::Model::ClinSeq::Command::ImportSnvsIndels;

use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq::Util qw(:all);

class Genome::Model::ClinSeq::Command::ImportSnvsIndels {
    is => 'Command::V2',
    has_input => [
        wgs_build => { 
              is => 'Genome::Model::Build::SomaticVariation',
              is_many => 0,
              is_optional => 1,
              doc => 'wgs somatic variation build(s) to get SNVs indels',
        },
        exome_build => { 
              is => 'Genome::Model::Build::SomaticVariation',
              is_many => 0,
              is_optional => 1,
              doc => 'exome somatic variation build(s) to get SNVs indels',
        },
        cancer_annotation_db => {
              is => 'Genome::Db::Tgi::CancerAnnotation',
              doc => 'cancer annotation db',
        },
        filter_mt => {
            is => 'Boolean',
            doc => 'Do not import MT variants',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
    ],
    has_output => [
        wgs_snv_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
        exome_snv_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
        wgs_exome_snv_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
        wgs_indel_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
        exome_indel_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
        wgs_exome_indel_file => {
            is => 'FilesystemPath',
            is_optional => 1,
        },
    ],
    doc => 'gather and reformat SNVs and Indels from a wgs somatic build, an exome somatic build or both',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq import-snvs-indels --outdir=/tmp/ --wgs-build=129396794 --exome-build=129396799  --filter-mt

EOS
}

sub help_detail {
    return <<EOS
summarize snvs and indels for a  wgs somatic build, an exome somatic build or both

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
  unless ($self->wgs_build || $self->exome_build) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['wgs_build','exome_build'],
	                                          desc => "Must define either a wgs or exome build (or both)",
                                          );
  }

  return @errors;
}

sub execute {
  my $self = shift;
  my $outdir = $self->outdir;
  my $wgs_build = $self->wgs_build;
  my $exome_build = $self->exome_build;
  my $filter_mt = $self->filter_mt;
  my $cancer_annotation_db = $self->cancer_annotation_db;

  unless ($outdir =~ /\/$/){
    $outdir .= "/";
  }

  #Get Entrez and Ensembl data for gene name mappings
  my $entrez_ensembl_data = &loadEntrezEnsemblData(-cancer_db => $cancer_annotation_db);

  my $snv_dir = $outdir . "snv/";
  Genome::Sys->create_directory($snv_dir);

  my $indel_dir = $outdir . "indel/";
  Genome::Sys->create_directory($indel_dir);

  #Define variant effect type filters
  #TODO: Allow different filters to be used as a parameter
  my $snv_filter = "missense|nonsense|splice_site|splice_region|rna";
  my $indel_filter = "in_frame_del|in_frame_ins|frame_shift_del|frame_shift_ins|splice_site_ins|splice_site_del|rna";

  #Define the dataset: WGS SNV, WGS indel, Exome SNV, Exome indel
  my %dataset;
  if ($wgs_build){
    my $snv_wgs_dir = $snv_dir . "wgs/";
    Genome::Sys->create_directory($snv_wgs_dir);
    my $indel_wgs_dir = $indel_dir . "wgs/";
    Genome::Sys->create_directory($indel_wgs_dir);
    my $effects_dir = $wgs_build->data_directory . "/effects/";
    $dataset{'1'}{data_type} = "wgs";
    $dataset{'1'}{var_type} = "snv";
    $dataset{'1'}{effects_dir} = $effects_dir;
    $dataset{'1'}{t1_hq_annotated} = "snvs.hq.tier1.v1.annotated";
    $dataset{'1'}{t1_hq_annotated_top} = "snvs.hq.tier1.v1.annotated.top";
    $dataset{'1'}{compact_file} = "$snv_wgs_dir"."snvs.hq.tier1.v1.annotated.compact.tsv";
    $dataset{'1'}{aa_effect_filter} = $snv_filter;
    $dataset{'1'}{target_dir} = $snv_wgs_dir;

    $dataset{'2'}{data_type} = "wgs";
    $dataset{'2'}{var_type} = "indel";
    $dataset{'2'}{effects_dir} = $effects_dir;
    $dataset{'2'}{t1_hq_annotated} = "indels.hq.tier1.v1.annotated";
    $dataset{'2'}{t1_hq_annotated_top} = "indels.hq.tier1.v1.annotated.top";
    $dataset{'2'}{compact_file} = "$indel_wgs_dir"."indels.hq.tier1.v1.annotated.compact.tsv";
    $dataset{'2'}{aa_effect_filter} = $indel_filter;
    $dataset{'2'}{target_dir} = $indel_wgs_dir;
  }
  if ($exome_build){
    my $snv_exome_dir = $snv_dir . "exome/";
    Genome::Sys->create_directory($snv_exome_dir);
    my $indel_exome_dir = $indel_dir . "exome/";
    Genome::Sys->create_directory($indel_exome_dir);
    my $effects_dir = $exome_build->data_directory . "/effects/";
    $dataset{'3'}{data_type} = "exome";
    $dataset{'3'}{var_type} = "snv";
    $dataset{'3'}{effects_dir} = $effects_dir;
    $dataset{'3'}{t1_hq_annotated} = "snvs.hq.tier1.v1.annotated";
    $dataset{'3'}{t1_hq_annotated_top} = "snvs.hq.tier1.v1.annotated.top";
    $dataset{'3'}{compact_file} = "$snv_exome_dir"."snvs.hq.tier1.v1.annotated.compact.tsv";
    $dataset{'3'}{aa_effect_filter} = $snv_filter;
    $dataset{'3'}{target_dir} = $snv_exome_dir;

    $dataset{'4'}{data_type} = "exome";
    $dataset{'4'}{var_type} = "indel";
    $dataset{'4'}{effects_dir} = $effects_dir;
    $dataset{'4'}{t1_hq_annotated} = "indels.hq.tier1.v1.annotated";
    $dataset{'4'}{t1_hq_annotated_top} = "indels.hq.tier1.v1.annotated.top";
    $dataset{'4'}{compact_file} = "$indel_exome_dir"."indels.hq.tier1.v1.annotated.compact.tsv";
    $dataset{'4'}{aa_effect_filter} = $indel_filter;
    $dataset{'4'}{target_dir} = $indel_exome_dir;
  }

  my %data_merge;
  foreach my $ds (sort {$a <=> $b} keys %dataset){
    my %data_out;

    #Make a copy of the high quality .annotated and .annotated.top files
    my $data_type = $dataset{$ds}{data_type};
    my $var_type = $dataset{$ds}{var_type};
    my $effects_dir = $dataset{$ds}{effects_dir};
    my $t1_hq_annotated = $dataset{$ds}{t1_hq_annotated};
    my $t1_hq_annotated_top = $dataset{$ds}{t1_hq_annotated_top};
    my $compact_file = $dataset{$ds}{compact_file};
    my $aa_effect_filter = $dataset{$ds}{aa_effect_filter};
    my $target_dir = $dataset{$ds}{target_dir};

    my $new_annotated_file = "$target_dir$t1_hq_annotated".".tsv";
    my $new_annotated_top_file = "$target_dir$t1_hq_annotated_top".".tsv";
    Genome::Sys->copy_file("$effects_dir$t1_hq_annotated",$new_annotated_file);
    Genome::Sys->copy_file("$effects_dir$t1_hq_annotated_top", $new_annotated_top_file);

    #Get a column count on the file and use this to determine the correct header for the annotated variant file
    my @input_headers; 
    #If empty file, just set header otherwise check number of columns before setting header
    if (-z $new_annotated_file){
      @input_headers = qw (chr start stop ref_base var_base var_type gene_name transcript_id species transcript_source transcript_version strand transcript_status var_effect_type coding_pos          aa_change ucsc_cons domain all_domains deletion_substructures transcript_error default_gene_name gene_name_source ensembl_gene_id);
    }else{
      my $col_count = 0;
      open (TMP, $new_annotated_file) || die "\n\nCould not open variant annotation file: $new_annotated_file\n\n";
      while(<TMP>){
      chomp($_);
      my @line = split("\t", $_);
      $col_count = scalar(@line);
      }
      close(TMP);

      if ($col_count == 24){
        @input_headers = qw (chr start stop ref_base var_base var_type gene_name transcript_id species transcript_source transcript_version strand transcript_status var_effect_type coding_pos aa_change ucsc_cons domain all_domains deletion_substructures transcript_error default_gene_name gene_name_source ensembl_gene_id);
      }else{
        $self->error_message("Unexpected column count ($col_count) found in SNV/INDEL file - ClinSeq is not compatible with old annotation format...");
        die();
      }
    }

    #Get AA changes from full .annotated file
    my %aa_changes;
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
      headers => \@input_headers,
      input => "$new_annotated_file",
      separator => "\t",
    );
    while (my $data = $reader->next) {
      my $coord = $data->{chr} .':'. $data->{start} .'-'. $data->{stop};
      $data->{coord} = $coord;
      #Apply the AA effect filter
      unless ($data->{var_effect_type} =~ /$aa_effect_filter/){
        next();
      }
      #Apply the MT/chrM filter
      if ($filter_mt){
        my $chr = $data->{chr};
        if ($chr =~ /^MT$|^chrMT$|^M$|^chrM$/i){
          next();
        }
      }
      $aa_changes{$coord}{$data->{aa_change}}=1;
    }

    #Get compact SNV info from the '.top' file but grab the complete list of AA changes from the '.annotated' file
    $reader = Genome::Utility::IO::SeparatedValueReader->create(
      headers => \@input_headers,
      input => "$new_annotated_top_file",
      separator => "\t",
    );

    while (my $data = $reader->next){
      my $coord = $data->{chr} .':'. $data->{start} .'-'. $data->{stop};
      #Apply the AA effect filter
      unless ($data->{var_effect_type} =~ /$aa_effect_filter/){
        next();
      }
      #Apply the MT/chrM filter
      if ($filter_mt){
        my $chr = $data->{chr};
        if ($chr =~ /^MT$|^chrMT$|^M$|^chrM$/i){
          next();
        }
      }
      my %aa = %{$aa_changes{$coord}};
      my $aa_string = join(",", sort keys %aa);
      $data_out{$coord}{gene_name} = $data->{gene_name};
      $data_merge{$var_type}{$coord}{gene_name} = $data->{gene_name};
      $data_out{$coord}{ensembl_gene_id} = $data->{ensembl_gene_id};
      $data_merge{$var_type}{$coord}{ensembl_gene_id} = $data->{ensembl_gene_id};
      $data_out{$coord}{aa_changes} = $aa_string;
      $data_merge{$var_type}{$coord}{aa_changes} = $aa_string;
      $data_out{$coord}{ref_base} = $data->{ref_base};
      $data_merge{$var_type}{$coord}{ref_base} = $data->{ref_base};
      $data_out{$coord}{var_base} = $data->{var_base};
      $data_merge{$var_type}{$coord}{var_base} = $data->{var_base};

      #Make note of the datatype (wgs or exome) this variant was called by...
      if ($data_type eq "wgs"){
        $data_merge{$var_type}{$coord}{wgs} = 1;
      }
      if ($data_type eq "exome"){
        $data_merge{$var_type}{$coord}{exome} = 1;
      }

      #Attempt to fix the gene name:
      my $fixed_gene_name = &fixGeneName('-gene'=>$data->{gene_name}, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-verbose'=>0);
      $data_out{$coord}{mapped_gene_name} = $fixed_gene_name;
      $data_merge{$var_type}{$coord}{mapped_gene_name} = $fixed_gene_name;
    }

    #Print out the resulting list, sorting on fixed gene name
    open (OUT, ">$compact_file") || die "\n\nCould not open output file: $compact_file\n\n";
    print OUT "coord\tgene_name\tmapped_gene_name\tensembl_gene_id\taa_changes\tref_base\tvar_base\n";
    foreach my $coord (sort {$data_out{$a}->{mapped_gene_name} cmp $data_out{$b}->{mapped_gene_name}} keys %data_out){
      print OUT "$coord\t$data_out{$coord}{gene_name}\t$data_out{$coord}{mapped_gene_name}\t$data_out{$coord}{ensembl_gene_id}\t$data_out{$coord}{aa_changes}\t$data_out{$coord}{ref_base}\t$data_out{$coord}{var_base}\n";
    }
    close(OUT);

    #Store the path for this output file as an output of the command ...
    $self->wgs_snv_file($compact_file) if ($data_type eq "wgs" && $var_type eq "snv");
    $self->exome_snv_file($compact_file) if ($data_type eq "exome" && $var_type eq "snv");
    $self->wgs_indel_file($compact_file) if ($data_type eq "wgs" && $var_type eq "indel");
    $self->exome_indel_file($compact_file) if ($data_type eq "exome" && $var_type eq "indel");
  }

  #If both WGS and Exome data were present, print out a data merge for SNVs and Indels
  if ($self->wgs_build && $self->exome_build){
    my $snv_wgs_exome_dir = $snv_dir . "wgs_exome/"; 
    Genome::Sys->create_directory($snv_wgs_exome_dir);
    my $indel_wgs_exome_dir = $indel_dir . "wgs_exome/";
    Genome::Sys->create_directory($indel_wgs_exome_dir);

    my $snv_merge_file = "$snv_wgs_exome_dir"."snvs.hq.tier1.v1.annotated.compact.tsv";
    my $indel_merge_file = "$indel_wgs_exome_dir"."indels.hq.tier1.v1.annotated.compact.tsv";

    open (OUT, ">$snv_merge_file") || die "\n\nCould not open output file: $snv_merge_file\n\n";
    print OUT "coord\tgene_name\tmapped_gene_name\tensembl_gene_id\taa_changes\tref_base\tvar_base\twgs_called\texome_called\n";
    my %data_out = %{$data_merge{'snv'}};
    foreach my $coord (sort {$data_out{$a}->{mapped_gene_name} cmp $data_out{$b}->{mapped_gene_name}} keys %data_out){
      my $wgs_called = 0;
      if (defined($data_out{$coord}{wgs})){ $wgs_called = 1; }
      my $exome_called = 0;
      if (defined($data_out{$coord}{exome})){ $exome_called = 1; }
      print OUT "$coord\t$data_out{$coord}{gene_name}\t$data_out{$coord}{mapped_gene_name}\t$data_out{$coord}{ensembl_gene_id}\t$data_out{$coord}{aa_changes}\t$data_out{$coord}{ref_base}\t$data_out{$coord}{var_base}\t$wgs_called\t$exome_called\n";
    }
    close(OUT);
    $self->wgs_exome_snv_file($snv_merge_file);
    
    open (OUT, ">$indel_merge_file") || die "\n\nCould not open output file: $indel_merge_file\n\n";
    print OUT "coord\tgene_name\tmapped_gene_name\tensembl_gene_id\taa_changes\tref_base\tvar_base\twgs_called\texome_called\n";
    %data_out = %{$data_merge{'indel'}};
    foreach my $coord (sort {$data_out{$a}->{mapped_gene_name} cmp $data_out{$b}->{mapped_gene_name}} keys %data_out){
      my $wgs_called = 0;
      if (defined($data_out{$coord}{wgs})){ $wgs_called = 1; }
      my $exome_called = 0;
      if (defined($data_out{$coord}{exome})){ $exome_called = 1; }
      print OUT "$coord\t$data_out{$coord}{gene_name}\t$data_out{$coord}{mapped_gene_name}\t$data_out{$coord}{ensembl_gene_id}\t$data_out{$coord}{aa_changes}\t$data_out{$coord}{ref_base}\t$data_out{$coord}{var_base}\t$wgs_called\t$exome_called\n";
    }
    close(OUT);
    $self->wgs_exome_indel_file($indel_merge_file);
  }

  return 1;
}

1;


