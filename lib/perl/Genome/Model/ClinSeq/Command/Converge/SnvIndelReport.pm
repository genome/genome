package Genome::Model::ClinSeq::Command::Converge::SnvIndelReport;
use strict;
use warnings;
use Genome;
use Data::Dumper;
use List::MoreUtils qw/ uniq /;
use Genome::Info::IUB;
use Spreadsheet::WriteExcel;

class Genome::Model::ClinSeq::Command::Converge::SnvIndelReport {
    is => ['Genome::Model::ClinSeq::Command::Converge::Base',
           'Genome::Model::ClinSeq::Util'],
    has_input => [
        outdir => {
               is => 'FilesystemPath',
               doc => 'Directory where output files will be written',
        },
    ],
    has_optional_input => [
        target_gene_list => {
                is => 'FilesystemPath',
                doc => 'Genes of interest to be highlighted (e.g. AML RMG list).  Tab delimited file with ENSGs in first column',
                default => '/gscmnt/sata132/techd/mgriffit/aml_trios/AML_RMG_List-Ensembl_Gene_IDs.tsv',
        },
        target_gene_list_name => {
                is => 'Text',
                doc => 'Human readable name used to refer to the target gene list',
                default => 'AML_RMG',
        },
        variant_filter_list => {
                is => 'FileSystemPath',
                doc => 'Arbitrary list of variants to filter out.  In "chr start stop ref var" format (space or tab separated). File may have a header or not. Indels should be indicated with a "-" as done in the output of this script.',
        },
        tiers => {
                is => 'Text',
                default => 'tier1,tier2,tier3,tier4',
                doc => 'Which tiers of snvs and indels to allow',
        },
        _annotation_build_id => {
              is => 'Text',
              default => 'd00a39c84382427fa0efdec3229e8f5f',
        },
        annotation_build => {
              is => 'Genome::Model::Build',
              id_by => '_annotation_build_id',
              doc => 'Desired reference annotation build',
        },
        max_normal_vaf => {
              is => 'Number',
              default => 10,
              doc => 'Variants with a normal VAF greater than this (in any normal sample) will be filtered out'
        },
        min_tumor_var_count => {
              is => 'Number',
              default => 3,
              doc => 'Variants with a tumor var count less than this (in any tumor sample) will be filtered out',
        },
        min_tumor_vaf => {
              is => 'Number',
              default => 2.5,
              doc => 'Variants with a tumor VAF less than this (in any tumor sample) will be filtered out',
        },
        min_coverage => {
              is => 'Number',
              default => 20,
              doc => 'Variants with coverage less than this (in any sample) will be filtered out',
        },
        max_gmaf => {
              is => 'Number',
              default => 2.5,
              doc => 'Variants with a GMAF from 1000 genomes higher than this will be filtered out',
        },
        trv_type_filter => {
              is => 'Text',
              default => '-|3_prime_flanking_region|3_prime_untranslated_region|5_prime_flanking_region|5_prime_untranslated_region|intronic|silent|rna|splice_region',
              doc => 'Variants with these transcript variant types will be filtered out of the final coding result.',
        },
        min_reads_per_lib => {
              is => 'Number',
              default => 1,
              doc => 'In per library read counting, the minimum variant supporting reads for a library to be said to support a variant',
        },
        min_tumor_var_supporting_libs => {
              is => 'Number',
              default => 0,
              doc => 'In per library analysis, variants with less than this number of tumor libraries supporting will be filtered out',
        },
        min_snv_caller_count => {
              is => 'Number',
              default => 1,
              doc => 'SNV variants called by fewer this number of SNV variant callers will be filtered out',
        },
        min_indel_caller_count => {
              is => 'Number',
              default => 1,
              doc => 'INDEL variants called by fewer this number of INDEL variant callers will be filtered out',
        },
        per_library => {
              is => 'Boolean',
              doc => 'Do per library bam-readcounting and generate associated statistics, summaries and figures'
        },
        summarize => {
              is => 'Boolean',
              doc => 'Summarize the SnvIndel report.',
              default => 1,
        },
        test => {
              is => 'Number',
              doc => 'Only import this many variants (for testing purposes)',
        },
        chromosome => {
              is => 'Text',
              doc => 'Limit analysis to variants on this chromosome only',
        },
        clean => {
              is => 'Boolean',
              doc => 'Remove intermediate files from output if you are running in a pipeline',
        },
        tmp_space => {
              is => 'Boolean',
              doc => 'Perform all file I/O in /tmp and copy results to outdir at the end',
        },
        _wgs_snv_variant_sources_file => {
              is => 'FilesystemPath',
        },
        _exome_snv_variant_sources_file => {
              is => 'FilesystemPath',
        },
        _wgs_indel_variant_sources_file => {
              is => 'FilesystemPath',
        },
        _exome_indel_variant_sources_file => {
              is => 'FilesystemPath',
        },
    ],
    has_param => [
        lsf_resource => {
            value => q{-R 'select[mem>12000] rusage[mem=12000]' -M 12000000},
        },
    ],
    doc => 'converge SNV and InDels from multiple clin-seq builds, annotate, bam-readcount, etc. and summarize into a single spreadsheet',
};

sub help_synopsis {
  return <<EOS

genome model clin-seq converge snv-indel-report --builds='id in ["4b7539bb10cc4b9c97577cf11f4c79a2","cdca0edf526c4fe193d3054627a5871b"]' --outdir=/tmp/snv_indel_report/

genome model clin-seq converge snv-indel-report --builds='model.model_groups.id=9d0fcdca2b5d4f4385b83d2f75addac4,is_last_complete=1' --outdir=/tmp/snv_indel_report/

genome model clin-seq converge snv-indel-report --builds='model_groups.id=9d0fcdca2b5d4f4385b83d2f75addac4,is_last_complete=1' --outdir=/tmp/snv_indel_report/

genome model clin-seq converge snv-indel-report --builds='model.id in ["279f50e35d2b479ea3c32486eafd4ad4","7143119a93984056ae3f32c88c9ac2a1"],is_last_complete=1' --outdir=/tmp/snv_indel_report/

EOS
}

sub help_detail {
  return <<EOS

Create a summary spreadsheet of SNVs and Indels for a set of clin-seq builds (e.g. an AML normal, day0_tumor, day30_tumor trio)

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

  if ($self->target_gene_list){
    unless (-e $self->target_gene_list) {
      push @errors, UR::Object::Tag->create(
                                            type => 'error',
                                            properties => ['target_gene_list'],
                                            desc => "File: " . $self->target_gene_list . " not found",
                                          );
    }
  }
  if ($self->tiers){
    my @tiers = split(",", $self->tiers);
    my $error = 0;
    $error++ unless (scalar(@tiers) >= 1);
    foreach my $tier (@tiers){
      $error++ unless ($tier =~ /^tier\d+$/);
    }
    if ($error){
      push @errors, UR::Object::Tag->create(
                                              type => 'error',
                                              properties => ['tiers'],
                                              desc => "Tiers: " . $self->tiers . " has incorrect format, see default",
                                            );
    }
  }
  if ($self->variant_filter_list){
    unless (-e $self->variant_filter_list) {
      push @errors, UR::Object::Tag->create(
                                              type => 'error',
                                              properties => ['variant_filter_list'],
                                              desc => "File: " . $self->variant_filter_list . " not found",

                                           );
    }
  }

  return @errors;
}


sub execute {
  my $self = shift;

  #Add trailing '/' to outdir if needed
  unless ($self->outdir =~ /\/$/){
    my $outdir = $self->outdir . "/";
    $self->outdir($outdir);
  }

  #If the user specified, perform all file I/O in /tmp and copy result over at the end
  my $original_outdir = $self->outdir;
  if ($self->tmp_space){
    my $outdir = Genome::Sys->create_temp_directory;
    $outdir .= "/" unless ($outdir =~ /\/$/);
    $self->outdir($outdir);
  }

  #Get the case name for this set of builds, if more than one are found, warn the user
  my $case_name = $self->get_case_name;
  $self->status_message("Producing report for individual: $case_name");

  #Gather variants for the tiers specified by the user from each build. Note which build each came from.
  #Get these from the underlying somatic-variation builds.
  #Annotate all variants (gmt annotate transcript-variants --help)
  my @clinseq_builds = $self->builds;
  my (%somatic_builds, %rnaseq_builds);
  foreach my $clinseq_build(@clinseq_builds) {
    $clinseq_build->resolve_somatic_builds(\%somatic_builds);
    $clinseq_build->resolve_rnaseq_builds(\%rnaseq_builds);
  }
  my $annotation_build_name = $self->annotation_build->name;
  $self->status_message("Using annotation build: $annotation_build_name");
  my $bed_dir = $self->outdir . "bed_files/";
  my $result = $self->gather_variants('-somatic_builds'=>\%somatic_builds, '-bed_dir'=>$bed_dir);
  my $variants = $result->{'variants'};
  my $header = $result->{'header'};

  #If no variants were found, warn the user and end here end here
  unless (keys %{$variants}){
    my $rm_cmd = "rm -fr $bed_dir";
    $self->warning_message("no variants found, no snv-indel report will be created");
    Genome::Sys->shellcmd(cmd => $rm_cmd);
    exit 1;
  }

  #Make note of which variants correspond to $self->target_gene_list
  $self->intersect_target_gene_list('-variants'=>$variants);

  #Pull in source variant caller info from the clin-seq build results
  $self->get_variant_caller_sources('-variants'=>$variants);

  #Print a grand table of all annotated variants with extra annotation columns (tier, target gene list)
  my $grand_anno_file = $self->print_grand_anno_table('-variants'=>$variants, '-header'=>$header);

  #TODO: Annotate with information from COSMIC


  #Identify underlying reference-alignments, sample names, sample common names, and timepoints (if available)
  my $align_builds = $self->get_ref_align_builds(
    '-somatic_builds'=>\%somatic_builds,
    '-rnaseq_builds'=>\%rnaseq_builds);
  my @prefixes = $self->get_header_prefixes('-align_builds'=>$align_builds);
  #Get bam-readcounts for all positions for all BAM files
  my $grand_anno_count_file = $self->add_read_counts(
    '-align_builds'=>$align_builds,
    '-anno_file'=>$grand_anno_file,
    '-prefixes'=>\@prefixes);

  #Parse the BAM read count info and gather minimal data needed to apply filters:
  #Max normal VAF, Min Coverage
  $self->parse_read_counts('-align_builds'=>$align_builds, '-grand_anno_count_file'=>$grand_anno_count_file, '-variants'=>$variants);

  #Get per-library bam-readcounts for all positions for all BAM files
  my $grand_anno_per_lib_count_file;
  if ($self->per_library){
    $grand_anno_per_lib_count_file = $self->add_per_lib_read_counts('-align_builds'=>$align_builds, '-anno_file'=>$grand_anno_count_file);
  }
  #Parse the per lib BAM read count info
  my $per_lib_header;
  if ($self->per_library){
    $per_lib_header = $self->parse_per_lib_read_counts('-align_builds'=>$align_builds, '-grand_anno_per_lib_count_file'=>$grand_anno_per_lib_count_file, '-variants'=>$variants);
  }

  #Apply arbitrary variant filter list (a list of variants supplied in a file that are to be removed)
  if ($self->variant_filter_list){
    $self->apply_filter_list('-variants'=>$variants);
  }

  #Apply some automatic filters:
  $self->apply_variant_filters('-variants'=>$variants);


  #TODO: Make note of which variants lie within a particular set of regions of interest (e.g. nimblegen v3 + AML RMG).  Have option to filter these out
  # - Allow for wingspan to be added to these regions

  #TODO: Filter out variants falling within known false positive regions/genes?

  #TODO: Limit analysis to variants *called in* a particular tumor only (e.g. day0 tumor) - rather than taking the union of calls from all samples

  #TODO: Add an additional filter that uses the false postive filter: 'gmt validation identify-outliers'

  #TODO: Add additional filters for low VAF variants (e.g. if VAF < 10%, require at least 3 callers, support from multiple libraries, etc.)

  #TODO: Get Exome Sequencing Project MAF and filter on this value as well (use existing gmt for this)


  #Write out final tsv files (filtered and unfiltered), a clean version with useless columns removed, and an Excel spreadsheet version of the final file
  my $result_files = $self->print_final_files('-variants'=>$variants, '-grand_anno_count_file'=>$grand_anno_count_file, '-case_name'=>$case_name, '-align_builds'=>$align_builds, '-per_lib_header'=>$per_lib_header);

  #Also for certain sample combinations create custom visualizations using R
  $self->create_plots('-result_files'=>$result_files, '-align_builds'=>$align_builds, '-case_name'=>$case_name);

  #TODO: If there are two tumors being converged, Produce SciClone plots showing tumor1 vs. tumor2 VAF and with gene names labelled


  #Print out a table of subject names for reference
  my $subject_table_file = $self->print_subject_table('-align_builds'=>$align_builds);


  #If the user specified the --clean option, remove intermediate files from the output dir
  if ($self->clean){
    my $rm_cmd = "rm -fr $bed_dir";
    Genome::Sys->shellcmd(cmd => $rm_cmd);
  }

  #If the user specified, perform all file I/O in /tmp and copy result over at the end
  if ($self->tmp_space){
    my $cp_cmd = "cp -fr " . $self->outdir . "* $original_outdir";
    Genome::Sys->shellcmd(cmd => $cp_cmd);
  }

  if($self->summarize) {
    my $summarize = Genome::Model::ClinSeq::Command::Converge::SummarizeSnvIndelReport->create(
      outdir => $original_outdir,
      min_mq => $self->min_quality_score,
      min_bq => $self->min_base_quality,
      filtered_report => $result_files->{final_filtered_clean_tsv},
      unfiltered_report => $result_files->{final_unfiltered_clean_tsv},
    );
    $summarize->execute();
  }

  return 1;
};


sub fixIUB{
    my ($ref,$var) = @_;
    my @vars = Genome::Info::IUB->variant_alleles_for_iub($ref,$var);
    return @vars;
}


sub print_subject_table{
  my $self = shift;
  my %args = @_;
  my $align_builds = $args{'-align_builds'};

  my $outfile = $self->outdir . "subjects_legend.txt";
  my $out_fh = Genome::Sys->open_file_for_writing($outfile);
  print $out_fh "name\tprefix\tday\ttimepoint_position\tsample_type\torder\ttissue_desc\ttissue_label\textraction_type\textraction_label\n";
  foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys %{$align_builds}){
    my $order = $align_builds->{$name}->{order};
    my $prefix = $align_builds->{$name}->{prefix};
    my $day = $align_builds->{$name}->{day};
    my $timepoint_position = $align_builds->{$name}->{timepoint_position};
    my $sample_type = $align_builds->{$name}->{sample_common_name};
    my $tissue_desc = $align_builds->{$name}->{tissue_desc};
    my $tissue_label = $align_builds->{$name}->{tissue_label};
    my $extraction_type = $align_builds->{$name}->{extraction_type};
    my $extraction_label = $align_builds->{$name}->{extraction_label};
    print $out_fh "$name\t$prefix\t$day\t$timepoint_position\t$sample_type\t$order\t$tissue_desc\t$tissue_label\t$extraction_type\t$extraction_label\n";
  }
  close($out_fh);

  return $outfile;
}


sub gather_variants{
  my $self = shift;
  my %args = @_;
  my $somatic_builds = $args{'-somatic_builds'};
  my $bed_dir = $args{'-bed_dir'};

  my $tiers_list = $self->tiers;
  my @tiers = split(",", $tiers_list);

  mkdir ($bed_dir) unless (-e $bed_dir && -d $bed_dir);

  #Get bed files to be processed
  my %bed_files;
  foreach my $somatic_build_id(sort keys %{$somatic_builds}){
    my $somatic_build = $somatic_builds->{$somatic_build_id}->{build};
    my $somatic_build_type = $somatic_builds->{$somatic_build_id}->{type};
    my $somatic_build_dir = $somatic_build->data_directory;

    foreach my $tier (@tiers){
      my $snvs_file = "$somatic_build_dir/effects/snvs.hq.novel.$tier.v2.bed";
      my $indels_file = "$somatic_build_dir/effects/indels.hq.novel.$tier.v2.bed";
      unless (-e $snvs_file && -e $indels_file){
        die $self->error_message("Could not find expected file:\n$snvs_file\n$indels_file");
      }
      $bed_files{$snvs_file}{somatic_build_id} = $somatic_build_id;
      $bed_files{$snvs_file}{var_type} = "snv";
      $bed_files{$snvs_file}{data_type} = $somatic_build_type;
      $bed_files{$snvs_file}{tier} = $tier;
      $bed_files{$snvs_file}{vcf_file} = "$somatic_build_dir/variants/snvs.annotated.vcf.gz";
      $bed_files{$indels_file}{somatic_build_id} = $somatic_build_id;
      $bed_files{$indels_file}{var_type} = "indel";
      $bed_files{$indels_file}{data_type} = $somatic_build_type;
      $bed_files{$indels_file}{tier} = $tier;
      $bed_files{$indels_file}{vcf_file} = "$somatic_build_dir/variants/indels.detailed.vcf.gz";
    }
  }

  #Make a copy of each desired .bed file.  Perform basic clean-ups and sanity checks
  foreach my $file (sort keys %bed_files){
    my $c = 0;
    my $var_type = $bed_files{$file}{var_type};
    my $data_type = $bed_files{$file}{data_type};
    my $tier = $bed_files{$file}{tier};
    my $somatic_build_id = $bed_files{$file}{somatic_build_id};
    my $new_file = $bed_dir . "$somatic_build_id" . "_$data_type" . "_$var_type" . "_$tier" . "clean.bed";
    $self->status_message("\nProcessing $var_type $data_type $tier file: $file");

    if (-e $new_file){
      $self->warning_message("Using pre-generated file: $new_file");
    }else{
      my $var_fh = Genome::Sys->open_file_for_reading($file);
      my $new_fh = Genome::Sys->open_file_for_writing($new_file);
      while(<$var_fh>){
        $c++;
        chomp($_);
        my @line = split("\t", $_);
        my ($chr, $start, $end, $var)  = ($line[0], $line[1], $line[2], $line[3]);
        if ($self->chromosome){next unless ($chr eq $self->chromosome);}
        #Make indel format consistent
        $var =~ s/\*/0/g;
        print $new_fh "$chr\t$start\t$end\t$var\n";
        if ($self->test){last if $c > $self->test;}
      }
      close($var_fh);
      close($new_fh);
    }
    $bed_files{$file}{clean_bed_file} = $new_file;

    #Run the annnotator on each bed file and add RSIDs and GMAFs for all variants
    #gmt annotate transcript-variants --variant-bed-file='' --output-file='' --annotation-filter=top --reference-transcripts='NCBI-human.ensembl/74_37'
    my $anno_file = $bed_files{$file}{clean_bed_file} . ".anno";
    if (-e $anno_file){
      $self->warning_message("Using pre-generated file: $anno_file");
    }else{
      if (-s $bed_files{$file}{clean_bed_file}){
        my $annotate_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
              variant_bed_file=>$bed_files{$file}{clean_bed_file},
              output_file=>$anno_file,
              annotation_filter=>'top',
              reference_transcripts=>$self->annotation_build->name,
            );
        my $r1 = $annotate_cmd->execute();
        die $self->error_message("annotation cmd unsuccessful") unless ($r1);
      }else{
        my $touch_cmd = "touch $anno_file";
        Genome::Sys->shellcmd(cmd => $touch_cmd);
      }
    }
    $bed_files{$file}{anno_file} = $anno_file;

    #Get RSIDs and GMAFs for all variants
    #gmt annotate add-rsid --anno-file='' --output-file='' --vcf-file=''
    my $vcf_file = $bed_files{$file}{vcf_file};
    my $rsid_file = $anno_file . ".rsid";

    if (-e $rsid_file){
      $self->warning_message("Using pre-generated file: $rsid_file");
    }else{
      my $rsid_cmd = Genome::Model::Tools::Annotate::AddRsid->create(
            anno_file=>$anno_file,
            vcf_file=>$vcf_file,
            output_file=>$rsid_file,
          );
      my $r2 = $rsid_cmd->execute();
      die $self->error_message("rsid cmd unsuccessful") unless ($r2);
    }
    $bed_files{$file}{rsid_file} = $rsid_file;
  }

  #Parse all variants into a single hash (keyed on $chr_$start_$end_$ref_$var)
  my $header;
  my %headers;
  my %variants;

  foreach my $file (sort keys %bed_files){
    my $tier = $bed_files{$file}{tier};
    my $rsid_file = $bed_files{$file}{rsid_file};
    my $data_type = $bed_files{$file}{data_type};
    my $var_fh = Genome::Sys->open_file_for_reading($rsid_file);
    while(<$var_fh>){
      chomp($_);
      if ($_ =~ /^chromosome\_name/){
        $header = $_;
        $headers{$header}=1;
        next;
      }
      my @line = split("\t", $_);
      my ($chr, $start, $stop, $ref, $var) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
      my $variant_type = $line[5];
      my $trv_type = $line[13];
      my $ensembl_gene_id = $line[23];
      my $v = $chr . "_$start" . "_$stop" . "_$ref" . "_$var";
      $variants{$v}{anno_line} = $_;
      $variants{$v}{tier} = $tier;
      $variants{$v}{ensembl_gene_id} = $ensembl_gene_id;
      $variants{$v}{variant_type} = $variant_type;
      $variants{$v}{trv_type} = $trv_type;
      $variants{$v}{filtered} = "";
      if(defined $variants{$v}{data_type}) {
        unless ($variants{$v}{data_type} =~ /$data_type/) {
          $variants{$v}{data_type} .= $data_type . ",";
        }
      } else {
        $variants{$v}{data_type} = $data_type . ",";
      }
    }
    close($var_fh);
  }

  #Make sure all headers came out identical
  my $header_count = keys %headers;
  die $self->error_message("found more than one different header in anno files to be merged") if $header_count > 1;

  my %result;
  $result{'variants'} = \%variants;
  $result{'header'} = $header;

  return(\%result);
}


sub print_grand_anno_table{
  my $self = shift;
  my %args = @_;
  my $variants = $args{'-variants'};
  my $header = $args{'-header'};

  #Also produce a joined annotated file with a single header and all snvs and indels combined
  my $grand_anno_file = $self->outdir . "variants.all.anno";
  if (-e $grand_anno_file){
    $self->warning_message("using pre-generated file: $grand_anno_file");
  }else{
    my $anno_fh = Genome::Sys->open_file_for_writing($grand_anno_file);
    my $target_gene_list_name = $self->target_gene_list_name;
    my $new_header = "$header\ttier\t$target_gene_list_name";
    print $anno_fh "$new_header\n";
    foreach my $v (sort keys %{$variants}){
      print $anno_fh "$variants->{$v}->{anno_line}\t$variants->{$v}->{tier}\t$variants->{$v}->{target_gene_list_match}\n";
    }
    close($anno_fh);
  }

  #delete the anno line now that it is not needed
  foreach my $v (keys %{$variants}){
    delete $variants->{$v}->{anno_line};
  }

  return($grand_anno_file);
}


sub intersect_target_gene_list{
  my $self = shift;
  my %args = @_;
  my $variants = $args{'-variants'};

  my $target_gene_list = $self->target_gene_list;
  my $target_gene_list_name = $self->target_gene_list_name;

  my %genes;
  my $genes_fh = Genome::Sys->open_file_for_reading($target_gene_list);
  while(<$genes_fh>){
    chomp($_);
    my @line = split("\t", $_);
    my $ensg = $line[0];
    $genes{$ensg}{name} = $line[1];
  }
  close($genes_fh);

  foreach my $v (keys %{$variants}){
    my $ensembl_gene_id = $variants->{$v}->{ensembl_gene_id};
    $variants->{$v}->{target_gene_list_match} = 0;
    if ($genes{$ensembl_gene_id}){
      $variants->{$v}->{target_gene_list_match} = 1;
    }
  }

  return;
}


sub get_variant_caller_sources{
  my $self = shift;
  my %args = @_;
  my $variants = $args{'-variants'};
  my @builds = $self->builds;

  #Initialize all variants with a list of source callers
  foreach my $v (keys %{$variants}){
    my %tmp;
    $variants->{$v}->{variant_source_callers} = \%tmp;
  }

  my @files;
  if ($self->_wgs_snv_variant_sources_file || $self->_exome_snv_variant_sources_file || $self->_wgs_indel_variant_sources_file || $self->_exome_indel_variant_sources_file){
    die $self->error_message("Do not specify variant source caller files directly if you are processing multiple clin-seq builds") if (scalar(@builds) > 1);
    push(@files, $self->_wgs_snv_variant_sources_file) if ($self->_wgs_snv_variant_sources_file);
    push(@files, $self->_exome_snv_variant_sources_file) if ($self->_exome_snv_variant_sources_file);
    push(@files, $self->_wgs_indel_variant_sources_file) if ($self->_wgs_indel_variant_sources_file);
    push(@files, $self->_exome_indel_variant_sources_file) if ($self->_exome_indel_variant_sources_file);
  }else{
    foreach my $build (@builds){
      my $common_name = $build->common_name;
      my $data_dir = $build->data_directory;
      my $base_path = $data_dir . "/" . $common_name . "/variant_source_callers";
      my $exome_base_path = $base_path . "/exome/";
      my $wgs_base_path = $base_path . "/wgs/";
      my @files_list = ($exome_base_path."snv_sources.tsv", $wgs_base_path."snv_sources.tsv", $exome_base_path."indel_sources.tsv", $wgs_base_path."indel_sources.tsv");
      push(@files, @files_list);
    }
  }

  my %columns;
  foreach my $file (@files){
    $self->status_message("processing variant caller sources file: $file");
    next unless (-e $file);
    my $source_fh = Genome::Sys->open_file_for_reading($file);
    while(<$source_fh>){
      chomp($_);
      my @line = split("\t", $_);
      if ($_ =~ /^coord/){
        my $p = 0;
        foreach my $col (@line){
          $columns{$col}{p} = $p;
          $p++;
        }
        next;
      }
      my ($chr, $start, $end, $variant, $caller_string, $tier) = ($line[$columns{'chr'}{p}], $line[$columns{'start'}{p}], $line[$columns{'end'}{p}], $line[$columns{'variant'}{p}], $line[$columns{'callers'}{p}], $line[$columns{'tier'}{p}]);
      $variant =~ s/\*/\-/g;
      $variant =~ s/0/\-/g;
      my ($ref, $var) = split(/\//, $variant);
      my @caller_list;
      if ($caller_string){
        @caller_list = split(",", $caller_string);
      }

      #Convert bed to anno
      if ($variant =~ /^[\-0\*]/){ #indel INS
        $end = $end+1;
      } else { #indel DEL or SNV
        $start = $start+1;
      }

      if ($file =~ /snv\_sources/){
        my @vars = fixIUB($ref, $var);
        foreach my $var_fixed (@vars){
          my $v = "$chr"."_$start"."_$end"."_$ref"."_"."$var_fixed";
          if (defined($variants->{$v})){
            my $callers = $variants->{$v}->{variant_source_callers};
            foreach my $caller (@caller_list){
              $callers->{$caller}=1;
            }
          }else{
            die $self->error_message("could not match variant to hash: $v $chr $tier") unless (($self->test) || !($self->tiers =~ /$tier/) || ($self->chromosome && ($chr ne $self->chromosome)));
          }
        }
      }elsif ($file =~ /indel\_sources/){
        my $v = "$chr"."_$start"."_$end"."_$ref"."_"."$var";
        if (defined($variants->{$v})){
          my $callers = $variants->{$v}->{variant_source_callers};
          foreach my $caller (@caller_list){
            $callers->{$caller}=1;
          }
        }else{
          die $self->error_message("could not match variant to hash: $v $chr $tier") unless (($self->test) || !($self->tiers =~ /$tier/) || ($self->chromosome && ($chr ne $self->chromosome)));
        }
      }else{
        die $self->error_message("file not recognized as snv or indel: $file");
      }
    }
    close ($source_fh);
  }

  #Initialize all variants with a list of source callers
  foreach my $v (keys %{$variants}){
    my %callers = %{$variants->{$v}->{variant_source_callers}};
    my @caller_list = keys %callers;
    my $caller_count = scalar(@caller_list);
    my $caller_string = join(",", @caller_list);
    if ($caller_string){
      $variants->{$v}->{variant_source_callers} = $caller_string;
      $variants->{$v}->{variant_source_caller_count} = $caller_count;
    }else{
      $variants->{$v}->{variant_source_callers} = "NA";
      $variants->{$v}->{variant_source_caller_count} = "NA";
    }
  }

  return;
}


sub add_per_lib_read_counts{
  my $self = shift;
  my %args = @_;
  my $align_builds = $args{'-align_builds'};
  my $grand_anno_file = $args{'-anno_file'};

  my @bam_files;
  foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
    push(@bam_files, $align_builds->{$name}->{bam_path});
  }
  my $bam_list = join(",", @bam_files);

  #Get the reference fasta
  my $reference_build = $self->resolve_clinseq_reference_build;
  my $reference_fasta = $reference_build->full_consensus_path('fa');

  #Determine header prefixes to use. In order of preference if all are unique: (time_points, samples, names)
  my @prefixes;
  foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
    push(@prefixes, $align_builds->{$name}->{prefix});
  }
  my $header_prefixes = join(",", @prefixes);

  my $output_file = $self->outdir . "/variants.all.anno.per.library.readcounts";
  if (-e $output_file){
    $self->warning_message("using pre-generated per-library bam read-count file: $output_file");
  }else{
    my $add_count_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
            bam_files=>$bam_list,
            genome_build=>$reference_fasta,
            output_file=>$output_file,
            variant_file=>$grand_anno_file,
            header_prefixes=>$header_prefixes,
            per_library=>1,
            bam_readcount_version => $self->bam_readcount_version,
          );
    my $r = $add_count_cmd->execute();
    die $self->error_message("per-lane add-readcounts cmd unsuccessful") unless ($r);
  }

  return ($output_file);
}


sub parse_read_counts{
  my $self = shift;
  my %args = @_;
  my $align_builds = $args{'-align_builds'};
  my $grand_anno_count_file = $args{'-grand_anno_count_file'};
  my $variants = $args{'-variants'};

  my %columns;
  my $l = 0;
  my $var_fh = Genome::Sys->open_file_for_reading($grand_anno_count_file);
  while(<$var_fh>){
    chomp($_);
    my @line = split("\t", $_);
    if ($l == 0){
      my $c = 0;
      foreach my $column (@line){
        $columns{$column}{c} = $c;
        $c++;
      }
      #Check for neccessary columns
      my $gmaf_colname = "GMAF";
      die $self->error_message("could not find expected column ($gmaf_colname) in file $grand_anno_count_file") unless($columns{$gmaf_colname});
      foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
        my $prefix = $align_builds->{$name}->{prefix};
        my $ref_count_colname = $prefix . "_ref_count";
        my $var_count_colname = $prefix . "_var_count";
        my $vaf_colname = $prefix . "_VAF";
        die $self->error_message("could not find expected columns ($ref_count_colname $var_count_colname $vaf_colname) in file $grand_anno_count_file") unless($columns{$ref_count_colname} && $columns{$var_count_colname} && $columns{$vaf_colname});
      }
      $l++;
      next;
    }
    $l++;

    my ($chr, $start, $stop, $ref, $var) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
    my $v = $chr . "_$start" . "_$stop" . "_$ref" . "_$var";
    die $self->error_message("parsed a variant that is not defined in the variant hash") unless $variants->{$v};

    #Get max normal VAF (across all 'normal' samples) and min coverage (across all samples)
    #If there are multiple data types for a single sample.  Only consider the highest coverage for that sample when determining $min_coverage_observed
    #For example if there is exome AND WGS data we want the min coverage criteria to applied to either of these (whichever is higher)
    my %samples;
    my $max_normal_vaf_observed = 0;
    my $max_tumor_var_count_observed = 0;
    my $max_tumor_vaf_observed = 0;
    my $na_found = 0;
    my @covs;
    foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
      my $prefix = $align_builds->{$name}->{prefix};
      my $sample_name = $align_builds->{$name}->{sample_name};
      my $sample_common_name = $align_builds->{$name}->{sample_common_name};
      my $ref_count_colname = $prefix . "_ref_count";
      my $var_count_colname = $prefix . "_var_count";
      my $vaf_colname = $prefix . "_VAF";

      my $vaf = "NA";
      $vaf = $line[$columns{$vaf_colname}{c}] if (defined($line[$columns{$vaf_colname}{c}]));

      if ($vaf eq "NA"){
        $na_found = 1;
        push(@covs, "NA");
        next;
      }
      if ($sample_common_name =~ /normal/i){
        my $normal_vaf = $line[$columns{$vaf_colname}{c}];
        $max_normal_vaf_observed = $normal_vaf if ($normal_vaf > $max_normal_vaf_observed);
      }else{
        my $tumor_var_count = $line[$columns{$var_count_colname}{c}];
        $max_tumor_var_count_observed = $tumor_var_count if ($tumor_var_count > $max_tumor_var_count_observed);
        my $tumor_vaf = $line[$columns{$vaf_colname}{c}];
        $max_tumor_vaf_observed = $tumor_vaf if ($tumor_vaf > $max_tumor_vaf_observed);
      }
      my $coverage = $line[$columns{$ref_count_colname}{c}] + $line[$columns{$var_count_colname}{c}];
      push(@covs, $coverage);
      if (defined($samples{$sample_name})){
        $samples{$sample_name}{coverage} = $coverage if ($coverage > $samples{$sample_name}{coverage});
      }else{
        $samples{$sample_name}{prefix} = $prefix;
        $samples{$sample_name}{coverage} = $coverage;
      }
    }
    my $min_coverage_observed = 'inf';
    foreach my $sample_name (keys %samples){
      my $prefix = $samples{$sample_name}{prefix};
      my $coverage = $samples{$sample_name}{coverage};
      #don't apply min_coverage on rnaseq, the transcript might not be expressed.
      $min_coverage_observed = $coverage if ($coverage < $min_coverage_observed and $prefix !~ /rnaseq/);
    }

    if ($na_found){
      $max_normal_vaf_observed = "NA";
      $max_tumor_vaf_observed ="NA";
      $min_coverage_observed = "NA";
    }

    $variants->{$v}->{max_normal_vaf_observed} = $max_normal_vaf_observed;
    $variants->{$v}->{max_tumor_var_count_observed} = $max_tumor_var_count_observed;
    $variants->{$v}->{max_tumor_vaf_observed} = $max_tumor_vaf_observed;
    $variants->{$v}->{min_coverage_observed} = $min_coverage_observed;
    $variants->{$v}->{coverages} = \@covs;

    #Get the GMAF
    my $gmaf_colname = "GMAF";
    my $gmaf_string = $line[$columns{$gmaf_colname}{c}];
    my $gmaf = "NA";
    if ($gmaf_string =~ /^GMAF\=(.*)$|^\-$/){
      $gmaf = $1 if defined($1);
    }else{
      die $self->error_message("invalid GMAF string: $gmaf_string");
    }
    $variants->{$v}->{gmaf} = $gmaf;
  }
  close($var_fh);

  return;
}

sub parse_per_lib_read_counts{
  my $self = shift;
  my %args = @_;
  my $align_builds = $args{'-align_builds'};
  my $grand_anno_per_lib_count_file = $args{'-grand_anno_per_lib_count_file'};
  my $variants = $args{'-variants'};

  foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
    my $prefix = $align_builds->{$name}->{prefix};
  }

  my @per_lib_header;
  my %columns;
  my %libs;
  my $l = 0;
  my $var_fh = Genome::Sys->open_file_for_reading($grand_anno_per_lib_count_file);
  while(<$var_fh>){
    chomp($_);
    my @line = split("\t", $_);
    my $tumor_var_supporting_libs = 0;
    if ($l == 0){
      my $c = 0;
      foreach my $column (@line){
        $columns{$column}{c} = $c;
        $c++;
      }

      #Get the per lib names used
      foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
        my $prefix = $align_builds->{$name}->{prefix};

        foreach my $column (keys %columns){
          if ($column =~ /$prefix\_(\S+)\_ref\_count$/){
            my $lib = $1;
            if (defined($align_builds->{$name}->{libs})){
              my $libs = $align_builds->{$name}->{libs};
              $libs->{$lib}->{rep} = (keys %{$libs}) + 1;
              my $per_lib_ref_col = $prefix . "_" . $lib . "_ref_count";
              my $per_lib_var_col = $prefix . "_" . $lib . "_var_count";
              my $per_lib_vaf_col = $prefix . "_" . $lib . "_VAF";
              unless ($columns{$per_lib_ref_col} && $columns{$per_lib_var_col} && $columns{$per_lib_vaf_col}){
                die "Could not find expected columns: $per_lib_ref_col $per_lib_var_col $per_lib_ref_col";
              }
            }else{
              my %tmp;
              $tmp{$lib}{rep} = 1;
              $align_builds->{$name}->{libs} = \%tmp;
            }
          }
        }
      }

      #Set up the header for the per-lib count section
      foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
        my $prefix = $align_builds->{$name}->{prefix};
        my $libs = $align_builds->{$name}->{libs};
        foreach my $lib (sort {$libs->{$a}->{rep} <=> $libs->{$b}->{rep}} keys %{$libs}){
          my $rep = $libs->{$lib}->{rep};
          my $per_lib_ref_col = $prefix . "_rep" . $rep . "_ref_count";
          push(@per_lib_header, $per_lib_ref_col);
          my $per_lib_var_col = $prefix . "_rep" . $rep . "_var_count";
          push(@per_lib_header, $per_lib_var_col);
          my $per_lib_vaf_col = $prefix . "_rep" . $rep . "_VAF";
          push(@per_lib_header, $per_lib_vaf_col);
        }
      }

      $l++;
      next;
    }
    $l++;

    my ($chr, $start, $stop, $ref, $var) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
    my $v = $chr . "_$start" . "_$stop" . "_$ref" . "_$var";
    die $self->error_message("parsed a variant that is not defined in the variant hash") unless $variants->{$v};

    my @per_lib_counts;
    foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
      my $prefix = $align_builds->{$name}->{prefix};
      my $libs = $align_builds->{$name}->{libs};
      foreach my $lib (sort {$libs->{$a}->{rep} <=> $libs->{$b}->{rep}} keys %{$libs}){
        my $rep = $libs->{$lib}->{rep};
        my $per_lib_ref_col = $prefix . "_" . $lib . "_ref_count";
        my $per_lib_var_col = $prefix . "_" . $lib . "_var_count";
        my $per_lib_vaf_col = $prefix . "_" . $lib . "_VAF";

        if (defined($line[$columns{$per_lib_ref_col}{c}])){
          push(@per_lib_counts, $line[$columns{$per_lib_ref_col}{c}]);
        }else{
          push(@per_lib_counts, "NA");
        }

        if (defined($line[$columns{$per_lib_var_col}{c}])){
          push(@per_lib_counts, $line[$columns{$per_lib_var_col}{c}]);
        }else{
          push(@per_lib_counts, "NA");
        }

        if (defined($line[$columns{$per_lib_vaf_col}{c}])){
          push(@per_lib_counts, $line[$columns{$per_lib_vaf_col}{c}]);
        }else{
          push(@per_lib_counts, "NA");
        }

        #Count tumor libraries that support the variant with at least N variant read counts
        unless ($prefix =~ /normal/i){
          if (defined($line[$columns{$per_lib_var_col}{c}])){
            if ($line[$columns{$per_lib_var_col}{c}] =~ /\d+/){
              $tumor_var_supporting_libs++ if ($line[$columns{$per_lib_var_col}{c}] >= $self->min_reads_per_lib);
            }
          }
        }
      }
    }
    $variants->{$v}->{per_lib_counts} = \@per_lib_counts;
    $variants->{$v}->{tumor_var_supporting_libs} = $tumor_var_supporting_libs;
  }
  close($var_fh);
  my $per_lib_header = join("\t", @per_lib_header);

  return $per_lib_header;
}


sub apply_filter_list{
  my $self = shift;
  my %args = @_;
  my $variants = $args{'-variants'};

  my $filtered_variants = 0;
  my $var_fh = Genome::Sys->open_file_for_reading($self->variant_filter_list);
  while(<$var_fh>){
    chomp($_);
    next if ($_ =~ /chr/);
    if ($_ =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
      my ($chr, $start, $stop, $ref, $var) = ($1, $2, $3, $4, $5);
      my $v = $chr . "_$start" . "_$stop" . "_$ref" . "_$var";
      die $self->error_message("parsed a variant that is not defined in the variant hash") unless $variants->{$v};
      $variants->{$v}->{filtered} = "Filter_List";
      $filtered_variants++;
    }
  }

  $self->status_message("Filtered $filtered_variants variants because they were specified in the variant filter list: " . $self->variant_filter_list);

  close($var_fh);

  return;
}


sub apply_variant_filters{
  my $self = shift;
  my %args = @_;
  my $variants = $args{'-variants'};
  my $max_normal_vaf = $self->max_normal_vaf;
  my $min_tumor_var_count = $self->min_tumor_var_count;
  my $min_tumor_vaf = $self->min_tumor_vaf;
  my $min_coverage = $self->min_coverage;
  my $max_gmaf = $self->max_gmaf;
  my $min_tumor_var_supporting_libs = $self->min_tumor_var_supporting_libs;
  my $min_snv_caller_count = $self->min_snv_caller_count;
  my $min_indel_caller_count = $self->min_indel_caller_count;

  foreach my $v (keys %{$variants}){
    my $max_normal_vaf_observed = $variants->{$v}->{max_normal_vaf_observed};
    my $max_tumor_var_count_observed = $variants->{$v}->{max_tumor_var_count_observed};
    my $max_tumor_vaf_observed = $variants->{$v}->{max_tumor_vaf_observed};
    my $min_coverage_observed = $variants->{$v}->{min_coverage_observed};
    my $gmaf = $variants->{$v}->{gmaf};
    my $tumor_var_supporting_libs = $variants->{$v}->{tumor_var_supporting_libs} if defined($variants->{$v}->{tumor_var_supporting_libs});
    my $variant_source_caller_count = $variants->{$v}->{variant_source_caller_count};

    #Normal VAF filter
    if ($max_normal_vaf_observed =~ /\d+/){
      $variants->{$v}->{filtered} .= "Max_Normal_VAF," if ($max_normal_vaf_observed > $max_normal_vaf);
    }elsif($max_normal_vaf_observed eq "NA"){
      $variants->{$v}->{filtered} .= "Max_Normal_VAF,";
    }

    #Tumor var count filter
    if ($max_tumor_var_count_observed =~ /\d+/){
      $variants->{$v}->{filtered} .= "Min_Tumor_Var_Count," if ($max_tumor_var_count_observed < $min_tumor_var_count);
    }elsif($max_tumor_var_count_observed eq "NA"){
      $variants->{$v}->{filtered} .= "Min_Tumor_Var_Count,";
    }

    #Tumor VAF filter
    if ($max_tumor_vaf_observed =~ /\d+/){
      $variants->{$v}->{filtered} .= "Min_Tumor_VAF," if ($max_tumor_vaf_observed < $min_tumor_vaf);
    }elsif($max_tumor_vaf_observed eq "NA"){
      $variants->{$v}->{filtered} .= "Min_Tumor_VAF,";
    }

    #Coverage filter
    if ($min_coverage_observed =~ /\d+/){
      $variants->{$v}->{filtered} .= "Min_Coverage," if ($min_coverage_observed < $min_coverage);
    }elsif($min_coverage_observed eq "NA"){
      $variants->{$v}->{filtered} .= "Min_Coverage,";
    }

    #GMAF filter
    if ($gmaf =~ /\d+/){
      $variants->{$v}->{filtered} .= "Max_GMAF," if (($gmaf*100) > $max_gmaf);
    }

    #Min library support filter
    if (defined($tumor_var_supporting_libs) && $min_tumor_var_supporting_libs){
      $variants->{$v}->{filtered} .= "Min_Library_Support," if ($tumor_var_supporting_libs < $min_tumor_var_supporting_libs);
    }

    #Min variant caller count filters (variant_source_caller_count)
    if ($variants->{$v}->{variant_type} =~ /SNP/i){
      #Min SNV variant caller count filter
      $variants->{$v}->{filtered} .= "Min_SNV_Caller_Count," if ($variant_source_caller_count < $min_snv_caller_count);
    }elsif($variants->{$v}->{variant_type} =~ /INS|DEL/i){
      #Min INDEL variant caller count filter
      $variants->{$v}->{filtered} .= "Min_INDEL_Caller_Count," if ($variant_source_caller_count < $min_indel_caller_count);
    }else{
      die $self->error_message("unrecognized variant type: $variants->{$v}->{variant_type}");
    }

    #If not filtered by any category, don't filter. Keep this at the end.
    if($variants->{$v}->{filtered} eq "") {
      $variants->{$v}->{filtered} = 0;
    }
  }

  return;
}

sub get_result_files {
  my $self = shift;
  my $case_name = shift;
  my $result_files;

  #Write out final tsv files (filtered and unfiltered), a clean version with useless columns removed, and an Excel spreadsheet version of the final file
  my $final_unfiltered_tsv = $self->outdir . "$case_name" . "_final_unfiltered.tsv"; #OUT1
  my $final_unfiltered_clean_tsv = $self->outdir . "$case_name" . "_final_unfiltered_clean.tsv"; #OUT1b
  my $final_filtered_tsv = $self->outdir . "$case_name" . "_final_filtered.tsv"; #OUT2
  my $final_filtered_clean_tsv = $self->outdir . "$case_name" . "_final_filtered_clean.tsv"; #OUT3
  my $final_filtered_coding_clean_tsv = $self->outdir . "$case_name" . "_final_filtered_coding_clean.tsv"; #OUT4
  my $final_filtered_clean_xls = $self->outdir . "$case_name" . "_final_filtered_clean.xls"; #OUT5
  my $final_filtered_coding_clean_xls = $self->outdir . "$case_name" . "_final_filtered_coding_clean.xls"; #OUT6

  #Store the result files paths and pass out to be used in the visualization step
  $result_files->{final_unfiltered_tsv} = $final_unfiltered_tsv;
  $result_files->{final_unfiltered_clean_tsv} = $final_unfiltered_clean_tsv;
  $result_files->{final_filtered_tsv} = $final_filtered_tsv;
  $result_files->{final_filtered_clean_tsv} = $final_filtered_clean_tsv;
  $result_files->{final_filtered_coding_clean_tsv} = $final_filtered_coding_clean_tsv;
  $result_files->{final_filtered_clean_xls} = $final_filtered_clean_xls;
  $result_files->{final_filtered_coding_clean_xls} = $final_filtered_coding_clean_xls;

  return $result_files;
}

sub print_final_files{
  my $self = shift;
  my %args = @_;
  my $variants = $args{'-variants'};
  my $grand_anno_count_file = $args{'-grand_anno_count_file'};
  my $case_name = $args{'-case_name'};
  my $align_builds = $args{'-align_builds'};
  my $per_lib_header = $args{'-per_lib_header'};
  my $trv_type_filter = $self->trv_type_filter;
  my $result_files = $self->get_result_files($case_name);

  my $anno_fh = Genome::Sys->open_file_for_reading($grand_anno_count_file);
  my $final_unfiltered_fh = Genome::Sys->open_file_for_writing($result_files->{final_unfiltered_tsv});
  my $final_unfiltered_clean_fh = Genome::Sys->open_file_for_writing($result_files->{final_unfiltered_clean_tsv});
  my $final_filtered_fh = Genome::Sys->open_file_for_writing($result_files->{final_filtered_tsv});
  my $final_filtered_clean_fh = Genome::Sys->open_file_for_writing($result_files->{final_filtered_clean_tsv});
  my $final_filtered_coding_clean_fh = Genome::Sys->open_file_for_writing($result_files->{final_filtered_coding_clean_tsv});

  my @skip = qw (gene_name transcript_species transcript_source transcript_version transcript_status c_position ucsc_cons domain all_domains deletion_substructures transcript_error gene_name_source);
  my %skip_columns;
  foreach my $name (@skip){
    $skip_columns{$name}=1;
  }
  my @include_col_pos;

  my $header = 1;
  my $l = 0;
  my %columns;
  while(<$anno_fh>){
    chomp($_);

    #Remove ugly GMAF=$val from annotation lines to tidy up output
    $_ =~ s/GMAF\=//g;

    my @line = split("\t", $_);

    #BAM readcount can give empty cells for sites not counted (e.g. indels that are larger than the max size allowed).  Replace empty cells with 'NA'
    my @tmp_line;
    foreach my $val (@line){
      $val = "NA" unless (defined($val));
      push(@tmp_line, $val);
    }
    @line = @tmp_line;

    if ($header){
      $header = 0;
      my $c = 0;
      foreach my $col (@line){
        $columns{$col}{c} = $c;
        $c++;
      }
      foreach my $col_name (sort {$columns{$a}{c} <=> $columns{$b}{c}} keys %columns){
        unless($skip_columns{$col_name}){
          push(@include_col_pos, $columns{$col_name}{c});
        }
      }

      #Print headers for each out file
      my $header_extension = "min_coverage_observed\tmax_normal_vaf_observed\tmax_tumor_vaf_observed\tvariant_source_callers\tvariant_source_caller_count\tdata_type\tfiltered";
      my $full_header = "$_"."\t$header_extension";
      $full_header .= "\t$per_lib_header" if $per_lib_header;
      print $final_unfiltered_fh "$full_header\n";
      print $final_filtered_fh "$full_header\n";

      my @include_values = @line[@include_col_pos];
      my $include_values_string = join("\t", @include_values);
      my $short_header = "$include_values_string"."\t$header_extension";
      $short_header .= "\t$per_lib_header" if $per_lib_header;
      print $final_unfiltered_clean_fh "$short_header\n";
      print $final_filtered_clean_fh "$short_header\n";
      print $final_filtered_coding_clean_fh "$short_header\n";
      next;
    }

    my ($chr, $start, $stop, $ref, $var) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
    my $v = $chr . "_$start" . "_$stop" . "_$ref" . "_$var";
    die $self->error_message("parsed a variant that is not defined in the variant hash") unless $variants->{$v};

    my @per_lib_counts = @{$variants->{$v}->{per_lib_counts}} if defined($variants->{$v}->{per_lib_counts});
    my $per_lib_count_line = join("\t", @per_lib_counts) if defined($variants->{$v}->{per_lib_counts});
    $variants->{$v}->{data_type} =~ s/,$//;
    $variants->{$v}->{filtered} =~ s/,$//;

    my $line_extension = "$variants->{$v}->{min_coverage_observed}\t$variants->{$v}->{max_normal_vaf_observed}\t$variants->{$v}->{max_tumor_vaf_observed}\t$variants->{$v}->{variant_source_callers}\t$variants->{$v}->{variant_source_caller_count}\t$variants->{$v}->{data_type}\t$variants->{$v}->{filtered}";
    my $full_line = "$_\t$line_extension";
    $full_line .= "\t$per_lib_count_line" if defined($per_lib_count_line);

    print $final_unfiltered_fh "$full_line\n";
    unless ($variants->{$v}->{filtered}){
      print $final_filtered_fh "$full_line\n";
    }

    my @include_values = @line[@include_col_pos];
    my $include_values_string = join("\t", @include_values);

    my $short_line = "$include_values_string"."\t$line_extension";
    $short_line .= "\t$per_lib_count_line" if defined($per_lib_count_line);

    print $final_unfiltered_clean_fh "$short_line\n";
    unless ($variants->{$v}->{filtered}){
      print $final_filtered_clean_fh "$short_line\n";
    }

    #apply a transcript variant type filter to define the 'coding' result
    my $trv_type = $variants->{$v}->{trv_type};
    unless ($variants->{$v}->{filtered}){
      unless ($trv_type_filter =~ /$trv_type/i){
        print $final_filtered_coding_clean_fh "$short_line\n";
      }
    }
    $l++;
  }
  close($anno_fh);
  close($final_unfiltered_fh);
  close($final_unfiltered_clean_fh);
  close($final_filtered_fh);
  close($final_filtered_clean_fh);
  close($final_filtered_coding_clean_fh);

  # convert master table to excel
  my $final_filtered_clean_workbook  = Spreadsheet::WriteExcel->new("$result_files->{final_filtered_clean_xls}");
  my $final_filtered_clean_worksheet = $final_filtered_clean_workbook->add_worksheet();
  my $in_fh = Genome::Sys->open_file_for_reading($result_files->{final_filtered_clean_tsv});
  my $row=0;
  while(<$in_fh>){
    chomp($_);
    if ($row == 0){
      $_ =~ s/\_/ /g;
    }
    my @F = split("\t",$_);
    for(my $i=0;$i<@F;$i++){
      $final_filtered_clean_worksheet->write($row, $i, $F[$i]);
    }
    $row++;
  }
  close($in_fh);
  $final_filtered_clean_workbook->close();

  my $final_filtered_coding_clean_workbook  = Spreadsheet::WriteExcel->new("$result_files->{final_filtered_coding_clean_xls}");
  my $final_filtered_coding_clean_worksheet = $final_filtered_coding_clean_workbook->add_worksheet();
  $in_fh = Genome::Sys->open_file_for_reading($result_files->{final_filtered_coding_clean_tsv});
  $row=0;
  while(<$in_fh>){
    chomp($_);
    if ($row == 0){
      $_ =~ s/\_/ /g;
    }
    my @F = split("\t",$_);
    for(my $i=0;$i<@F;$i++){
      $final_filtered_coding_clean_worksheet->write($row, $i, $F[$i]);
    }
    $row++;
  }
  close($in_fh);
  $final_filtered_coding_clean_workbook->close();

  return($result_files);
}

sub create_plots{
  my $self = shift;
  my %args = @_;
  my $result_files = $args{'-result_files'};
  my $align_builds = $args{'-align_builds'};
  my $case_name = $args{'-case_name'};

  my $final_filtered_clean_tsv = $result_files->{final_filtered_clean_tsv};
  my $final_filtered_coding_clean_tsv = $result_files->{final_filtered_coding_clean_tsv};

  #Get the header for the file to be fed into R to determine per-lib VAF column positions
  my $tmp_fh = Genome::Sys->open_file_for_reading($final_filtered_clean_tsv);
  my $header = <$tmp_fh>;
  close($tmp_fh);
  chomp($header);
  my @cols = split("\t", $header);
  my $p = 0;
  my %columns;
  foreach my $col (@cols){
    $columns{$col}{p} = $p;
    $p++;
  }

  #Set the R script that will process output from this perl script
  #TODO: Right now, this script will only work for a very particular situation (normal, dayX_tumor, dayY_tumor) - Make it more flexible ...
  my $build_count = keys %{$align_builds};
  my @prefixes;
  my @combined_vaf_cols;
  my %rep_vaf_cols;
  my @timepoint_names;
  my @timepoint_positions;
  my @sample_types;
  #Resolve the names of the samples, and columns containing per_library read counts
  foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys %{$align_builds}){
    my $prefix = $align_builds->{$name}->{prefix};
    push(@prefixes, $prefix);
    my $day = $align_builds->{$name}->{day};
    push(@timepoint_names, $day);
    my $timepoint_position = $align_builds->{$name}->{timepoint_position};
    push(@timepoint_positions, $timepoint_position);
    my $sample_type = $align_builds->{$name}->{sample_common_name};
    push(@sample_types, $sample_type);
    my $combined_vaf_col_name = $prefix . "_VAF";
    push(@combined_vaf_cols, $columns{$combined_vaf_col_name}{p} + 1);

    my $libs = $align_builds->{$name}->{libs};
    if ($self->per_library){
      #If replicate library analysis is being run define the replicate VAF columns
      foreach my $lib (sort {$libs->{$a}->{rep} <=> $libs->{$b}->{rep}} keys %{$libs}){
        my $rep = $libs->{$lib}->{rep};
        my $per_lib_vaf_col = $prefix . "_rep" . $rep . "_VAF";

        die $self->error_message("could not resolve column $per_lib_vaf_col in file $final_filtered_clean_tsv") unless ($columns{$per_lib_vaf_col});

        if (defined($rep_vaf_cols{$prefix})){
          my $cols = $rep_vaf_cols{$prefix}{cols};
          push (@{$cols}, $columns{$per_lib_vaf_col}{p} + 1);
        }else{
          my @cols;
          push (@cols, $columns{$per_lib_vaf_col}{p} + 1);
          $rep_vaf_cols{$prefix}{cols} = \@cols;
        }
      }
    }else{
      #If replicate library analysis is NOT being run simply use the combined VAF values
      my @cols;
      push (@cols, $columns{$combined_vaf_col_name}{p} + 1);
      $rep_vaf_cols{$prefix}{cols} = \@cols;
    }
  }

  my $r_script = __FILE__ . '.R';

  #Run the R script to generate some per-library readcount and VAF plots from the complete clean tsv file
  my @vaf_cols;
  foreach my $prefix (@prefixes){
    my @vaf_cols_sample = @{$rep_vaf_cols{$prefix}{cols}};
    my $vaf_cols_sample_string = join(" ", @vaf_cols_sample);
    push (@vaf_cols, "\"$vaf_cols_sample_string\"");
  }
  my $vaf_cols_string = join(" ", @vaf_cols);
  my $sample_types_string = join(" ", @sample_types);
  my $timepoint_names_string = join(" ", @timepoint_names);
  my $timepoint_positions_string = join(" ", @timepoint_positions);
  my $target_gene_list_name = $self->target_gene_list_name;

  if ($self->per_library){
    my $prefix_string = join(" ", @prefixes);
    my $combined_vaf_col_string = join(" ", @combined_vaf_cols);

    my $outdir1 = $self->outdir . "filtered_pdfs/";
    mkdir($outdir1);
    my $r_cmd1 = "$r_script $case_name $final_filtered_clean_tsv \"$prefix_string\" \"$combined_vaf_col_string\" \"$target_gene_list_name\" $outdir1 \"$sample_types_string\" \"$timepoint_names_string\" \"$timepoint_positions_string\" $vaf_cols_string";
    Genome::Sys->shellcmd(cmd => $r_cmd1);

    my $outdir2 = $self->outdir . "filtered_coding_pdfs/";
    mkdir($outdir2);
    my $r_cmd2 = "$r_script $case_name $final_filtered_coding_clean_tsv \"$prefix_string\" \"$combined_vaf_col_string\" \"$target_gene_list_name\" $outdir2 \"$sample_types_string\" \"$timepoint_names_string\" \"$timepoint_positions_string\" $vaf_cols_string";
    Genome::Sys->shellcmd(cmd => $r_cmd2);
  }
  return;
}



1;

