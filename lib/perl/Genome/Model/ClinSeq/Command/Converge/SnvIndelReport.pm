package Genome::Model::ClinSeq::Command::Converge::SnvIndelReport;
use strict;
use warnings;
use Genome;
use Data::Dumper;
use List::MoreUtils qw/ uniq /;
use Genome::Info::IUB;

class Genome::Model::ClinSeq::Command::Converge::SnvIndelReport {
    is => 'Genome::Model::ClinSeq::Command::Converge::Base',
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
        subject_labels_file => {
                is => 'FilesystemPath',
                doc => 'Use a custom subjects_legend file.  First run without this option, then copy the legend created, modify and then specify with this option. (use to change order, sample names, etc.)',
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
              default => 0,
              doc => 'In per library read counting, the minimum variant supporting reads for a library to be said to support a variant',
        },
        min_tumor_var_supporting_libs => {
              is => 'Number',
              default => 0,
              doc => 'In per library analysis, variants with less than this number of tumor libraries supporting will be filtered out',
        },
        test => {
              is => 'Number',
              doc => 'Only import this many variants (for testing purposes)',
        }

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

  return @errors;
}


sub execute {
  my $self = shift;
  my @builds = $self->builds;
  my $outdir = $self->outdir;
  $outdir .= "/" unless ($outdir =~ /\/$/);

  #Get human readable names hash, keyed on build id
  my $subject_labels = $self->resolve_clinseq_subject_labels;

  #Get the case name for this set of builds, if more than one are found, warn the user
  my $case_name = $self->get_case_name;
  $self->status_message("Producing report for individual: $case_name");

  #Print out a table of subject names for reference
  $self->print_subject_table('-subject_labels'=>$subject_labels);

  #Gather variants for the tiers specified by the user from each build. Note which build each came from.
  #Get these from the underlying somatic-variation builds.
  #Annotate all variants (gmt annotate transcript-variants --help)
  my $somatic_builds = $self->resolve_somatic_builds;
  my $annotation_build_name = $self->annotation_build->name;
  $self->status_message("Using annotation build: $annotation_build_name");
  my $result = $self->gather_variants('-somatic_builds'=>$somatic_builds);
  my $variants = $result->{'variants'};
  my $header = $result->{'header'};
  
  #TODO: Perform hotspot analysis of critical mutation sites
  #If these variants are not already in the list, inject them so that all subsequent steps are performed on these sites...
  #As a starting point create a database of these variants, separate into SNVs and Indels and then tier them

  #Make note of which variants correspond to $self->target_gene_list
  $self->intersect_target_gene_list('-variants'=>$variants);

  #Pull in source variant caller info from the clin-seq build results
  $self->get_variant_caller_sources('-variants'=>$variants);

  #Print a grand table of all annotated variants with extra annotation columns (tier, target gene list)
  my $grand_anno_file = $self->print_grand_anno_table('-variants'=>$variants, '-header'=>$header);

  #Identify underlying reference-alignments, sample names, sample common names, and timepoints (if available)
  my $align_builds = $self->get_ref_align_builds('-somatic_builds'=>$somatic_builds);

  #Get bam-readcounts for all positions for all BAM files
  my $grand_anno_count_file = $self->add_read_counts('-align_builds'=>$align_builds, '-grand_anno_file'=>$grand_anno_file);
  
  #TODO: This will need to be done at the level of replicate libraries once that readcounting tool is availble.  For now just use standard bam read counts
  #Example usage for per-library BAM readcounting from Dave Larson
  my $grand_anno_per_lib_count_file = $self->add_per_lib_read_counts('-align_builds'=>$align_builds, '-grand_anno_file'=>$grand_anno_count_file);

  #Parse the BAM read count info and gather minimal data needed to apply filters:
  #Max normal VAF, Min Coverage
  $self->parse_read_counts('-align_builds'=>$align_builds, '-grand_anno_count_file'=>$grand_anno_count_file, '-variants'=>$variants);

  #print Dumper $variants;
  #Parse the per lib BAM read count info
  my $per_lib_header = $self->parse_per_lib_read_counts('-align_builds'=>$align_builds, '-grand_anno_per_lib_count_file'=>$grand_anno_per_lib_count_file, '-variants'=>$variants);

  #Apply some automatic filters:
  #Normal VAF > $self->max_normal_vaf
  #Min coverage > $self->min_coverage (position must be covered at this level across all samples)
  #GMAF > $self->max_gmaf
  $self->apply_variant_filters('-variants'=>$variants);

  #TODO: Make note of which variants lie within a particular set of regions of interest (e.g. nimblegen v3 + AML RMG)
  #TODO: Have option to filter these out


  #TODO: Filter out variants falling within false positive regions/genes?


  #TODO: Write out final tsv files (filtered and unfiltered), a clean version with useless columns removed, and an Excel spreadsheet version of the final file
  $self->print_final_files('-variants'=>$variants, '-grand_anno_count_file'=>$grand_anno_count_file, '-case_name'=>$case_name, '-align_builds'=>$align_builds, '-per_lib_header'=>$per_lib_header);

  #Produce some visualizations for variants in the target gene list as well as all high quality variants, the higher VAF variants etc.
  #If there are two tumors, Produce SciClone plots showing tumor1 vs. tumor2 VAF and with gene names labelled


  print "\n\n";

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
  my $subject_labels = $args{'-subject_labels'};

  my $outfile = $self->outdir . "/subjects_legend.txt";
  open (OUT, ">$outfile") || die $self->error_message("Could not open output file: $outfile for writing");
  print OUT "build_id\tname\tname_abr\torder\n";
  foreach my $bid (sort {$subject_labels->{$a}->{order} <=> $subject_labels->{$b}->{order}} keys %{$subject_labels}){
    my $name = $subject_labels->{$bid}->{name};
    my $name_abr = $subject_labels->{$bid}->{name_abr};
    my $order = $subject_labels->{$bid}->{order};
    print OUT "$bid\t$name\t$name_abr\t$order\n";
  }

  close(OUT);

  return;
}

sub get_case_name{
  my $self = shift;
  my @builds = $self->builds;
  
  #First attempt to find a common name
  my %common_names;
  my $final_common_name;
  foreach my $build (@builds){
    my $name = $self->get_final_common_name('-clinseq_build'=>$build);
    $common_names{$name}=1 if $name;
    $final_common_name = $name;
  }
  my $ncn = keys %common_names;
  if ($ncn > 1){
    die $self->error_message("$ncn cases found among these builds, this tool is meant to operate on builds from a single individual"); 
  }

  #Second attempt to find a name
  my %names;
  my $final_name;
  foreach my $build (@builds){
    my $name = $self->get_final_name('-clinseq_build'=>$build);
    $names{$name}=1 if $name;
    $final_name = $name;
  }
  my $nn = keys %names;
  if ($nn > 1){
    die $self->error_message("$nn cases found among these builds, this tool is meant to operate on builds from a single individual"); 
  }

  my $resolved_name;
  if ($final_common_name){
    $resolved_name = $final_common_name;
  }elsif($final_name){
    $resolved_name = $final_name
  }else{
    die $self->error_message("could not find an individual common_name or name in these builds");
  }

  return $resolved_name;
}


sub gather_variants{
  my $self = shift;
  my %args = @_;
  my $somatic_builds = $args{'-somatic_builds'};

  my $tiers_list = $self->tiers;
  my @tiers = split(",", $tiers_list);

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

  my $bed_dir = $self->outdir . "/bed_files/";
  mkdir ($bed_dir) unless (-e $bed_dir && -d $bed_dir);

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
      open(VAR, $file) || die $self->error_message("Could not open var file: $file");
      open(NEW, ">$new_file") || die $self->error_message("Could not open new var file: $new_file");
      while(<VAR>){
        $c++;
        if ($self->test){last if $c > $self->test;}
        chomp($_);
        my @line = split("\t", $_);
        my ($chr, $start, $end, $var)  = ($line[0], $line[1], $line[2], $line[3]);
        #Make indel format consistent
        $var =~ s/\*/0/g;
        print NEW "$chr\t$start\t$end\t$var\n";
      }
      close(VAR);
      close(NEW);
    }
    $bed_files{$file}{clean_bed_file} = $new_file;

    #Run the annnotator on each bed file and add RSIDs and GMAFs for all variants
    #gmt annotate transcript-variants --variant-bed-file='' --output-file='' --annotation-filter=top --reference-transcripts='NCBI-human.ensembl/74_37'
    my $anno_file = $bed_files{$file}{clean_bed_file} . ".anno";
    if (-e $anno_file){
      $self->warning_message("Using pre-generated file: $anno_file");
    }else{
      my $annotate_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
            variant_bed_file=>$bed_files{$file}{clean_bed_file}, 
            output_file=>$anno_file,
            annotation_filter=>'top',
            reference_transcripts=>$self->annotation_build->name,
          );
      my $r1 = $annotate_cmd->execute();
      die $self->error_message("annotation cmd unsuccessful") unless ($r1);
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
    open (VAR, $rsid_file) || die $self->error_message("could not open file: $rsid_file");
    while(<VAR>){
      chomp($_);
      if ($_ =~ /^chromosome\_name/){
        $header = $_;
        $headers{$header}=1;
        next;
      }
      my @line = split("\t", $_);
      my ($chr, $start, $stop, $ref, $var) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
      my $ensembl_gene_id = $line[23];
      my $trv_type = $line[13];
      my $v = $chr . "_$start" . "_$stop" . "_$ref" . "_$var";
      $variants{$v}{anno_line} = $_;
      $variants{$v}{tier} = $tier;
      $variants{$v}{ensembl_gene_id} = $ensembl_gene_id;
      $variants{$v}{trv_type} = $trv_type;
    }
    close(VAR);
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
  my $grand_anno_file = $self->outdir . "/variants.all.anno";
  if (-e $grand_anno_file){
    $self->warning_message("using pre-generated file: $grand_anno_file");
  }else{
    open (ANNO, ">$grand_anno_file") || die $self->error_message("could not open grand anno file: $grand_anno_file");
    my $target_gene_list_name = $self->target_gene_list_name;
    my $new_header = "$header\ttier\t$target_gene_list_name";
    print ANNO "$new_header\n";
    foreach my $v (sort keys %{$variants}){
      print ANNO "$variants->{$v}->{anno_line}\t$variants->{$v}->{tier}\t$variants->{$v}->{target_gene_list_match}\n";
    }
    close(ANNO);
  }

  #delete the anno line now that it is not needed
  foreach my $v (keys %{$variants}){
    delete $variants->{$v}->{anno_line};
  }

  return($grand_anno_file);
}


sub get_ref_align_builds{
  my $self = shift;
  my %args = @_;
  my $somatic_builds = $args{'-somatic_builds'};

  my %ref_builds;

  my $sort_on_time_point = 0;

  foreach my $somatic_build_id (keys %{$somatic_builds}){
    my $build_type = $somatic_builds->{$somatic_build_id}->{type};
    my $somatic_build = $somatic_builds->{$somatic_build_id}->{build};
    my $normal_build = $somatic_build->normal_build;
    my $normal_subject_name = $normal_build->subject->name;
    my $tumor_build = $somatic_build->tumor_build;
    my $tumor_subject_name = $tumor_build->subject->name;
    my $normal_refalign_name = $normal_subject_name . "_$build_type"."_normal";
    my $tumor_refalign_name = $tumor_subject_name . "_$build_type"."_tumor";
    my $normal_bam_path = $normal_build->whole_rmdup_bam_file;
    my $tumor_bam_path = $tumor_build->whole_rmdup_bam_file;
    my @normal_timepoints = $normal_build->subject->attributes(attribute_label => "timepoint", nomenclature => "caTissue");
    my @tumor_timepoints = $tumor_build->subject->attributes(attribute_label => "timepoint", nomenclature => "caTissue");

    my $normal_time_point = "day0";
    if (@normal_timepoints){
      $normal_time_point = $normal_timepoints[0]->attribute_value;
      $normal_time_point =~ s/\s+//g;
      $sort_on_time_point = 1;
    }
    $normal_refalign_name .= "_$normal_time_point";

    my $tumor_time_point = "day0";
    if (@tumor_timepoints){
      $tumor_time_point = $tumor_timepoints[0]->attribute_value;
      $tumor_time_point =~ s/\s+//g;
      $sort_on_time_point = 1;
    }
    $tumor_refalign_name .= "_$tumor_time_point";

    $ref_builds{$normal_refalign_name}{type} = $build_type;
    $ref_builds{$normal_refalign_name}{sample_name} = $normal_subject_name;
    $ref_builds{$normal_refalign_name}{bam_path} = $normal_bam_path;
    $ref_builds{$normal_refalign_name}{time_point} = "normal_" . $normal_time_point;

    $ref_builds{$tumor_refalign_name}{type} = $build_type;
    $ref_builds{$tumor_refalign_name}{sample_name} = $tumor_subject_name;
    $ref_builds{$tumor_refalign_name}{bam_path} = $tumor_bam_path;
    $ref_builds{$tumor_refalign_name}{time_point} = "tumor_" . $tumor_time_point;
  }

  #Set an order on refalign builds (use time points if available, otherwise name)
  my $o = 0;
  if ($sort_on_time_point){
    foreach my $name (sort {$ref_builds{$a}->{time_point} cmp $ref_builds{$b}->{time_point}} keys %ref_builds){
      $o++;
      $ref_builds{$name}{order} = $o;
    }
  }else{
    foreach my $name (sort keys %ref_builds){
      $o++;
      $ref_builds{$name}{order} = $o;
    }
  }

  return(\%ref_builds);
}


sub intersect_target_gene_list{
  my $self = shift;
  my %args = @_;
  my $variants = $args{'-variants'};

  my $target_gene_list = $self->target_gene_list;
  my $target_gene_list_name = $self->target_gene_list_name;

  my %genes;
  open (GENES, $target_gene_list) || die $self->error_message("could not open target gene list: $target_gene_list");
  while(<GENES>){
    chomp($_);
    my @line = split("\t", $_);
    my $ensg = $line[0];
    $genes{$ensg}{name} = $line[1];
  }
  close(GENES);
 
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
  foreach my $build (@builds){
    my $common_name = $build->common_name;
    my $data_dir = $build->data_directory;
    my $base_path = $data_dir . "/" . $common_name . "/variant_source_callers";
    my $exome_base_path = $base_path . "/exome/";
    my $wgs_base_path = $base_path . "/wgs/";
    my @files_list = ($exome_base_path."snv_sources.tsv", $wgs_base_path."snv_sources.tsv", $exome_base_path."indel_sources.tsv", $wgs_base_path."indel_sources.tsv");
    push(@files, @files_list);
  }

  foreach my $file (@files){
    $self->status_message("processing variant caller sources file: $file");
    next unless (-e $file);
    open (SOURCES, $file) || die $self->error_message("could not open variant caller sources file: $file");
    while(<SOURCES>){
      chomp($_);
      next if ($_ =~ /^coord/);
      my @line = split("\t", $_);
      my ($chr, $start, $end, $variant, $caller_string) = ($line[1], $line[2], $line[3], $line[4], $line[7]);
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
            die $self->error_message("could not match variant to hash: $v");
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
          die $self->error_message("could not match variant to hash: $v");
        }
      }else{
        die $self->error_message("file not recognized as snv or indel: $file");
      }
    }
    close (SOURCES);
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


sub add_read_counts{
  my $self = shift;
  my %args = @_;
  my $align_builds = $args{'-align_builds'};
  my $grand_anno_file = $args{'-grand_anno_file'};

  my @bam_files;
  my @time_points;
  my @samples;
  my @names;
  foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
    push(@bam_files, $align_builds->{$name}->{bam_path});
    push(@time_points, $align_builds->{$name}->{time_point});
    push(@samples, $align_builds->{$name}->{sample_name});
    push(@names, $name);
  }
  my $bam_list = join(",", @bam_files);

  #Get the reference fasta
  my $reference_build = $self->resolve_clinseq_reference_build;
  my $reference_fasta = $reference_build->full_consensus_path('fa');

  #Determine header prefixes to use. In order of preference if all are unique: (time_points, samples, names)
  my @prefixes;
  my @unique_time_points = uniq @time_points;
  my @unique_samples = uniq @samples;
  my @unique_names = uniq @names;
  if (scalar(@unique_time_points) == scalar(@time_points)){
    @prefixes = @time_points;
  }elsif(scalar(@unique_samples) == scalar(@samples)){
    @prefixes = @samples;
  }elsif(scalar(@unique_names) == scalar(@names)){
    @prefixes = @names;
  }else{
    die $self->error_message("could not resolve unique prefixes for add-readcounts");
  }
  my $header_prefixes = join(",", @prefixes);

  #Record the header prefix chosen on the align_builds object
  foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
    my $prefix = shift @prefixes;
    $align_builds->{$name}->{prefix} = $prefix;
  }

  #gmt analysis coverage add-readcounts --bam-files=? --genome-build=? --output-file=? --variant-file=? [--header-prefixes=?] 
  my $output_file = $self->outdir . "/variants.all.anno.readcounts";
  if (-e $output_file){
    $self->warning_message("using pre-generated bam read count file: $output_file");
  }else{
    my $add_count_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
            bam_files=>$bam_list,
            genome_build=>$reference_fasta,
            output_file=>$output_file,
            variant_file=>$grand_anno_file,
            header_prefixes=>$header_prefixes,
          );
    my $r = $add_count_cmd->execute();
    die $self->error_message("add-readcounts cmd unsuccessful") unless ($r);
  }

  return ($output_file);
}


sub add_per_lib_read_counts{
  my $self = shift;
  my %args = @_;
  my $align_builds = $args{'-align_builds'};
  my $grand_anno_file = $args{'-grand_anno_file'};

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
    #perl -I /gscuser/dlarson/src/genome/lib/perl/ `which gmt` analysis coverage add-readcounts --bam-file /gscuser/dlarson/src/bam-readcount/test-data/test.bam --genome-build ~dlarson/src/somatic-snv-test-data/ref.fa --variant-file variants_wheader.csv --output-file add_readcount_new5.csv --per-library 
    
    my $per_lib_count_cmd = "`which gmt` analysis coverage add-readcounts --bam-files $bam_list --genome-build $reference_fasta --variant-file $grand_anno_file --output-file $output_file --header-prefixes $header_prefixes --per-library";
    $self->status_message("attempting per library read counting: \n$per_lib_count_cmd");
    Genome::Sys->shellcmd(cmd => $per_lib_count_cmd);

    #TODO: replace with proper method call once this gets pushed to master
    #my $add_count_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
    #        bam_files=>$bam_list,
    #        genome_build=>$reference_fasta,
    #        output_file=>$output_file,
    #        variant_file=>$grand_anno_file,
    #        header_prefixes=>$header_prefixes,
    #      );
    #my $r = $add_count_cmd->execute();
    #die $self->error_message("add-readcounts cmd unsuccessful") unless ($r);
  }

  return ($output_file);
}


sub parse_read_counts{
  my $self = shift;
  my %args = @_;
  my $align_builds = $args{'-align_builds'};
  my $grand_anno_count_file = $args{'-grand_anno_count_file'};
  my $variants = $args{'-variants'};

  my @prefixes;
  foreach my $name (sort {$align_builds->{$a}->{order} <=> $align_builds->{$b}->{order}} keys  %{$align_builds}){
    my $prefix = $align_builds->{$name}->{prefix};
    push(@prefixes, $prefix);
  }

  my %columns;
  my $l = 0;
  open (VAR, $grand_anno_count_file) || die $self->error_message("could not open var anno count file: $grand_anno_count_file");
  while(<VAR>){
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
      foreach my $prefix (@prefixes){
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
    my $max_normal_vaf_observed = 0;
    my $max_tumor_vaf_observed = 0;
    my $min_coverage_observed = 10000000000000000;
    my $na_found = 0;
    my @covs;
    foreach my $prefix (@prefixes){
      my $ref_count_colname = $prefix . "_ref_count";
      my $var_count_colname = $prefix . "_var_count";
      my $vaf_colname = $prefix . "_VAF";
      if ($line[$columns{$vaf_colname}{c}] eq "NA"){
        $na_found = 1;
        push(@covs, "NA");
        next;
      }
      if ($vaf_colname =~ /normal/){
        my $normal_vaf = $line[$columns{$vaf_colname}{c}];
        $max_normal_vaf_observed = $normal_vaf if ($normal_vaf > $max_normal_vaf_observed);
      }else{
        my $tumor_vaf = $line[$columns{$vaf_colname}{c}];
        $max_tumor_vaf_observed = $tumor_vaf if ($tumor_vaf > $max_tumor_vaf_observed);
      }
      my $coverage = $line[$columns{$ref_count_colname}{c}] + $line[$columns{$var_count_colname}{c}];
      push(@covs, $coverage);
      $min_coverage_observed = $coverage if ($coverage < $min_coverage_observed);
    }
    $max_normal_vaf_observed = "NA" if $na_found;
    $min_coverage_observed = "NA" if $na_found;

    $variants->{$v}->{max_normal_vaf_observed} = $max_normal_vaf_observed;
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
  close(VAR);

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
  open (VAR, $grand_anno_per_lib_count_file) || die $self->error_message("could not open var anno count file: $grand_anno_per_lib_count_file");
  while(<VAR>){
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
  my $per_lib_header = join("\t", @per_lib_header);

  return $per_lib_header;
}


sub apply_variant_filters{
  my $self = shift;
  my %args = @_;
  my $variants = $args{'-variants'};
  my $max_normal_vaf = $self->max_normal_vaf;
  my $min_tumor_vaf = $self->min_tumor_vaf;
  my $min_coverage = $self->min_coverage;
  my $max_gmaf = $self->max_gmaf;
  my $min_tumor_var_supporting_libs = $self->min_tumor_var_supporting_libs;
  
  foreach my $v (keys %{$variants}){
    $variants->{$v}->{filtered} = 0;
    my $max_normal_vaf_observed = $variants->{$v}->{max_normal_vaf_observed};
    my $max_tumor_vaf_observed = $variants->{$v}->{max_tumor_vaf_observed};
    my $min_coverage_observed = $variants->{$v}->{min_coverage_observed};
    my $gmaf = $variants->{$v}->{gmaf};
    my $tumor_var_supporting_libs = $variants->{$v}->{tumor_var_supporting_libs} if defined($variants->{$v}->{tumor_var_supporting_libs});
  
    #Normal VAF filter
    if ($max_normal_vaf_observed =~ /\d+/){
      $variants->{$v}->{filtered} = 1 if ($max_normal_vaf_observed > $max_normal_vaf);
    }elsif($max_normal_vaf_observed eq "NA"){
      $variants->{$v}->{filtered} = 1;
    }

    #Tumor VAF filter
    if ($max_tumor_vaf_observed =~ /\d+/){
      $variants->{$v}->{filtered} = 1 if ($max_tumor_vaf_observed < $min_tumor_vaf);
    }elsif($max_tumor_vaf_observed eq "NA"){
      $variants->{$v}->{filtered} = 1;
    }

    #Coverage filter
    if ($min_coverage_observed =~ /\d+/){
      $variants->{$v}->{filtered} = 1 if ($min_coverage_observed < $min_coverage);
    }elsif($min_coverage_observed eq "NA"){
      $variants->{$v}->{filtered} = 1;
    }
    
    #GMAF filter
    if ($gmaf =~ /\d+/){
      $variants->{$v}->{filtered} = 1 if (($gmaf*100) > $max_gmaf);
    }

    #Min library support filter
    if (defined($tumor_var_supporting_libs) && $min_tumor_var_supporting_libs){
      $variants->{$v}->{filtered} = 1 if ($tumor_var_supporting_libs < $min_tumor_var_supporting_libs);
    }
  }

  return;
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


  #Write out final tsv files (filtered and unfiltered), a clean version with useless columns removed, and an Excel spreadsheet version of the final file
  my $final_unfiltered_tsv = $self->outdir . "/$case_name" . "_final_unfiltered.tsv"; #OUT1
  my $final_filtered_tsv = $self->outdir . "/$case_name" . "_final_filtered.tsv"; #OUT2
  my $final_filtered_clean_tsv = $self->outdir . "/$case_name" . "_final_filtered_clean.tsv"; #OUT3
  my $final_filtered_coding_clean_tsv = $self->outdir . "/$case_name" . "_final_filtered_coding_clean.tsv"; #OUT4
  my $final_filtered_clean_xls = $self->outdir . "/$case_name" . "_final_filtered_clean.xls"; #OUT5
  my $final_filtered_coding_clean_xls = $self->outdir . "/$case_name" . "_final_filtered_coding_clean.xls"; #OUT6

  open(ANNO, $grand_anno_count_file) || die $self->error_message("could not open grand anno read counts file: $grand_anno_count_file");
  open(OUT1, ">$final_unfiltered_tsv") || die $self->error_message("could not open output file: $final_unfiltered_tsv");
  open(OUT2, ">$final_filtered_tsv") || die $self->error_message("could not open output file: $final_filtered_tsv");
  open(OUT3, ">$final_filtered_clean_tsv") || die $self->error_message("could not open output file: $final_filtered_clean_tsv");
  open(OUT4, ">$final_filtered_coding_clean_tsv") || die $self->error_message("could not open output file: $final_filtered_coding_clean_tsv");

  my @skip = qw (gene_name transcript_species transcript_source transcript_version transcript_status c_position ucsc_cons domain all_domains deletion_substructures transcript_error gene_name_source);
  my %skip_columns;
  foreach my $name (@skip){
    $skip_columns{$name}=1;
  }
  my @include_col_pos;

  my $header = 1;
  my %columns;
  while(<ANNO>){
    chomp($_);
    my @line = split("\t", $_);
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
      my $header_extension = "min_coverage_observed\tmax_normal_vaf_observed\tmax_tumor_vaf_observed\tvariant_source_callers\tvariant_source_caller_count\tfiltered";
      my $full_header = "$_"."\t$header_extension";
      $full_header .= "\t$per_lib_header" if $per_lib_header;
      print OUT1 "$full_header\n";
      print OUT2 "$full_header\n";

      my @include_values = @line[@include_col_pos];
      my $include_values_string = join("\t", @include_values);
      my $short_header = "$include_values_string"."\t$header_extension";
      $short_header .= "\t$per_lib_header" if $per_lib_header;
      print OUT3 "$short_header\n";
      print OUT4 "$short_header\n";
      next;
    }

    my ($chr, $start, $stop, $ref, $var) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
    my $v = $chr . "_$start" . "_$stop" . "_$ref" . "_$var";
    die $self->error_message("parsed a variant that is not defined in the variant hash") unless $variants->{$v};

    my @per_lib_counts = @{$variants->{$v}->{per_lib_counts}} if defined($variants->{$v}->{per_lib_counts});
    my $per_lib_count_line = join("\t", @per_lib_counts) if defined($variants->{$v}->{per_lib_counts});

    my $line_extension = "$variants->{$v}->{min_coverage_observed}\t$variants->{$v}->{max_normal_vaf_observed}\t$variants->{$v}->{max_tumor_vaf_observed}\t$variants->{$v}->{variant_source_callers}\t$variants->{$v}->{variant_source_caller_count}\t$variants->{$v}->{filtered}";  
    my $full_line = "$_\t$line_extension";
    $full_line .= "\t$per_lib_count_line" if defined($per_lib_count_line);

    print OUT1 "$full_line\n";
    unless ($variants->{$v}->{filtered}){
      print OUT2 "$full_line\n";
    }

    my @include_values = @line[@include_col_pos];
    my $include_values_string = join("\t", @include_values);
    my $short_line = "$include_values_string"."\t$line_extension";
    $short_line .= "\t$per_lib_count_line" if defined($per_lib_count_line);

    unless ($variants->{$v}->{filtered}){
      print OUT3 "$short_line\n";
    }

    #apply a transcript variant type filter to define the 'coding' result
    my $trv_type = $variants->{$v}->{trv_type};
    unless ($variants->{$v}->{filtered}){
      unless ($trv_type_filter =~ /$trv_type/i){
        print OUT4 "$short_line\n";
      }
    }

  }
  close(ANNO);
  close(OUT1);
  close(OUT2);
  close(OUT3);

  #Run the R scripts to generate some plots 
  my $r_script = __FILE__ . '.R';
  my $r_cmd = "$r_script $final_filtered_coding_clean_tsv \"normal_day0_VAF tumor_day0_VAF tumor_day30_VAF\" \"34 37 40\" \"43 46 49\" \"52 55 58\" " . $self->target_gene_list_name . " " . $self->outdir;
  print Dumper $r_cmd;
  Genome::Sys->shellcmd(cmd => $r_cmd);


  return;
}

1;

