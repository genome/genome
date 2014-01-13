package Genome::Model::ClinSeq::Command::Converge::SnvIndelReport;
use strict;
use warnings;
use Genome;
use Data::Dumper;
use List::MoreUtils qw/ uniq /;

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
  print "\nUsing annotation build: $annotation_build_name";
  my $result = $self->gather_variants('-somatic_builds'=>$somatic_builds);
  my $variants = $result->{'variants'};
  my $header = $result->{'header'};

  #Make note of which variants correspond to $self->target_gene_list
  $self->intersect_target_gene_list('-variants'=>$variants);

  #Print a grand table of all annotated variants with extra annotation columns (tier, target gene list)
  my $grand_anno_file = $self->print_grand_anno_table('-variants'=>$variants, '-header'=>$header);

  #Identify underlying reference-alignments, sample names, sample common names, and timepoints (if available)
  my $align_builds = $self->get_ref_align_builds('-somatic_builds'=>$somatic_builds);

  #Get bam-readcounts for all positions for all BAM files
  #TODO: This will need to be done at the level of replicate libraries once that readcounting tool is availble
  #TODO: For now just use standard bam read counts
  $self->add_read_counts('-align_builds'=>$align_builds, '-grand_anno_file'=>$grand_anno_file);
  
  #Apply some automatic filters
  #Normal VAF > $self->max_normal_vaf
  #Min coverage > $self->min_coverage (position must be covered at this level across all samples)

  #Make note of which variants lie within a particular set of regions of interest (e.g. nimblegen v3 + AML RMG)
  #Have option to filter these out

  #Filter out variants falling within false positive regions/genes?

  #Remove variants with disallowed amino acid effects (e.g. silent, rna, null, etc.) as defined in $self->effects_filter

  #Write out a tsv file and an Excel spreadsheet file

  print "\n\n";

  return 1;
};


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
  mkdir ($bed_dir);

  #Make a copy of each desired .bed file.  Perform basic clean-ups and sanity checks
  foreach my $file (sort keys %bed_files){
    my $c = 0;
    my $var_type = $bed_files{$file}{var_type};
    my $data_type = $bed_files{$file}{data_type};
    my $tier = $bed_files{$file}{tier};
    my $somatic_build_id = $bed_files{$file}{somatic_build_id};
    my $new_file = $bed_dir . "$somatic_build_id" . "_$data_type" . "_$var_type" . "_$tier" . "clean.bed";
    $self->status_message("\nProcessing $var_type $data_type $tier file: $file"); 
    
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
    $bed_files{$file}{clean_bed_file} = $new_file;

    #Run the annnotator on each bed file and add RSIDs and GMAFs for all variants
    #gmt annotate transcript-variants --variant-bed-file='' --output-file='' --annotation-filter=top --reference-transcripts='NCBI-human.ensembl/74_37'
    my $anno_file = $bed_files{$file}{clean_bed_file} . ".anno";
    my $annotate_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
          variant_bed_file=>$bed_files{$file}{clean_bed_file}, 
          output_file=>$anno_file,
          annotation_filter=>'top',
          reference_transcripts=>$self->annotation_build->name,
        );
    my $r1 = $annotate_cmd->execute();
    die $self->error_message("annotation cmd unsuccessful") unless ($r1);
    $bed_files{$file}{anno_file} = $anno_file;

    #Get RSIDs and GMAFs for all variants 
    #gmt annotate add-rsid --anno-file='' --output-file='' --vcf-file=''
    my $vcf_file = $bed_files{$file}{vcf_file};
    my $rsid_file = $anno_file . ".rsid";
    my $rsid_cmd = Genome::Model::Tools::Annotate::AddRsid->create(
          anno_file=>$anno_file,
          vcf_file=>$vcf_file,
          output_file=>$rsid_file,
        );
    my $r2 = $rsid_cmd->execute();
    die $self->error_message("rsid cmd unsuccessful") unless ($r2);
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
      my $v = $chr . "_$start" . "_$stop" . "_$ref" . "_$var";
      $variants{$v}{anno_line} = $_;
      $variants{$v}{tier} = $tier;
      $variants{$v}{ensembl_gene_id} = $ensembl_gene_id;
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
  open (ANNO, ">$grand_anno_file") || die $self->error_message("could not open grand anno file: $grand_anno_file");
  my $target_gene_list_name = $self->target_gene_list_name;
  my $new_header = "$header\ttier\t$target_gene_list_name";
  print ANNO "$new_header\n";
  foreach my $v (sort keys %{$variants}){
    print ANNO "$variants->{$v}->{anno_line}\t$variants->{$v}->{tier}\t$variants->{$v}->{target_gene_list_match}\n";
  }
  close(ANNO);

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
      $variants->{$v}->{target_gene_list_match} = 0;
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
  my $header_prefixes;
  my @unique_time_points = uniq @time_points;
  my @unique_samples = uniq @samples;
  my @unique_names = uniq @names;
  if (scalar(@unique_time_points) == scalar(@time_points)){
    $header_prefixes = join(",", @time_points);
  }elsif(scalar(@unique_samples) == scalar(@samples)){
    $header_prefixes = join(",", @samples);
  }elsif(scalar(@unique_names) == scalar(@names)){
    $header_prefixes = join(",", @names);
  }else{
    die $self->error_message("could not resolve unique prefixes for add-readcounts");
  }

  #gmt analysis coverage add-readcounts --bam-files=? --genome-build=? --output-file=? --variant-file=? [--header-prefixes=?] 
  my $output_file = $self->outdir . "/variants.all.anno.readcounts";
  my $add_count_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
          bam_files=>$bam_list,
          genome_build=>$reference_fasta,
          output_file=>$output_file,
          variant_file=>$grand_anno_file,
          header_prefixes=>$header_prefixes,
        );
  my $r = $add_count_cmd->execute();
  die $self->error_message("add-readcounts cmd unsuccessful") unless ($r);

  return ($output_file);
}


1;

