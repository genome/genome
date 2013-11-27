package Genome::Model::ClinSeq::Command::GetBamReadCountsMatrix;

#Written by Malachi Griffith and Scott Smith

#Load modules
use strict;
use warnings;
use Genome; 
use Data::Dumper;
use Genome;
use Genome::Info::IUB; 

class Genome::Model::ClinSeq::Command::GetBamReadCountsMatrix {
    is => 'Command::V2',
    has => [
        output_dir              => { is => 'FilesystemPath',
                                     doc => 'Directory where all output files will be written', },

        somatic_build_ids       => { is => 'Text',
                                     is_optional => 1,
                                     doc => 'Comma separated list of somatic variation builds.  Will be used to build a list of SNV positions to query', },

        ref_align_build_ids     => { is => 'Text',
                                     is_optional => 1,
                                     doc => "Comma separated list of reference alignment builds for which BAM readcounts will be generated", },

        somatic_labels          => { is => 'Text',
                                     is_optional => 1,
                                     doc => "Comma separated list of labels that correspond to the somatic variation builds (build ID will be used if not specified)", },

        ref_align_labels        => { is => 'Text',
                                     is_optional => 1,
                                     doc => "Comma separated list of labels that correspond to the reference alignment builds (build ID will be used if not specified)", },

        positions_file          => { is => 'FilesystemPath',
                                     is_optional => 1,
                                     doc => "File containing SNV positions of interest in the BED like format: \nchr\tstart\tend\tgenotype\ne.g. 1 100 101 C/S", },

        target_file_names       => { is => 'Text',
                                     is_optional => 1,
                                     doc => "Comma separated list of file names to gather SNV positions from (by default, all snv.hq files in the effects dir will be used).  e.g. snvs.hq.novel.tier1.v2.bed", },

        skip_mt                 => { is => 'Number', 
                                     is_optional => 1,
                                     doc => 'Whole genome sequence (WGS) somatic variation build', },

        no_fasta_check          => { is => 'Number',
                                     is_optional => 1,
                                     doc => 'Do not check that ref base matches expected in reference fasta', },
        max_positions           => { is => 'Number',
                                     is_optional => 1,
                                     doc => 'For debugging/testing purposes, limit counting to this many positions',
                                   },

    ],
    doc => 'This script attempts to get read counts for an arbitrary set of positions from an arbitrary list of reference alignment BAMs',
};


sub sub_command_category { 'pipeline' }

sub help_detail {
    return <<EOS

 This script takes a list of somatic variation builds, gathers the unique set of SNV positions for those and 
 obtains BAM read counts for the underlying reference alignment builds or an explicitly provided list of them.
 Note: Do NOT use for Indels!  SNVs only.
 Note: The output will be a matrix file for read counts and one for VAFs.  One column will be written for each somatic variation build ID specified.

EOS
}

sub help_synopsis {
  return <<EOS
  genome model clin-seq get-bam-read-counts-matrix \
    --somatic-build-ids='126851693,126847207' \
    --somatic-labels='DefaultSomatic_PNC4,StrelkaSomatic_PNC4' \
    --output-dir='/tmp/bam_readcount_matrix/' \ 
    --skip-mt=1
EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless ($self->somatic_build_ids || ($self->positions_file && $self->ref_align_build_ids)) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => [qw/somatic_build_ids positions_file ref_align_build_ids/],
	                                          desc => 'must define either a list of somatic variation builds or supply a positions file in the correct form along with reference alignment builds',
                                          );
  }
  unless (-e $self->output_dir && -d $self->output_dir){
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['output_dir'],
	                                          desc => "output directory: ". $self->output_dir . " does not appear to be valid",
                                          );
  }
  if ($self->somatic_build_ids && $self->somatic_labels){
    my @b = split(",", $self->somatic_build_ids);
    my @l = split(",", $self->somatic_labels);
    my $b_count = scalar(@b);
    my $l_count = scalar(@l); 
    unless ($b_count == $l_count){
     push @errors, UR::Object::Tag->create(
                                            type => 'error',
                                            properties => [qw/somatic_build_ids somatic_labels/],
                                            desc => "number of somatic build ids and labels must match",
                                          );
    }
  }

  if ($self->ref_align_build_ids && $self->ref_align_labels){
    my @b = split(",", $self->ref_align_build_ids);
    my @l = split(",", $self->ref_align_labels);
    my $b_count = scalar(@b);
    my $l_count = scalar(@l);
    unless ($b_count == $l_count){
     push @errors, UR::Object::Tag->create(
                                            type => 'error',
                                            properties => [qw/ref_align_build_ids ref_align_labels/],
                                            desc => "number of ref align build ids and labels must match",
                                          );
    }
  }

  return @errors;
}

sub help_usage {
    my $self = shift;
    my $usage = $self->SUPER::help_usage(@_);
    $self->status_message("$usage");
}

sub execute {
  my $self = shift;
  
  eval "require Bio::DB::Sam";
  if ($@) {
      die "Failed to use the Bio::DB::Sam module.  Use /usr/bin/perl (5.10 or greater!) instead of /gsc/bin/perl.:\n$@";
  }

  #Get the target file names to be examined
  my @target_file_names = $self->get_target_file_names;

  #Get all SNV positions and write to a positions file
  my $pos;
  if ($self->positions_file && $self->somatic_build_ids){
    $self->status_message("Warning: user specified positions file is overriding the list that would be generated by examining the somatic_build_ids");
    $pos = $self->format_positions_file;
  }elsif($self->positions_file){
    $pos = $self->format_positions_file;
  }else{
    $pos = $self->write_positions_file('-target_file_names'=>\@target_file_names);
  }

  my $pos_count = keys %{$pos};
  $self->status_message("Found $pos_count unique positions to obtain counts for");

  #If the user defined --max-positions, thin the list of target positions to this max
  if ($self->max_positions){
    $self->warning_message("Limiting counting to " . $self->max_positions . " positions as specified by --max_positions");
    my $c = 0;
    foreach my $p (sort keys %{$pos}){
      $c++;
      if ($c > $self->max_positions){
        delete $pos->{$p};
      }
    }
  }

  $pos_count = keys %{$pos};
  $self->status_message("\tProceeding to obtain counts for $pos_count positions");

  #Get reference alignment builds (from ref_align_build_ids or somatic_build_ids) and associated BAM files to be counted
  my $ref_align_builds = $self->get_ref_align_builds;

  #Get the reference fasta file from the underlying builds (ref_align_build_ids if they are available otherwise somatic_build_ids)
  my $ref_fasta = $self->get_ref_fasta('-ref_align_builds'=>$ref_align_builds);

  #Define a code reference needed for Bio::DB::Sam
  #Code reference needed for Bio::DB::Bam
  my $callback = sub {
    my ($tid,$pos,$pileups,$callback_data) = @_;
    my $data = $callback_data->[0];
    my $read_counts = $callback_data->[1];
    my $fai = $callback_data->[2];
    if ( ($pos == ($data->{start} - 1) ) ) {
      unless ($self->no_fasta_check){
        my $ref_base = $fai->fetch($data->{chr} .':'. $data->{start} .'-'. $data->{stop});
        unless ($data->{reference} eq $ref_base) {
          die("\n\nReference base " . $ref_base .' does not match '. $data->{reference} .' at '. $pos .' for chr '. $data->{chr} . '(tid = '. $tid . ')');
        }
      }
      for my $pileup ( @{$pileups} ) {
        my $alignment = $pileup->alignment;
        next if $pileup->indel or $pileup->is_refskip;
        my $qbase  = substr($alignment->qseq,$pileup->qpos,1);
        next if $qbase =~ /[nN]/;
        $read_counts->{$qbase}++;
      }
    }
  };


  #Perform the count on each BAM file for all positions in the positions file and write the output files
  $self->status_message("Getting BAM read counts for each alignment build");
  foreach my $bid (sort {$ref_align_builds->{$a}->{label} cmp $ref_align_builds->{$b}->{label}} keys %{$ref_align_builds}){
    my $label = $ref_align_builds->{$bid}->{label};
    my $bam_path = $ref_align_builds->{$bid}->{bam_file};
    $self->status_message("Processing: $label");

    #Get Bio:DB:Bam objects for this BAM file and reference fasta file
    my $bam = Bio::DB::Bam->open($bam_path);
    my $header = $bam->header;
    my $index = Bio::DB::Bam->index($bam_path);
    my $fai = Bio::DB::Sam::Fai->load($ref_fasta);

    foreach my $p (sort keys %{$pos}){
      my %data;
      my $data = \%data;
      $data->{chr} = $pos->{$p}->{chr};
      $data->{start} = $pos->{$p}->{start} + 1; #Add 1 because we are starting with BED format
      $data->{stop} = $pos->{$p}->{end}; 
      $data->{reference} = $pos->{$p}->{ref_base};
      $data->{variant} = $pos->{$p}->{var_base};
      my $seq_id = $data->{chr} .':'. $data->{start} .'-'. $data->{stop};
      my ($tid,$start,$end) = $header->parse_region($seq_id);
      #print "\n\nlabel: $label\tseq_id: $seq_id\ttid: $tid\tstart: $start\tend: $end";

      my %read_counts;
      #print "\n\ntid: $tid\tstart: $start\tend: $end\tref_base: $data->{reference}\tvar_base: $data->{variant}";
      $index->pileup($bam,$tid,$start,$end,$callback,[$data,\%read_counts, $fai]);
      $data->{A} = $read_counts{A} || 0;
      $data->{T} = $read_counts{T} || 0;
      $data->{C} = $read_counts{C} || 0;
      $data->{G} = $read_counts{G} || 0;
      #print "\n\tA: $data->{A}\tT: $data->{T}\tC: $data->{C}\tG: $data->{G}";

      my $total_rc =  $data->{A} + $data->{T} + $data->{C} + $data->{G};
      my $ref_rc = $data->{$pos->{$p}->{ref_base}};
      my $var_rc = $data->{$pos->{$p}->{var_base}};
      my $var_allele_frequency = 0;
      if ($total_rc){
        $var_allele_frequency = sprintf ("%.3f", (($var_rc / $total_rc)*100));
      }
      #print "\n\t\tTotalCount: $total_rc\tRefCount: $ref_rc\tVarCount: $var_rc\tVAF: $var_allele_frequency%";

      #Store resulting counts on the positions hash with a secondary key on the alignment build id
      $pos->{$p}->{$bid}->{rc} = $total_rc;
      $pos->{$p}->{$bid}->{vaf} = $var_allele_frequency;
    }
  }

  #Print out the final matrix files for total read count and VAF
  my @label_list;
  foreach my $bid (sort {$ref_align_builds->{$a}->{label} cmp $ref_align_builds->{$b}->{label}} keys %{$ref_align_builds}){
    my $label = $ref_align_builds->{$bid}->{label};
    push(@label_list, $label);
  }
  my $label_string = join ("\t", @label_list);

  my $header = "coord\tchr\tstart\tend\tgenotype\tref_base\tvar_base\t$label_string";

  my $rc_file_path = $self->output_dir . "/All_ReadCounts.tsv";
  my $vaf_file_path = $self->output_dir . "/All_VAFs.tsv";
  open (RC, ">$rc_file_path") || die "\n\nCould not open read counts outfile: $rc_file_path\n\n";
  open (VAF, ">$vaf_file_path") || die "\n\nCould not open VAF outfile: $vaf_file_path\n\n";
  print RC "$header\n";
  print VAF "$header\n";

  foreach my $p (sort keys %{$pos}){
    my @rc_list;
    my @vaf_list;
    foreach my $bid (sort {$ref_align_builds->{$a}->{label} cmp $ref_align_builds->{$b}->{label}} keys %{$ref_align_builds}){
      push(@rc_list, $pos->{$p}->{$bid}->{rc});
      push(@vaf_list, $pos->{$p}->{$bid}->{vaf});
    }
    my $rc_string = join("\t", @rc_list);
    my $vaf_string = join("\t", @vaf_list);

    my $annotation = "$p\t$pos->{$p}->{chr}\t$pos->{$p}->{start}\t$pos->{$p}->{end}\t$pos->{$p}->{genotype}\t$pos->{$p}->{ref_base}\t$pos->{$p}->{var_base}";

    print RC "$annotation\t$rc_string\n";
    print VAF "$annotation\t$vaf_string\n";
  }
  
  close(RC);
  close(VAF);

  $self->status_message("\n\nCOMPLETE\n\n");

  return 1;

}


#Get a list of target file names for SNV files in the effects dir of a somatic variation build dir
sub get_target_file_names{
  my $self = shift;
  my @file_names;
  if ($self->target_file_names){
    @file_names = split(",", $self->target_file_names);
  }else{
    #If not specified by the user, use all snvs.hq files that are expected in the 'effects' dir of a somatic variation build dir
    @file_names = qw (
                      snvs.hq.novel.tier1.v2.bed
                      snvs.hq.novel.tier2.v2.bed
                      snvs.hq.novel.tier3.v2.bed
                      snvs.hq.novel.tier4.v2.bed
                      snvs.hq.previously_detected.tier1.v2.bed
                      snvs.hq.previously_detected.tier2.v2.bed
                      snvs.hq.previously_detected.tier3.v2.bed
                      snvs.hq.previously_detected.tier4.v2.bed
                     );
  }
  return @file_names;
}

#If a list of somatic build IDs was supplied and a positions file was not, pull the unique union of all SNVs from the specified somatic variation builds
sub write_positions_file{
  my $self = shift;
  my %args = @_;
  my @target_file_names = @{$args{'-target_file_names'}};

  my @somatic_build_ids = split(",", $self->somatic_build_ids);

  my @somatic_labels = @somatic_build_ids;
  if ($self->somatic_labels){
    @somatic_labels = split(",", $self->somatic_labels);
  }

  my %pos;

  my $out_header = "coord\tchr\tstart\tend\tgenotype\tref_base\tvar_base";

  foreach my $somatic_build_id (@somatic_build_ids){
    my $somatic_label = shift @somatic_labels;
    my $somatic_file_path = $self->output_dir . "/" . $somatic_label . "_snvs.tsv";
    open (OUT, ">$somatic_file_path") || die "\n\nCould not open $somatic_file_path for writing\n\n";
    print OUT "$out_header\n";

    my $somatic_build = Genome::Model::Build->get($somatic_build_id);
    unless ($somatic_build){
      $self->error_message("Could not obtain somatic variation build object for specified id: $somatic_build_id");
      exit(1);
    }
    my $data_directory = $somatic_build->data_directory;
    my $effects_dir = $data_directory . "/effects/";
    unless (-e $effects_dir){
      $self->error_message("Could not find effects dir: $effects_dir for specified build id: $somatic_build_id");
      exit(1);
    }

    foreach my $target_file_name (@target_file_names){
      my $snv_file = $effects_dir . $target_file_name;
      unless (-e $snv_file){
        $self->error_message("Could not SNV file: $snv_file");
        exit(1);
      }

      open (SNV, "$snv_file") || die "\n\nCould not open SNV file: $snv_file\n\n";
      while (<SNV>){
        chomp($_);
        my @line = split("\t", $_);
        my $chr = $line[0];
        my $start = $line[1];
        my $end = $line[2];
        my $genotype = $line[3];

        if ($self->skip_mt && $chr =~ /^mt$|^m$/i){
          next();
        }

        my $ref_base;
        my $var_iub;
        if ($genotype =~ /(\w+)\/(\w+)/){
          $ref_base = $1;
          $var_iub = $2;
        }else{
          $self->error_message("Genotype in file: $snv_file not understood: $genotype");
          exit(1);
        }
        my $var_bases = Genome::Info::IUB->iub_to_string($var_iub);
        my $var_base;
        if ($var_bases =~ /^(\w)(\w)$/){
          if ($1 eq $ref_base){
            $var_base = $2;
          }else{
            $var_base = $1;
          }
        }else{
          $self->error_message("Could not interpret var bases from iub_to_string: $var_bases");
          exit(1);
        }
        #$self->status_message("$chr\t$start\t$end\t$genotype\t$ref_base\t$var_iub\t$var_bases\t$var_base");
        
        #Store this coord.  Key on chr:position and only store one instance of each (i.e. mutltiple allele will be collapsed)
        my $coord = "$chr:$start-$end";
        $pos{$coord}{chr} = $chr;
        $pos{$coord}{start} = $start;
        $pos{$coord}{end} = $end;
        $pos{$coord}{genotype} = $genotype;
        $pos{$coord}{ref_base} = $ref_base;
        $pos{$coord}{var_iub} = $var_iub;
        $pos{$coord}{var_base} = $var_base;
        print OUT "$coord\t$pos{$coord}{chr}\t$pos{$coord}{start}\t$pos{$coord}{end}\t$pos{$coord}{genotype}\t$pos{$coord}{ref_base}\t$pos{$coord}{var_base}\n";
      }
      close (SNV);
    }
    close (OUT);
  }

  my $file_path = $self->output_dir . "/All_snvs.tsv";
  open (OUT, ">$file_path") || die "\n\nCould not open $file_path for writing\n\n";
  print OUT "$out_header\n";
  foreach my $coord (sort keys %pos){
    print OUT "$coord\t$pos{$coord}{chr}\t$pos{$coord}{start}\t$pos{$coord}{end}\t$pos{$coord}{genotype}\t$pos{$coord}{ref_base}\t$pos{$coord}{var_base}\n";
  }
  close(OUT);

  return(\%pos);
}

sub format_positions_file{
  my $self = shift;

  my %pos;

  my $out_header = "coord\tchr\tstart\tend\tgenotype\tref_base\tvar_base";
  my $snv_file = $self->positions_file;
      
  unless (-e $snv_file){
    $self->error_message("Could not SNV file: $snv_file");
    exit(1);
  }

  open (SNV, "$snv_file") || die "\n\nCould not open SNV file: $snv_file\n\n";
  while (<SNV>){
    chomp($_);
    my @line = split("\t", $_);
    my $chr = $line[0];
    my $start = $line[1];
    my $end = $line[2];
    my $genotype = $line[3];

    if ($self->skip_mt && $chr =~ /^mt$|^m$/i){
      next();
    }

    my $ref_base;
    my $var_iub;
    if ($genotype =~ /(\w+)\/(\w+)/){
      $ref_base = $1;
      $var_iub = $2;
    }else{
      $self->error_message("Genotype in file: $snv_file not understood: $genotype");
      exit(1);
    }
    my $var_bases = Genome::Info::IUB->iub_to_string($var_iub);
    my $var_base;
    if ($var_bases =~ /^(\w)(\w)$/){
      if ($1 eq $ref_base){
        $var_base = $2;
      }else{
        $var_base = $1;
      }
    }else{
      $self->error_message("Could not interpret var bases from iub_to_string: $var_bases");
      exit(1);
    }
        
    #Store this coord.  Key on chr:position and only store one instance of each (i.e. mutltiple allele will be collapsed)
    my $coord = "$chr:$start-$end";
    $pos{$coord}{chr} = $chr;
    $pos{$coord}{start} = $start;
    $pos{$coord}{end} = $end;
    $pos{$coord}{genotype} = $genotype;
    $pos{$coord}{ref_base} = $ref_base;
    $pos{$coord}{var_iub} = $var_iub;
    $pos{$coord}{var_base} = $var_base;
  }
  close (SNV);

  my $file_path = $self->output_dir . "/All_snvs.tsv";
  open (OUT, ">$file_path") || die "\n\nCould not open $file_path for writing\n\n";
  print OUT "$out_header\n";
  foreach my $coord (sort keys %pos){
    print OUT "$coord\t$pos{$coord}{chr}\t$pos{$coord}{start}\t$pos{$coord}{end}\t$pos{$coord}{genotype}\t$pos{$coord}{ref_base}\t$pos{$coord}{var_base}\n";
  }
  close(OUT);

  return(\%pos);
}

sub get_ref_align_builds{
  my $self = shift;

  my %builds;

  if ($self->ref_align_build_ids){
    my @build_ids = split(",", $self->ref_align_build_ids);
    my @labels = @build_ids;
    if ($self->ref_align_labels){
      @labels = split(",", $self->ref_align_labels);
    }
    foreach my $build_id (@build_ids){
      my $ref_align_build = Genome::Model::Build->get($build_id);
      unless ($ref_align_build){
        $self->error_message("Could not obtain build object from ref align build id: $build_id");
        exit(1);
      }
      my $label = shift @labels;
      $builds{$build_id}{build} = $ref_align_build;
      $builds{$build_id}{label} = $label;
    }
  }elsif($self->somatic_build_ids){
    my @somatic_build_ids = split(",", $self->somatic_build_ids);
    #Note that multiple somatic builds (with different called SNVs) could lean on the same ref align builds
    foreach my $somatic_build_id (@somatic_build_ids){
      my $somatic_build = Genome::Model::Build->get($somatic_build_id);
      my $normal_build = $somatic_build->normal_build;
      my $normal_build_id = $normal_build->id;
      my $normal_label = $normal_build_id . "_normal";
      my $tumor_build = $somatic_build->tumor_build;
      my $tumor_build_id = $tumor_build->id;
      my $tumor_label = $tumor_build_id . "_tumor";
      my $common_name = '';
      my $subject = $somatic_build->subject;
      my $patient = $subject->patient;
      if ($patient){
        if ($patient->can("common_name")){
          if ($patient->common_name){
            $common_name = $patient->common_name;
            $normal_label = $common_name . "_" . $normal_label;
            $tumor_label = $common_name . "_" . $tumor_label;
          }
        }
        unless ($common_name){
          my $normal_subject_name = $normal_build->model->subject->name;
          my $tumor_subject_name = $tumor_build->model->subject->name;
          $normal_label = $normal_subject_name . "_" . $normal_label;
          $tumor_label = $tumor_subject_name . "_" . $tumor_label;
        }
      }
      $builds{$normal_build_id}{build} = $normal_build;
      $builds{$normal_build_id}{label} = $normal_label;
      $builds{$tumor_build_id}{build} = $tumor_build;
      $builds{$tumor_build_id}{label} = $tumor_label;
    }
  }

  #Get the BAM file for each alignment build
  foreach my $build_id (keys %builds){
    my $build = $builds{$build_id}{build};
    my $bam_file;
    if ($build->can("whole_rmdup_bam_file")){
      $bam_file = $build->whole_rmdup_bam_file;
    }elsif($build->can("alignment_result")){
      my $alignment_result = $build->alignment_result;
      if ($alignment_result->can("bam_file")){
        $bam_file = $alignment_result->bam_file;
      }
    }
    if ($bam_file){
      if (-e $bam_file){
        $builds{$build_id}{bam_file} = $bam_file;
      }else{
        $self->error_message("Found BAM file path but BAM file is missing for build: $build");
        exit(1);
      }
    }else{
      $self->error_message("Could not find BAM file for build: $build");
      exit(1);
    }
  }

  return \%builds;
}

sub get_ref_fasta{
  my $self = shift;
  my %args = @_;
  my $ref_align_builds = $args{'-ref_align_builds'};
  my $ref_fasta_path;

  my %refs;
  foreach my $build_id (keys %{$ref_align_builds}){
    my $ref_align_build = $ref_align_builds->{$build_id}->{build};
    my $m = $ref_align_build->model;
    if ($m->can("reference_sequence_build")){
      my $rb = $m->reference_sequence_build;
      my $rb_name = $rb->name;
      $ref_fasta_path = $rb->full_consensus_path('fa');
      $refs{$rb_name}=1;
    }
  }
  my $rb_count = keys %refs;
  if ($rb_count > 1){
    $self->error_message("Found multiple reference fasta builds, aborting:");
    print Dumper %refs;
    exit(1);
  }
  if ($rb_count == 0){
    $self->error_message("Could not find the reference fasta path, aborting");
    exit(1);
  }

  return $ref_fasta_path;
}



