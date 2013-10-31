Genome::Model::ClinSeq::Command::GetBamReadCounts->class; #loads the class
package Genome::Model::ClinSeq::Command::GetBamReadCounts;

#Written by Malachi Griffith and Scott Smith

#Load modules
use strict;
use warnings;
use Genome; 
use Data::Dumper;
use Genome;
use Genome::Model::ClinSeq::Util qw(:all);
use Genome::Model::ClinSeq::RnaSeqUtil qw(:all);

class Genome::Model::ClinSeq::Command::GetBamReadCounts {
    is => 'Command::V2',
    has => [
        positions_file          => { is => "FilesystemPath",
                                     doc => "File containing SNV positions of interest and ref/var bases\n"
                                            . "  (e.g. 5:112176318-112176318	APC	APC	p.R1676T	G	C)" },

        wgs_som_var_build       => { is => 'Genome::Model::Build::SomaticVariation', is_optional => 1,
                                     doc => 'Whole genome sequence (WGS) somatic variation build' },

        exome_som_var_build     => { is => 'Genome::Model::Build::SomaticVariation', is_optional => 1,
                                     doc => 'Exome capture sequence somatic variation build' },
        
        rna_seq_normal_build    => { is => "Genome::Model::Build", is_optional => 1,
                                     doc => "RNA-seq model id for normal" },

        rna_seq_tumor_build     => { is => "Genome::Model::Build", is_optional => 1,
                                     doc => 'RNA-seq model id for tumor' },

        output_file             => { is => 'FilesystemPath',
                                     doc => 'File where output will be written (input file values with read counts appended)', },

        verbose                 => { is => 'Number', is_optional => 1,
                                     doc => 'To display more output, set this to 1.' },

        no_fasta_check          => { is => 'Number', is_optional => 1,
                                     doc => 'To prevent checking of the reported reference base against the reference genome fasta set --no_fasta_check=1 [Not recommended!]' },

        cancer_annotation_db    => { is => 'Genome::Db::Tgi::CancerAnnotation',
                                     doc => 'database of cancer annotation' },

    ],
    doc => 'This script attempts to get read counts, frequencies and gene expression values for a series of genome positions',
};


sub sub_command_category { 'pipeline' }

sub help_detail {
    return <<EOS
 This script takes an input file with SNV positions and determines: reference and variant allele read counts, frequencies, RNA-seq gene expression levels
 Up to three pairs of BAMs can be specified (WGS tumor+normal, Exome tumor+normal, RNA tumor+normal)
 Results will be appended as new columns in the input file.
 Notes: Do NOT use for Indels!  SNVs only.
EOS
}

sub help_synopsis {
  return <<EOS
  genome model clin-seq get-bam-read-counts \
    --positions-file=snvs.hq.tier1.v1.annotated.compact.tsv \
    --wgs-som-var-build='2880644349' \
    --exome-som-var-build='2880732183' \
    --rna-seq-tumor-build='2880693923' \
    --output-file=snvs.hq.tier1.v1.annotated.compact.readcounts.tsv
EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless (($self->wgs_som_var_build || $self->exome_som_var_build || $self->rna_seq_normal_build || $self->rna_seq_tumor_build)) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => [qw/wgs_som_var_build exome_som_var_build rna_seq_normal_build rna_seq_tumor_build/],
	                                          desc => 'at least one of the four build types must be specified!'
                                          );
  }
  unless (-e $self->positions_file) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['positions_file'],
	                                          desc => "Positions file: " . $self->positions_file . " not found",
      );
  }
  return @errors;
}

sub help_usage {
    my $self = shift;
    my $usage = $self->SUPER::help_usage(@_);
    return $usage;
}

sub execute {
  my $self = shift;
  
  eval "require Bio::DB::Sam";
  if ($@) {
      die "Failed to use the Bio::DB::Sam module.  Use /usr/bin/perl (5.10 or greater!) instead of /gsc/bin/perl.:\n$@";
  }

  my $positions_file = $self->positions_file; 
  my $wgs_som_var_build = $self->wgs_som_var_build;
  my $exome_som_var_build = $self->exome_som_var_build;
  my $rna_seq_normal_build = $self->rna_seq_normal_build;
  my $rna_seq_tumor_build = $self->rna_seq_tumor_build;
  my $no_fasta_check = $self->no_fasta_check;
  my $output_file = $self->output_file;
  my $verbose = $self->verbose;
  my $cancer_annotation_db = $self->cancer_annotation_db;

  #Build a map of ensembl transcript ids to gene ids and gene names
  my $annotation_build;
  my $reference_build;
  my @builds = ($wgs_som_var_build, $exome_som_var_build, $rna_seq_normal_build, $rna_seq_tumor_build);
  foreach my $build (@builds){
    if ($build){
      $reference_build = $build->reference_sequence_build;
      my $model = $build->model;
      if ($model->can("annotation_reference_build")){
        $annotation_build = $model->annotation_reference_build;
      }
      if ($model->can("annotation_build")){
        $annotation_build = $model->annotation_build;
      }
    }
  }
  $self->error_message("Could not resolve annotation build from input builds") unless $annotation_build;
  die $self->error_message unless $annotation_build;
  $self->error_message("Could not resolve reference sequence build from input builds") unless $reference_build;
  die $self->error_message unless $reference_build;

  my $reference_build_id = $reference_build->id;
  my $annotation_build_name = $annotation_build->name;
  my $annotation_data_dir = $annotation_build->data_directory;
  my $transcript_info_path = $annotation_data_dir . "/annotation_data/rna_annotation/$reference_build_id-transcript_info.tsv";
  my $gtf_path = $annotation_build->annotation_file('gtf',$reference_build_id);
  unless (defined($gtf_path)) {
    die $self->error_message("'There is no annotation GTF file defined for annotation_reference_transcripts build: ". $annotation_build->__display_name__);
  }
  unless (-e $transcript_info_path) {
    die $self->error_message("'There is no transcript info file for annotation_reference_transcripts build: ". $annotation_build->__display_name__);
  }

  my $ensembl_map = &loadEnsemblMap('-gtf_path'=>$gtf_path, '-transcript_info_path'=>$transcript_info_path);

  #Get Entrez and Ensembl data for gene name mappings
  my $entrez_ensembl_data = &loadEntrezEnsemblData(-cancer_db => $cancer_annotation_db);

  #Import SNVs from the specified file
  my $result = &importPositions('-positions_file'=>$positions_file);
  my $snvs = $result->{'snvs'};
  my $snv_header = $result->{'header'};
  #print Dumper $result;

  #Get BAM file paths from build IDs.  Perform sanity checks
  my $data;
  $data = &getFilePaths_Genome(
    '-wgs_som_var_model_id'     => ($wgs_som_var_build      ? $wgs_som_var_build->model_id : undef), 
    '-exome_som_var_model_id'   => ($exome_som_var_build    ? $exome_som_var_build->model_id : undef), 
    '-rna_seq_normal_model_id'  => ($rna_seq_normal_build   ? $rna_seq_normal_build->model_id : undef), 
    '-rna_seq_tumor_model_id'   => ($rna_seq_tumor_build    ? $rna_seq_tumor_build->model_id : undef)
   );

  #For each mutation get BAM read counts for a tumor/normal pair of BAM files
  foreach my $bam (sort {$a <=> $b} keys %{$data}){
    my $data_type = $data->{$bam}->{data_type};
    my $sample_type = $data->{$bam}->{sample_type};
    my $bam_path = $data->{$bam}->{bam_path};
    my $ref_fasta = $data->{$bam}->{ref_fasta};
    my $snv_count = keys %{$snvs};

    if ($verbose){
	$self->status_message("\n\nSNV count = $snv_count\n$data_type\n$sample_type\n$bam_path\n$ref_fasta\n")
    }
    my $counts = &getBamReadCounts('-snvs'=>$snvs, '-data_type'=>$data_type, '-sample_type'=>$sample_type, '-bam_path'=>$bam_path, '-ref_fasta'=>$ref_fasta, '-verbose'=>$verbose, '-no_fasta_check'=>$no_fasta_check);
    $data->{$bam}->{read_counts} = $counts;
  }


  #Get the FPKM and calculate a percentile value from the RNAseq build dir - do this for tumor and normal if available
  foreach my $bam (sort {$a <=> $b} keys %{$data}){
    my $data_type = $data->{$bam}->{data_type};
    my $sample_type = $data->{$bam}->{sample_type};
    unless ($data_type eq "RNAseq"){
      next();
    }
    my $build_dir = $data->{$bam}->{build_dir};
    my $exp = &getExpressionValues('-snvs'=>$snvs, '-build_dir'=>$build_dir, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-ensembl_map'=>$ensembl_map, '-verbose'=>$verbose);
    $data->{$bam}->{gene_expression} = $exp;
  }


  #Create an output file that is the same as the input file with new columns appended:
  #All of the following are optional
  #1.) WGS Normal Ref Count, WGS Normal Var Count, WGS Normal Var Frequency
  #2.) WGS Tumor Ref Count, WGS Tumor Var Count, WGS Tumor Var Frequency
  #3.) Exome Normal Ref Count, Exome Normal Var Count, Exome Normal Var Frequency
  #4.) Exome Tumor Ref Count, Exome Tumor Var Count, Exome Tumor Var Frequency
  #5.) RNAseq Normal Ref Count, RNAseq Normal Var Count, RNAseq Normal Var Frequency - usually not available
  #6.) RNAseq Normal Gene FPKM, RNAseq Normal Gene Percentile
  #7.) RNAseq Tumor Ref Count, RNAseq Tumor Var Count, RNAseq Tumor Var Frequency
  #8.) RNAseq Tumor Gene FPKM, RNAseq Tumor Gene Percentile

  my %new_snv;
  foreach my $bam (sort {$a <=> $b} keys %{$data}){
    my $data_type = $data->{$bam}->{data_type};
    my $sample_type = $data->{$bam}->{sample_type};
    my $read_counts = $data->{$bam}->{read_counts};

    my $new_header = "\t$data_type"."_"."$sample_type"."_ref_rc\t"."$data_type"."_"."$sample_type"."_var_rc\t"."$data_type"."_"."$sample_type"."_VAF";
    $snv_header .= $new_header;
    foreach my $snv_pos (keys %{$read_counts}){
      my $total_rc = $read_counts->{$snv_pos}->{total_rc};
      my $ref_rc = $read_counts->{$snv_pos}->{ref_rc};
      my $var_rc = $read_counts->{$snv_pos}->{var_rc};
      my $var_allele_frequency = $read_counts->{$snv_pos}->{var_allele_frequency};
      if ($new_snv{$snv_pos}){
        $new_snv{$snv_pos}{read_count_string} .= "\t$ref_rc\t$var_rc\t$var_allele_frequency";
      }else{
        $new_snv{$snv_pos}{read_count_string} = "\t$ref_rc\t$var_rc\t$var_allele_frequency";
      }
    }

    if (defined($data->{$bam}->{gene_expression})){
      my $gene_exp = $data->{$bam}->{gene_expression};
      my $new_header = "\t$data_type"."_"."$sample_type"."_gene_FPKM\t"."$data_type"."_"."$sample_type"."_gene_FPKM_percentile";
      $snv_header .= $new_header;
      foreach my $snv_pos (keys %{$gene_exp}){
        my $fpkm = $gene_exp->{$snv_pos}->{FPKM};
        my $percentile = $gene_exp->{$snv_pos}->{percentile};
        my $rank = $gene_exp->{$snv_pos}->{rank};
        if ($new_snv{$snv_pos}){
          $new_snv{$snv_pos}{read_count_string} .= "\t$fpkm\t$percentile";
	      }else{
	        $new_snv{$snv_pos}{read_count_string} = "\t$fpkm\t$percentile";
	      }
      }
    }
  }

  open (OUT, ">$output_file") || die "\n\nCould not open output file: $output_file\n\n";
  print OUT "$snv_header\n";
  foreach my $snv_pos (sort {$snvs->{$a}->{order} <=> $snvs->{$b}->{order}} keys %{$snvs}){
    my $read_count_string = $new_snv{$snv_pos}{read_count_string};
    print OUT "$snvs->{$snv_pos}->{line}"."$read_count_string\n";  
  }
  close (OUT);

  if ($verbose){print "\n\n";}

  return 1;
}

#########################################################################################################################################
#Import SNVs from the specified file                                                                                                    #
#########################################################################################################################################
sub importPositions{
  my %args = @_;
  my $self = shift;
  my $infile = $args{'-positions_file'};
  my %result;
  my %s;

  my $header = 1;
  my $header_line;
  my %columns;
  my $order = 0;
  open (SNV, "$infile") || die "\n\nCould not open input SNV file: $infile\n\n";
  while(<SNV>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header == 1){
      my $p = 0;
      foreach my $head (@line){
        $columns{$head}{position} = $p;
        $p++;
      }
      $header = 0;
      $header_line = $_;
      #Make sure all neccessary columns are defined
      unless (defined($columns{'coord'}) && defined($columns{'mapped_gene_name'}) && defined($columns{'ref_base'}) && defined($columns{'var_base'}) && defined($columns{'ensembl_gene_id'})){
        die $self->error_message("\n\nRequired column missing from file: $infile (need: coord, mapped_gene_name, ref_base, var_base, ensembl_gene_id)");
      }
      next();
    }
    $order++;
    my $coord = $line[$columns{'coord'}{position}];
    $s{$coord}{order} = $order;
    $s{$coord}{mapped_gene_name} = $line[$columns{'mapped_gene_name'}{position}];
    $s{$coord}{ensembl_gene_id} = $line[$columns{'ensembl_gene_id'}{position}];
    $s{$coord}{ref_base} = $line[$columns{'ref_base'}{position}];
    $s{$coord}{var_base} = $line[$columns{'var_base'}{position}];
    $s{$coord}{line} = $_;

    if ($coord =~ /(\S+)\:(\d+)\-(\d+)/){
      $s{$coord}{chr} = $1;
      $s{$coord}{start} = $2;
      $s{$coord}{end} = $3;
    }else{
      die $self->error_message("\n\nCoord: $coord not understood\n\n");
    }

  }
  close(SNV);
  $result{'snvs'} = \%s;
  $result{'header'} = $header_line;
  return(\%result);
}


#########################################################################################################################################
#getFilePaths_Genome - Get file paths from model IDs                                                                                    #
#########################################################################################################################################
sub getFilePaths_Genome{
  my %args = @_;
  my $self = shift;
  my $wgs_som_var_model_id = $args{'-wgs_som_var_model_id'};
  my $exome_som_var_model_id = $args{'-exome_som_var_model_id'};
  my $rna_seq_normal_model_id = $args{'-rna_seq_normal_model_id'};
  my $rna_seq_tumor_model_id = $args{'-rna_seq_tumor_model_id'};

  my %d;

  my $b = 0;
  #WGS tumor normal BAMs
  if ($wgs_som_var_model_id){
    my $wgs_som_var_model = Genome::Model->get($wgs_som_var_model_id);
    if ($wgs_som_var_model){
      my $wgs_som_var_build = $wgs_som_var_model->last_succeeded_build;
      if ($wgs_som_var_build){

        #... /genome/lib/perl/Genome/Model/Build/SomaticVariation.pm
        my $reference_build = $wgs_som_var_build->reference_sequence_build;
        my $reference_fasta_path = $reference_build->full_consensus_path('fa');
        my $reference_display_name = $reference_build->__display_name__;
        my $build_dir = $wgs_som_var_build->data_directory ."/";
        $b++;
        $d{$b}{build_dir} = $build_dir;
        $d{$b}{data_type} = "WGS";
        $d{$b}{sample_type} = "Normal";
        $d{$b}{bam_path} = $wgs_som_var_build->normal_bam;
        $d{$b}{ref_fasta} = $reference_fasta_path;
        $d{$b}{ref_name} = $reference_display_name;
        $b++;
        $d{$b}{build_dir} = $build_dir;
        $d{$b}{data_type} = "WGS";
        $d{$b}{sample_type} = "Tumor";
        $d{$b}{bam_path} = $wgs_som_var_build->tumor_bam;
        $d{$b}{ref_fasta} = $reference_fasta_path;
        $d{$b}{ref_name} = $reference_display_name;
      }else{
        die $self->error_message("\n\nA WGS model ID was specified, but a successful build could not be found!\n\n");
      }
    }else{
      die $self->error_message("\n\nA WGS model ID was specified, but it could not be found!\n\n");
    }
  }

  #Exome tumor normal BAMs
  if ($exome_som_var_model_id){
    my $exome_som_var_model = Genome::Model->get($exome_som_var_model_id);
    if ($exome_som_var_model){
      my $exome_som_var_build = $exome_som_var_model->last_succeeded_build;
      if ($exome_som_var_build){
        #... /genome/lib/perl/Genome/Model/Build/SomaticVariation.pm
        my $reference_build = $exome_som_var_build->reference_sequence_build;
        my $reference_fasta_path = $reference_build->full_consensus_path('fa');
        my $reference_display_name = $reference_build->__display_name__;
        my $build_dir = $exome_som_var_build->data_directory ."/";
        $b++;
        $d{$b}{build_dir} = $build_dir;
        $d{$b}{data_type} = "Exome";
        $d{$b}{sample_type} = "Normal";
        $d{$b}{bam_path} = $exome_som_var_build->normal_bam;
        $d{$b}{ref_fasta} = $reference_fasta_path;
        $d{$b}{ref_name} = $reference_display_name;
        $b++;
        $d{$b}{build_dir} = $build_dir;
        $d{$b}{data_type} = "Exome";
        $d{$b}{sample_type} = "Tumor";
        $d{$b}{bam_path} = $exome_som_var_build->tumor_bam;
        $d{$b}{ref_fasta} = $reference_fasta_path;
        $d{$b}{ref_name} = $reference_display_name;
      }else{
        die $self->error_message("\n\nA Exome model ID was specified, but a successful build could not be found!\n\n");
      }
    }else{
      die $self->error_message("\n\nA Exome model ID was specified, but it could not be found!\n\n");
    }
  }

  #RNAseq normal BAM
  if ($rna_seq_normal_model_id){
    my $rna_seq_model = Genome::Model->get($rna_seq_normal_model_id);
      if ($rna_seq_model){
      my $rna_seq_build = $rna_seq_model->last_succeeded_build;
      if ($rna_seq_build){
        my $reference_build = $rna_seq_model->reference_sequence_build;
        my $reference_fasta_path = $reference_build->full_consensus_path('fa');
        my $reference_display_name = $reference_build->__display_name__;
        my $build_dir = $rna_seq_build->data_directory ."/";
        $b++;
        $d{$b}{build_dir} = $build_dir;
        $d{$b}{data_type} = "RNAseq";
        $d{$b}{sample_type} = "Normal";
        my $alignment_result = $rna_seq_build->alignment_result;
        $d{$b}{bam_path} = $alignment_result->bam_file;
        $d{$b}{ref_fasta} = $reference_fasta_path;
        $d{$b}{ref_name} = $reference_display_name;
      }else{
        die $self->error_message("\n\nAn RNA-seq model ID was specified, but a successful build could not be found!\n\n");
      }
    }else{
      die $self->error_message("\n\nAn RNA-seq model ID was specified, but it could not be found!\n\n");
    }
  }

  #RNAseq tumor BAM
  if ($rna_seq_tumor_model_id){
    my $rna_seq_model = Genome::Model->get($rna_seq_tumor_model_id);
      if ($rna_seq_model){
      my $rna_seq_build = $rna_seq_model->last_succeeded_build;
      if ($rna_seq_build){
        my $reference_build = $rna_seq_model->reference_sequence_build;
        my $reference_fasta_path = $reference_build->full_consensus_path('fa');
        my $reference_display_name = $reference_build->__display_name__;
        my $build_dir = $rna_seq_build->data_directory ."/";
        $b++;
        $d{$b}{build_dir} = $build_dir;
        $d{$b}{data_type} = "RNAseq";
        $d{$b}{sample_type} = "Tumor";
        my $alignment_result = $rna_seq_build->alignment_result;
        $d{$b}{bam_path} = $alignment_result->bam_file;
        $d{$b}{ref_fasta} = $reference_fasta_path;
        $d{$b}{ref_name} = $reference_display_name;
      }else{
        die $self->error_message("\n\nAn RNA-seq model ID was specified, but a successful build could not be found!\n\n");
      }
    }else{
      die $self->error_message("\n\nAn RNA-seq model ID was specified, but it could not be found!\n\n");
    }
  }

  #Make sure the same reference build was used to create all BAM files!
  my $test_ref_name = $d{1}{ref_name};
  foreach my $b (keys %d){
    my $ref_name = $d{$b}{ref_name};
    unless ($ref_name eq $test_ref_name){
      print Dumper %d;
      die $self->error_message("\n\nOne or more of the reference build names used to generate BAMs did not match\n\n");
    }
  }

  return(\%d)
}


#########################################################################################################################################
#getBamReadCounts                                                                                                                       #
#########################################################################################################################################
sub getBamReadCounts{
  my %args = @_;
  my $self = shift;
  my $snvs = $args{'-snvs'};
  my $data_type = $args{'-data_type'};
  my $sample_type = $args{'-sample_type'};
  my $bam_path = $args{'-bam_path'};
  my $ref_fasta = $args{'-ref_fasta'};
  my $verbose = $args{'-verbose'};
  my $no_fasta_check = $args{'-no_fasta_check'};
  my %c;

  #Code reference needed for Bio::DB::Bam
  my $callback = sub {
    my ($tid,$pos,$pileups,$callback_data) = @_;
    my $data = $callback_data->[0];
    my $read_counts = $callback_data->[1];
    my $fai = $callback_data->[2];

    if ( ($pos == ($data->{start} - 1) ) ) {
      #print STDERR 'PILEUP:'. $data->{chr} ."\t". $tid ."\t". $pos ."\t". $data->{start} ."\t". $data->{stop}."\n";
      unless ($no_fasta_check){
        my $ref_base = $fai->fetch($data->{chr} .':'. $data->{start} .'-'. $data->{stop});
        unless ($data->{reference} eq $ref_base) {
          #print RED, "\n\nReference base " . $ref_base .' does not match expected '. $data->{reference} .' at postion '. $pos .' for chr '. $data->{chr} . '(tid = '. $tid . ')' . "\n$bam_path", RESET;
          die $self->error_message("\n\nReference base " . $ref_base .' does not match expected '. $data->{reference} .' at postion '. $pos .' for chr '. $data->{chr} . '(tid = '. $tid . ')' . "\n$bam_path");
        }
      }
      for my $pileup ( @{$pileups} ) {
        my $alignment = $pileup->alignment;
        #Skip indels or skip regions
        next if $pileup->indel or $pileup->is_refskip;
        my $qbase  = substr($alignment->qseq,$pileup->qpos,1);
        next if $qbase =~ /[nN]/;
        $read_counts->{$qbase}++;
      }
    }
  };

  #Get Bio:DB:Bam objects for this BAM file and reference fasta file
  my $bam = Bio::DB::Bam->open($bam_path);
  my $header = $bam->header;
  my $index = Bio::DB::Bam->index($bam_path);
  my $fai = Bio::DB::Sam::Fai->load($ref_fasta);

  #my $name_arrayref = $header->target_name;
  #my %chr_tid;
  #for (my $i = 0; $i < scalar@{$name_arrayref}; $i++){
  #  my $seq_id = $name_arrayref->[$i];
  #  $chr_tid{$seq_id}=$i;
  #}
  #print Dumper %chr_tid;

  foreach my $snv_pos (keys %{$snvs}){
    my %data;
    my $data = \%data;
    $data->{chr} = $snvs->{$snv_pos}->{chr};
    $data->{start} = $snvs->{$snv_pos}->{start};
    $data->{stop} = $snvs->{$snv_pos}->{end};
    $data->{reference} = $snvs->{$snv_pos}->{ref_base};
    $data->{variant} = $snvs->{$snv_pos}->{var_base};
    my $seq_id = $data->{chr} .':'. $data->{start} .'-'. $data->{stop};
    my ($tid,$start,$end) = $header->parse_region($seq_id);
    #$tid = $chr_tid{$data->{chr}};
    #print "\n\nseq_id: $seq_id\ttid: $tid\tstart: $start\tend: $end";

    my %read_counts;
    if ($verbose){print "\n\n$sample_type\t$data_type\t$snv_pos\ttid: $tid\tstart: $start\tend: $end\tref_base: $data->{reference}\tvar_base: $data->{variant}";}
    $index->pileup($bam,$tid,$start,$end,$callback,[$data,\%read_counts, $fai]);

    $data->{A} = $read_counts{A} || 0;
    $data->{T} = $read_counts{T} || 0;
    $data->{C} = $read_counts{C} || 0;
    $data->{G} = $read_counts{G} || 0;

    if ($verbose){print "\n\tA: $data->{A}\tT: $data->{T}\tC: $data->{C}\tG: $data->{G}";}

    #Store ref read count, var read count, var allele frequency
    my $total_rc =  $data->{A} + $data->{T} + $data->{C} + $data->{G};
    my $ref_rc = $data->{$snvs->{$snv_pos}->{ref_base}};
    my $var_rc = $data->{$snvs->{$snv_pos}->{var_base}};
    my $var_allele_frequency = 0;
    if ($total_rc){
      $var_allele_frequency = sprintf ("%.3f", (($var_rc / $total_rc)*100));
    }
    if ($verbose){print "\n\t\tTotalCount: $total_rc\tRefCount: $ref_rc\tVarCount: $var_rc\tVAF: $var_allele_frequency%";}

    $c{$snv_pos}{total_rc} = $total_rc;
    $c{$snv_pos}{ref_rc} = $ref_rc;
    $c{$snv_pos}{var_rc} = $var_rc;
    $c{$snv_pos}{var_allele_frequency} = $var_allele_frequency;

  }

  return(\%c);
}


#########################################################################################################################################
#getExpressionValues                                                                                                                    #
#########################################################################################################################################
sub getExpressionValues{
  my %args = @_;
  my $class = __PACKAGE__;
  my $snvs = $args{'-snvs'};
  my $build_dir = $args{'-build_dir'};
  my $verbose = $args{'-verbose'};
  my $entrez_ensembl_data = $args{'-entrez_ensembl_data'};
  my $ensembl_map = $args{'-ensembl_map'};

  if ($verbose){ $class->status_message("\nGetting expression data from: $build_dir");}

  my %e;

  #Import FPKM values from the gene-level expression file created by merging the isoforms of each gene
  my $isoforms_infile = "$build_dir"."expression/isoforms.fpkm_tracking";
  my $merged_fpkm = &mergeIsoformsFile('-infile'=>$isoforms_infile, '-entrez_ensembl_data'=>$entrez_ensembl_data, '-ensembl_map'=>$ensembl_map, '-verbose'=>$verbose);
  
  #Calculate the ranks and percentiles for all genes
  my $rank = 0;
  my $gene_count = keys %{$merged_fpkm};
  foreach my $ensembl_gene_id (sort {$merged_fpkm->{$b}->{FPKM} <=> $merged_fpkm->{$a}->{FPKM}} keys %{$merged_fpkm}){
    $rank++;
    $merged_fpkm->{$ensembl_gene_id}->{rank} = $rank;
    $merged_fpkm->{$ensembl_gene_id}->{percentile} = sprintf("%.3f", (($rank/$gene_count)*100));
  }

  foreach my $snv_pos (keys %{$snvs}){
    my $mapped_gene_name = $snvs->{$snv_pos}->{mapped_gene_name};
    my $ensembl_gene_id = $snvs->{$snv_pos}->{ensembl_gene_id};
    if ($merged_fpkm->{$ensembl_gene_id}){
      $e{$snv_pos}{FPKM} = $merged_fpkm->{$ensembl_gene_id}->{FPKM};
      $e{$snv_pos}{rank} = $merged_fpkm->{$ensembl_gene_id}->{rank};
      $e{$snv_pos}{percentile} = $merged_fpkm->{$ensembl_gene_id}->{percentile};
    }else{
      #Unmappable
      $e{$snv_pos}{FPKM} = "NA";
      $e{$snv_pos}{rank} = "NA";
      $e{$snv_pos}{percentile} = "NA";
    }
  }

  return(\%e);
}


