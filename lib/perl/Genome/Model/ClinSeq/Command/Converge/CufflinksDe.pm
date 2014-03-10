package Genome::Model::ClinSeq::Command::Converge::CufflinksDe;
use strict;
use warnings;
use Genome;
use Data::Dumper;

class Genome::Model::ClinSeq::Command::Converge::CufflinksDe {
    is => 'Genome::Model::ClinSeq::Command::Converge::Base',
    doc => 'converge cufflinks differential expression results from mutiple clinseq builds',
    has_input => [
      fc_cutoff => {
          is => 'Number',
          is_optional => 1,
          doc => 'Apply a log2 fold-change cutoff for a DE result to be considered',
      },
    ],
};

sub help_synopsis {
  return <<EOS

genome model clin-seq converge cufflinks-de --builds='id in [133577030,133611960]' --outdir=/tmp/converge_de/

genome model clin-seq converge cufflinks-de --builds='model.model_groups.id=50714,is_last_complete=1' --outdir=/tmp/converge_de/

genome model clin-seq converge cufflinks-de --builds='model_groups.id=50714,is_last_complete=1' --outdir=/tmp/converge_de/

EOS
}

sub execute {
  my $self = shift;
  my @builds = $self->builds;
  my $outdir = $self->outdir;
  $outdir .= "/" unless ($outdir =~ /\/$/);

  #Produce output files of the general format (and do this for up, down, up-or-down, complete-matrix):
  #ensembl_id	*gene_annotations* SubjectCount	SubjectList	*de values, one column per subject*
  
  #Target out files:
  #de.genes.up.tsv, de.genes.down.tsv, de.genes.updown.tsv, de.genes.matrix.tsv
  #de.transcripts.up.tsv, de.transcripts.down.tsv, de.transcripts.updown.tsv, de.transcripts.matrix.tsv

  #Get human readable names hash, keyed on build id
  my $labels = $self->resolve_clinseq_subject_labels;

  my @feature_types = qw /genes transcripts/;
  my @event_types = qw /up down de/;
  my @biotypes = ('', 'coding.');
  
  foreach my $feature_type (@feature_types){

    #Get the master DE matrix for this feature type (genes/transcripts)
    my $target = "rnaseq/cufflinks_differential_expression/$feature_type/case_vs_control.tsv";
    my $files = $self->get_clinseq_files('-target'=>$target);

    my $de_matrix_outfile = $outdir . $feature_type . ".case_vs_control.matrix.tsv";
    my $de_matrix = $self->get_de_matrix('-outfile'=>$de_matrix_outfile, '-files'=>$files, '-labels'=>$labels);


    foreach my $event_type (@event_types){
      foreach my $biotype (@biotypes){

        #Produce output for each sub-category but get the complete list of DE values for all samples from the master DE matrix

        #Get DE up/down/de genes/transcripts files hash, keyed on build id
        my $target = "rnaseq/cufflinks_differential_expression/$feature_type/case_vs_control.$biotype"."hq.$event_type.tsv";
        my $files = $self->get_clinseq_files('-target'=>$target);

        #Create converged file for this set of DE files
        my $outfile = $outdir . $feature_type . ".case_vs_control.$biotype" . "hq.$event_type.tsv";
        $self->converge_de_data('-outfile'=>$outfile, '-files'=>$files, '-labels'=>$labels, '-de_matrix'=>$de_matrix);
      }
    }
  }

  return 1;
}


sub get_de_matrix{
  my $self = shift;
  my %args = @_;
  my $outfile = $args{'-outfile'};
  my $files = $args{'-files'};
  my $labels = $args{'-labels'};
  my %data;
  
  $self->debug_message("Creating de matrix: $outfile");

  my $gene_col;
  foreach my $build_id (keys %{$files}){
    my $label = $labels->{$build_id}->{name};
    my $path = $files->{$build_id}->{path};

    open (DATA, $path) || die $self->error_message("Could not open file: $path");
    my $header = 1;
    my %columns;
    while(<DATA>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header){
        my $p = 0;
        foreach my $column (@line){
          $columns{$column}{p} = $p;
          $p++;
        }
        $header = 0;
        unless (defined($columns{'tracking_id'}) && defined($columns{'mapped_gene_name'}) && (defined($columns{'gene_id'}) || defined($columns{'ensg_name'})) && defined($columns{'biotype'}) && defined($columns{'case_vs_control_log2_de'})){
          die $self->error_message("Required column missing from $path (tracking_id mapped_gene_name gene_id ensg_name biotype case_vs_control_log2_de)");
        }
        if (defined($columns{'gene_id'})){
          $gene_col = 'gene_id';
        }elsif(defined($columns{'ensg_name'})){
          $gene_col = 'ensg_name';
        }
        next;
      }
      my $tracking_id = $line[$columns{'tracking_id'}{p}];
      my $mapped_gene_name = $line[$columns{'mapped_gene_name'}{p}];
      my $gene_id = $line[$columns{$gene_col}{p}];
      my $biotype = $line[$columns{'biotype'}{p}];
      my $case_vs_control_log2_de = $line[$columns{'case_vs_control_log2_de'}{p}];

      if ($data{$tracking_id}){
        my $subjects = $data{$tracking_id}{subjects};
        $subjects->{$label}->{de} = $case_vs_control_log2_de;
      }else{
        $data{$tracking_id}{mapped_gene_name} = $mapped_gene_name;
        $data{$tracking_id}{gene_id} = $gene_id;
        $data{$tracking_id}{biotype} = $biotype;
        my %subjects;
        $subjects{$label}{de} = $case_vs_control_log2_de;
        $data{$tracking_id}{subjects} = \%subjects;
      }
    }
    close(DATA);
  }

  my @labels_list;
  foreach my $build_id (sort {$labels->{$a}->{name} cmp $labels->{$b}->{name}} keys %{$labels}){
    push(@labels_list, $labels->{$build_id}->{name});
  }
  my $labels_list_string = join("\t", @labels_list);

  open (OUT, ">$outfile") || die $self->error_message("Could not open $outfile for writing");
  my $header_line = "tracking_id\tmapped_gene_name\t$gene_col\tbiotype\t$labels_list_string\n";
  print OUT $header_line;
  foreach my $tracking_id (keys %data){
    my $subjects = $data{$tracking_id}{subjects};
    my @data;
    foreach my $label (@labels_list){
      push(@data, $subjects->{$label}->{de});
    }
    my $data_string = join("\t", @data);
    print OUT "$tracking_id\t$data{$tracking_id}{mapped_gene_name}\t$data{$tracking_id}{gene_id}\t$data{$tracking_id}{biotype}\t$data_string\n";
  }
  close (OUT);

  return (\%data);
}

sub converge_de_data{
  my $self = shift;
  my %args = @_;
  my $outfile = $args{'-outfile'};
  my $files = $args{'-files'};
  my $labels = $args{'-labels'};
  my $de_matrix = $args{'-de_matrix'};
  my %data;
  
  $self->debug_message("Creating matrix: $outfile");

  #Example input
  # tracking_id mapped_gene_name  gene_id biotype locus case_fpkm case_fpkm_conf_hi case_fpkm_conf_lo case_fpkm_status  control_fpkm  control_fpkm_conf_hi  control_fpkm_conf_lo  control_fpkm_status case_vs_control_log2_de case_vs_control_de_lq case_vs_control_de_hq

  #Desired output
  #ensembl_id	*gene_annotations* SubjectCount	SubjectList	*de values, one column per subject*
    
  my $gene_col;
  foreach my $build_id (keys %{$files}){
    my $label = $labels->{$build_id}->{name};
    my $path = $files->{$build_id}->{path};

    open (DATA, $path) || die $self->error_message("Could not open file: $path");
    my $header = 1;
    my %columns;
    while(<DATA>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header){
        my $p = 0;
        foreach my $column (@line){
          $columns{$column}{p} = $p;
          $p++;
        }
        $header = 0;
        unless (defined($columns{'tracking_id'}) && defined($columns{'mapped_gene_name'}) && (defined($columns{'gene_id'}) || defined($columns{'ensg_name'})) && defined($columns{'biotype'}) && defined($columns{'case_vs_control_log2_de'})){
          die $self->error_message("Required column missing from $path (tracking_id mapped_gene_name gene_id ensg_name biotype case_vs_control_log2_de)");
        }
        if (defined($columns{'gene_id'})){
          $gene_col = 'gene_id';
        }elsif(defined($columns{'ensg_name'})){
          $gene_col = 'ensg_name';
        }
        next;
      }
      my $tracking_id = $line[$columns{'tracking_id'}{p}];
      my $mapped_gene_name = $line[$columns{'mapped_gene_name'}{p}];
      my $gene_id = $line[$columns{$gene_col}{p}];
      my $biotype = $line[$columns{'biotype'}{p}];
      my $case_vs_control_log2_de = $line[$columns{'case_vs_control_log2_de'}{p}];

      #Apply fold change cutoff is specified by the user
      if ($self->fc_cutoff){
        next unless (abs($case_vs_control_log2_de) >= $self->fc_cutoff);
      }

      if ($data{$tracking_id}){
        my $subjects = $data{$tracking_id}{subjects};
        $subjects->{$label} = 1;
      }else{
        $data{$tracking_id}{mapped_gene_name} = $mapped_gene_name;
        $data{$tracking_id}{gene_id} = $gene_id;
        $data{$tracking_id}{biotype} = $biotype;
        my %subjects;
        $subjects{$label} = 1;
        $data{$tracking_id}{subjects} = \%subjects;
      }
    }
    close(DATA);
  }

  my @labels_list;
  foreach my $build_id (sort {$labels->{$a}->{name} cmp $labels->{$b}->{name}} keys %{$labels}){
    push(@labels_list, $labels->{$build_id}->{name});
  }
  my $labels_list_string = join("\t", @labels_list);

  open (OUT, ">$outfile") || die $self->error_message("Could not open $outfile for writing");
  my $header_line = "tracking_id\tmapped_gene_name\t$gene_col\tbiotype\tde_subject_count\tde_subject_list\t$labels_list_string\n";
  print OUT $header_line;
  foreach my $tracking_id (keys %data){
    my $subjects = $data{$tracking_id}{subjects};
    my @subject_list = keys %{$subjects};
    my @subject_list_sorted = sort {$a cmp $b} @subject_list;
    my $de_subject_list_string = join(",", @subject_list_sorted);
    my $de_subject_count = scalar(@subject_list);

    my @data;
    my $complete_subjects = $de_matrix->{$tracking_id}->{subjects};
    foreach my $label (@labels_list){
      push(@data, $complete_subjects->{$label}->{de});
    }
    my $data_string = join("\t", @data);
    print OUT "$tracking_id\t$data{$tracking_id}{mapped_gene_name}\t$data{$tracking_id}{gene_id}\t$data{$tracking_id}{biotype}\t$de_subject_count\t$de_subject_list_string\t$data_string\n";

  }
  close (OUT);


  return;
}

1;
