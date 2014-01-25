package Genome::Model::ClinSeq::Command::Converge::AllEvents;
use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq::Util qw(:all);

class Genome::Model::ClinSeq::Command::Converge::AllEvents {
    is => 'Genome::Model::ClinSeq::Command::Converge::Base',
    has_input => [
        outdir => { 
               is => 'FilesystemPath',
               doc => 'Directory where output files will be written',
        },
    ],
    has_optional_input => [
        snv_label => {
                is => 'Text',
                doc => 'Consider SNVs and label affected samples with this string.  e.g. S, SNV, etc.',
        },
        indel_label => {
                is => 'Text',
                doc => 'Consider Indels and label affected samples with this string.  e.g. I, Indel, etc.',
        },
        cnv_gain_label => {
                is => 'Text',
                doc => 'Consider DNA CNV gains (amplifications) and label affected samples with this string.  e.g. A, CnGain, etc.',
        },
        cnv_loss_label => {
                is => 'Text',
                doc => 'Consider DNA CNV losses (deletions) and label affected samples with this string.  e.g. D, CnLoss, etc.',
        },
        de_up_label => {
                is => 'Text',
                doc => 'Consider RNA DE gains and label affected samples with this string.  e.g. G, RnaUp, etc.',
        },
        de_down_label => {
                is => 'Text',
                doc => 'Consider RNA DE losses and label affected samples with this string.  e.g. L, RnaUp, etc.',
        },
        sv_fusion_label => {
                is => 'Text',
                doc => 'Consider predicted SV fusion gene pairs and label affected samples with this string.  e.g. T, SV, etc.',
        },
        tophat_outlier_label => {
                is => 'Text',
                doc => 'Consider top N% tophat junction expressed genes and label affected samples with this string. e.g. J, JO, etc.',
        },
        cufflinks_outlier_label => {
                is => 'Text',
                doc => 'Consider top N% cufflinks expressed genes and label affected samples with this string.  e.g. C, CO, etc.',
        },
        target_gene_list => {
                is => 'FilesystemPath',
                doc => 'Limit consideration to genes in this list (a TSV file containing Ensembl Gene IDs of the form ENSGXXXXXXXXX in the first column)',
        },
        ignore_gene_list => {
                is => 'FilesystemPath',
                doc => 'Regardless of anything else do not allow any genes in this list (Ensembl Gene IDs)',
        },
        subject_labels_file => {
                is => 'FilesystemPath',
                doc => 'Use a custom subjects_legend file.  First run without this option, then copy the legend created, modify and then specify with this option. (use to change order, sample names, etc.)',
        },
        event_labels_file => {
                is => 'FilesystemPath',
                doc => 'Use a custom events_legend file.  First run without this option, then copy the legend created, modify and then specify with this option. (use to change event labels, colors, etc.)',
        },
        max_genes => {
                is => 'Number',
                default => 50,
                doc => 'Maximum number of genes that will be displayed (rows in the heatmap)',
        },
    ],   
    doc => 'converge SNV, InDel, CNV, Exp, and DE results from mutiple clinseq builds into a single table',
};

sub help_synopsis {
  return <<EOS

genome model clin-seq converge all-events --builds='id in [133577030,133611960]' --outdir=/tmp/converge_all_events/

genome model clin-seq converge all-events --builds='model.model_groups.id=50714,is_last_complete=1' --outdir=/tmp/converge_all_events/

genome model clin-seq converge all-events --builds='model_groups.id=50714,is_last_complete=1' --outdir=/tmp/converge_all_events/

Example event selections:

--snv-label=S --indel-label=I --cnv-gain-label=A --cnv-loss-label=D --de-up-label=G --de-down-label=L --sv-fusion-label=T --tophat-outlier-label=J --cufflinks-outlier-label=C

EOS
}

sub help_detail {
  return <<EOS

A master converge script to help create a heatmap of genes, affected by each event type in each sample in a clinseq model group

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

  if ($self->ignore_gene_list){
    unless (-e $self->ignore_gene_list) {
      push @errors, UR::Object::Tag->create(
                                            type => 'error',
                                            properties => ['ignore_gene_list'],
                                            desc => "File: " . $self->ignore_gene_list . " not found",
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

  #Print out a table of subject names for later use
  $self->print_subject_table('-subject_labels'=>$subject_labels);

  #If the user specified an alternative subject table file, override the $subject_labels object here
  if ($self->subject_labels_file){
    $subject_labels = $self->override_clinseq_subject_labels('-old_subject_labels'=>$subject_labels);
  }

  #Check event labels to make sure they are all unique (no duplicates)
  my $event_labels = $self->check_event_labels;

  #If the user specied an alternative event label file, overwrite the events legend file that will be fed into R
  if ($self->event_labels_file){
    $self->override_event_labels;
  }

  #Get the reference sequence build common to all underlying models of all clinseq builds
  my $reference_build = $self->resolve_clinseq_reference_build;
  my $reference_build_id = $reference_build->id;

  #Get the reference annotation build common to all underlying models of all clinseq builds
  my $annotation_build = $self->resolve_clinseq_annotation_build;
  my $annotation_build_name = $annotation_build->name;
  my $annotation_data_dir = $annotation_build->data_directory;
  my $transcript_info_path = $annotation_data_dir . "/annotation_data/rna_annotation/$reference_build_id-transcript_info.tsv";
  my $gtf_path = $annotation_build->annotation_file('gtf',$reference_build_id);
  $self->status_message("Getting transcript to gene and gene name mappings from annotation build: $annotation_build_name");
  unless (defined($gtf_path)) {
    $self->error_message("'There is no annotation GTF file defined for annotation_reference_transcripts build: ". $annotation_build->__display_name__);
    die $self->error_message;
  }
  unless (-e $transcript_info_path) {
    $self->error_message("'There is no transcript info file for annotation_reference_transcripts build: ". $annotation_build->__display_name__);
    die $self->error_message;
  }
  $self->status_message("\t$transcript_info_path");

  #Get Ensembl gene ID to name mappings from the annotation build
  my $tmap = $self->loadEnsemblMap('-gtf_path'=>$gtf_path, '-transcript_info_path'=>$transcript_info_path);

  #Create a new ensembl map that is keyed on gene ID
  my %gmap;
  foreach my $tid (keys %{$tmap}){
    my $gid = $tmap->{$tid}->{ensg_id};
    $gmap{$gid}{ensg_name} = $tmap->{$tid}->{ensg_name};
  }
  my $gmap = \%gmap;

  #Gather all events of all types and store in a single data structure.
  #gene -> sample -> event type
  #Do not allow genes that can not be found by name/id in the ensembl_map
  $self->status_message("\nBegin aggregating events of all types\n");
  my %targets;

  foreach my $clinseq_build (@builds){
    my $clinseq_build_id = $clinseq_build->id;
    my $label = $subject_labels->{$clinseq_build_id}->{name_abr};
    $self->status_message("Gathering events for clinseq build: $clinseq_build_id (" . $label . ")");

    #SNVs
    if ($self->snv_label){
      my $file;
      if ($clinseq_build->wgs_build and $clinseq_build->exome_build){
        $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.tsv');
      }elsif ($clinseq_build->wgs_build){
        $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'snv/wgs/snvs.hq.tier1.v1.annotated.compact.tsv');
      }elsif ($clinseq_build->exome_build){
        $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'snv/exome/snvs.hq.tier1.v1.annotated.compact.tsv');
      }
      $self->add_events('-genes'=>\%targets, '-gmap'=>$gmap, '-file'=>$file, '-column_name'=>'ensembl_gene_id', '-sample_label'=>$label, '-event_label'=>$self->snv_label);
    }

    #INDELS
    if ($self->indel_label){
      my $file;
      if ($clinseq_build->wgs_build and $clinseq_build->exome_build){
        $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'indel/wgs_exome/indels.hq.tier1.v1.annotated.compact.tsv');
      }elsif ($clinseq_build->wgs_build){
        $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'indel/wgs/indels.hq.tier1.v1.annotated.compact.tsv');
      }elsif ($clinseq_build->exome_build){
        $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'indel/exome/indels.hq.tier1.v1.annotated.compact.tsv');
      }
      $self->add_events('-genes'=>\%targets, '-gmap'=>$gmap, '-file'=>$file, '-column_name'=>'ensembl_gene_id', '-sample_label'=>$label, '-event_label'=>$self->indel_label);
    }

    #CNV-GAINS
    if ($self->cnv_gain_label){
      if ($clinseq_build->wgs_build){
        my $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'cnv/cnview/cnv.All_genes.amp.tsv');
        $self->add_events('-genes'=>\%targets, '-gmap'=>$gmap, '-file'=>$file, '-column_name'=>'gene_id', '-sample_label'=>$label, '-event_label'=>$self->cnv_gain_label);
      }
    }

    #CNV-LOSSES
    if ($self->cnv_loss_label){
      if ($clinseq_build->wgs_build){
        my $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'cnv/cnview/cnv.All_genes.del.tsv');
        $self->add_events('-genes'=>\%targets, '-gmap'=>$gmap, '-file'=>$file, '-column_name'=>'gene_id', '-sample_label'=>$label, '-event_label'=>$self->cnv_loss_label);
      }
    }

    #DE-UP
    if ($self->de_up_label){
      if ($clinseq_build->normal_rnaseq_build and $clinseq_build->tumor_rnaseq_build){
        my $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'rnaseq/cufflinks_differential_expression/genes/case_vs_control.coding.hq.up.tsv');
        $self->add_events('-genes'=>\%targets, '-gmap'=>$gmap, '-file'=>$file, '-column_name'=>'tracking_id', '-sample_label'=>$label, '-event_label'=>$self->de_up_label);
      }
    }

    #DE-DOWN
    if ($self->de_down_label){
      if ($clinseq_build->normal_rnaseq_build and $clinseq_build->tumor_rnaseq_build){
        my $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'rnaseq/cufflinks_differential_expression/genes/case_vs_control.coding.hq.down.tsv');
        $self->add_events('-genes'=>\%targets, '-gmap'=>$gmap, '-file'=>$file, '-column_name'=>'tracking_id', '-sample_label'=>$label, '-event_label'=>$self->de_down_label);
      }
    }

    #SV-FUSION
    #TODO: SVs are not based on our standard annotations yet :(.  Will have to deal with these as a special case
    if ($self->sv_fusion_label){
      if ($clinseq_build->wgs_build){
        $self->warning_message("SV annotations are not consistent with other event types!");
        my $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'sv/CandidateSvCodingFusions.tsv');
        $self->add_sv_events('-genes'=>\%targets, '-gmap'=>$gmap, '-file'=>$file, '-sample_label'=>$label, '-event_label'=>$self->sv_fusion_label);
      }
    }

    #TOPHAT JUNCTION GENE OUTLIER
    if ($self->tophat_outlier_label){
      if ($clinseq_build->tumor_rnaseq_build){
        my $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'rnaseq/tumor/tophat_junctions_absolute/Junction.GeneExpression.topnpercent.tsv');
        $self->add_events('-genes'=>\%targets, '-gmap'=>$gmap, '-file'=>$file, '-column_name'=>'ensg_id', '-sample_label'=>$label, '-event_label'=>$self->tophat_outlier_label);
      }
    }

    #CUFFLINKS GENE OUTLIER
    if ($self->cufflinks_outlier_label){
      if ($clinseq_build->tumor_rnaseq_build){
        my $file = $self->get_clinseq_file('-clinseq_build'=>$clinseq_build, '-target'=>'rnaseq/tumor/cufflinks_expression_absolute/isoforms_merged/isoforms.merged.fpkm.expsort.top1percent.tsv');
        $self->add_events('-genes'=>\%targets, '-gmap'=>$gmap, '-file'=>$file, '-column_name'=>'tracking_id', '-sample_label'=>$label, '-event_label'=>$self->cufflinks_outlier_label);
      }
    }
  }

  #The final table can optionally be limited to only those genes in a pre-approved list
  #If the user supplied a target-gene-list, make sure all these genes are in the target list (even if no events were found)
  #Also remove any genes that are not defined in the target-gene-list
  #This way an entry will be produced for every gene in this list even if no events were found for some of the genes
  if ($self->target_gene_list){
    $self->target_specific_genes('-genes'=>\%targets, '-gmap'=>$gmap);
  }

  #The final table can optionally be filtered to exclude a genes in a user supplied list
  #TODO: If the user supplied a ignore-gene-list, remove the genes in this from consideration (even if some events were found)
  if ($self->ignore_gene_list){
    $self->remove_ignore_genes('-genes'=>\%targets);
  }

  #Produce a table that summarizes the molecular events of each gene, broken down by event type for each clinseq subject
  #Include summary columns that provide samples counts affected by: SNVs, Indels, ..., Any Event Type
  #One gene per row
  #Each sample column will contain a list of observed event types, each indicated by a single letter, word, etc.
  #The user will specify these as options. Indicating it will trigger that kind of event consideration

  $self->print_final_table('-genes'=>\%targets, '-gmap'=>$gmap, '-subject_labels'=>$subject_labels, '-event_labels'=>$event_labels);

  #Warn the user if the number of targets imported is more than can be displaye according to $max_genes
  my $target_count = keys %targets;
  my $max_genes = $self->max_genes;
  $self->warning_message("Imported $target_count genes with events but only $max_genes (--max-genes) will be displayed in the heatmap") if ($target_count > $max_genes);

  #Run the paired R tool to create the final heatmap visualizations from the results tables created
  my $sample_count = keys %{$subject_labels};
  my $all_events_r_script = __FILE__ . '.R';
  my $r_cmd = "$all_events_r_script " . $self->outdir . " $sample_count $max_genes";
  Genome::Sys->shellcmd(cmd => $r_cmd);

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


sub override_clinseq_subject_labels{
  my $self = shift;
  my %args = @_;
  my $old_subject_labels = $args{'-old_subject_labels'};
  my $subject_labels_file = $self->subject_labels_file;
  die $self->error_message("Could not find subject legends table file: $subject_labels_file") unless (-e $subject_labels_file);

  my %subject_labels;

  my $header = 1;
  open (IN, $subject_labels_file) || die $self->error_message("Could not open subject legends table file: $subject_labels_file");
  my %columns;
  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    next unless ($_ =~ /\w+/);
    if ($header){
      my $p = 0;
      foreach my $col (@line){
        $columns{$col}{p} = $p;
        $p++;
      }
      $header = 0;
      unless (defined($columns{'build_id'}) && defined($columns{'name'}) && defined($columns{'name_abr'}) && defined($columns{'order'})){
        die $self->error_message("Subject legends table file must be a TSV with the following columns/headers: build_id, name, name_abr, order");
      }
      next;
    }
    my $bid = $line[$columns{'build_id'}{p}];
    my $name = $line[$columns{'name'}{p}];
    my $name_abr = $line[$columns{'name_abr'}{p}];
    my $order = $line[$columns{'order'}{p}];
    $subject_labels{$bid}{name} = $name;
    $subject_labels{$bid}{name_abr} = $name_abr;
    $subject_labels{$bid}{order} = $order;
  }
  close (IN);

  #Make sure you have parity between these manually specified subject labels and those obtained from the database (i.e. the build ids match)
  unless (keys %subject_labels == keys %{$old_subject_labels}){
    die $self->error_message("The number of subject labels in the subject_labels_file does not match that found in the database.  Update this file!");
  }
  foreach my $bid (keys %subject_labels){
    unless ($old_subject_labels->{$bid}){
      die $self->error_message("The build ids in the subject_labels_file do not match those found in the database for your model-group.  Update this file!");
    }
  }

  return(\%subject_labels);
}


sub check_event_labels{
  my $self = shift;
  
  $self->status_message("Making sure all supplied event labels are unique");
  
  my @event_labs;
  push (@event_labs, $self->snv_label) if $self->snv_label;
  push (@event_labs, $self->indel_label) if $self->indel_label;
  push (@event_labs, $self->cnv_gain_label) if $self->cnv_gain_label;
  push (@event_labs, $self->cnv_loss_label) if $self->cnv_loss_label;
  push (@event_labs, $self->de_up_label) if $self->de_up_label;
  push (@event_labs, $self->de_down_label) if $self->de_down_label;
  push (@event_labs, $self->sv_fusion_label) if $self->sv_fusion_label;
  push (@event_labs, $self->tophat_outlier_label) if $self->tophat_outlier_label;
  push (@event_labs, $self->cufflinks_outlier_label) if $self->cufflinks_outlier_label;
  my $labs_given = scalar(@event_labs);

  my %event_labs;
  $event_labs{$self->snv_label}{name} = 'SNV' if $self->snv_label;
  $event_labs{$self->indel_label}{name} = 'InDel' if $self->indel_label;
  $event_labs{$self->cnv_gain_label}{name} = 'CNV Gain' if $self->cnv_gain_label;
  $event_labs{$self->cnv_loss_label}{name} = 'CNV Loss' if $self->cnv_loss_label;
  $event_labs{$self->de_up_label}{name} = 'DE Up' if $self->de_up_label;
  $event_labs{$self->de_down_label}{name} = 'DE Down' if $self->de_down_label;
  $event_labs{$self->sv_fusion_label}{name} = 'SV Fusion' if $self->sv_fusion_label;
  $event_labs{$self->tophat_outlier_label}{name} = 'Tophat Outlier' if $self->tophat_outlier_label;
  $event_labs{$self->cufflinks_outlier_label}{name} = 'Cufflinks Outlier' if $self->cufflinks_outlier_label;
  my $labs_found = keys %event_labs;

  unless ($labs_given == $labs_found){
    $self->error_message("All event labels supplied must be unique");
    $self->error_message("User supplied: " . "@event_labs");
    die $self->error_message("Aborting");
  }

  #Set default colors
  $event_labs{$self->snv_label}{color} = '#984EA3' if $self->snv_label;
  $event_labs{$self->indel_label}{color} = '#4DAF4A' if $self->indel_label;
  $event_labs{$self->cnv_gain_label}{color} = '#E41A1C' if $self->cnv_gain_label;
  $event_labs{$self->cnv_loss_label}{color} = '#377EB8' if $self->cnv_loss_label;
  $event_labs{$self->de_up_label}{color} = '#FF7F00' if $self->de_up_label;
  $event_labs{$self->de_down_label}{color} = '#A6CEE3' if $self->de_down_label;
  $event_labs{$self->sv_fusion_label}{color} = '#FFFF33' if $self->sv_fusion_label;
  $event_labs{$self->tophat_outlier_label}{color} = '#FB9A99' if $self->tophat_outlier_label;
  $event_labs{$self->cufflinks_outlier_label}{color} = '#FDBF6F' if $self->cufflinks_outlier_label;

  #Display a summary of labels and their names
  foreach my $label (sort {$event_labs{$a}->{name} cmp $event_labs{$b}->{name}} keys %event_labs){
    my $name = $event_labs{$label}{name};
    $self->status_message("\t$label -> $name");
  }
  $self->status_message("\n");

  #Store the labels and their names to a legend file
  my $legend_file = $self->outdir . "/events_legend.txt";
  my $c = 0;
  open (LEG, ">$legend_file") || die $self->error_message("Could not open events legend file for writing: $legend_file");
  print LEG "event_type\tevent_label\tindex\tcolor\n";
  print LEG "None\t"."-"."\t$c\t#FFFFFF\n"; #White for no hits
  foreach my $label (sort {$event_labs{$a}->{name} cmp $event_labs{$b}->{name}} keys %event_labs){
    my $name = $event_labs{$label}{name};
    $c++;
    my $color = $event_labs{$label}{color};
    print LEG "$name\t$label\t$c\t$color\n";
    $event_labs{$label}{index} = $c;
  }
  $c++;
  print LEG "Multiple\tvariable\t$c\t#000000\n"; #Black for multiple hits
  close(LEG);

  return \%event_labs;
}


sub override_event_labels{
  my $self = shift;
  my $event_labels_file = $self->event_labels_file; 
  die $self->error_message("Could not find event labels file: $event_labels_file") unless (-e $event_labels_file);
  my %event_labs;
  my $legend_file = $self->outdir . "/events_legend.txt";
  Genome::Sys->copy_file($event_labels_file, $legend_file);
    
  return;
}


sub add_events{
  my $self = shift;
  my %args = @_;
  my $genes = $args{'-genes'};
  my $gmap = $args{'-gmap'};
  my $file = $args{'-file'};
  my $column_name = $args{'-column_name'};
  my $sample_label = $args{'-sample_label'};
  my $event_label = $args{'-event_label'};
 
  my $path = $file->{path};
  my $columns = $file->{columns};

  unless (defined($columns->{$column_name})){
    die $self->error_message("Expected column $column_name not found in file: $path");
  }

  my $header = 1;
  open (IN, $path) || die $self->error_message("Could not open file: $path");
  while(<IN>){
    if ($header){
      $header = 0;
      next;
    }
    chomp($_);
    my @line = split("\t", $_);
    my $gid = $line[$columns->{$column_name}->{p}];

    #Make sure this gene is recognized
    unless ($gmap->{$gid}){
      $self->warning_message("gene_id: $gid not recognized from file ($path) - skipping");
      next;
    }

    #Add the event to the genes object
    $genes->{$gid}->{$sample_label}->{$event_label} = 1;

  }
  close(IN);

  return;
}


#TODO: Eliminate this non-generic sub-routine once the SV annotation issue is cleaned up
sub add_sv_events{
  my $self = shift;
  my %args = @_;
  my $genes = $args{'-genes'};
  my $gmap = $args{'-gmap'};
  my $file = $args{'-file'};
  my $sample_label = $args{'-sample_label'};
  my $event_label = $args{'-event_label'};
 
  my $path = $file->{path};
  my $columns = $file->{columns};

  unless (defined($columns->{'gene1'}) && defined($columns->{'gene2'})){
    die $self->error_message("Expected columns gene1 and gene2 not found in file: $path");
  }

  #Create a name to ensg map
  my %gmap_name;
  foreach my $ensg_id (sort keys %{$gmap}){
    my $ensg_name = $gmap->{$ensg_id}->{ensg_name};
    $gmap_name{$ensg_name}{ensg_id} = $ensg_id;
  }

  my $header = 1;
  open (IN, $path) || die $self->error_message("Could not open file: $path");
  while(<IN>){
    if ($header){
      $header = 0;
      next;
    }
    chomp($_);
    my @line = split("\t", $_);
    my $gname1 = $line[$columns->{'gene1'}->{p}];
    my $gname2 = $line[$columns->{'gene2'}->{p}];

    #Attempt to map this gene name to an ENSG
    if ($gmap_name{$gname1}){
      #Add the event to the genes object
      my $gid = $gmap_name{$gname1}{ensg_id};
      $genes->{$gid}->{$sample_label}->{$event_label} = 1;
    }
    if ($gmap_name{$gname2}){
      #Add the event to the genes object
      my $gid = $gmap_name{$gname2}{ensg_id};
      $genes->{$gid}->{$sample_label}->{$event_label} = 1;
    }

  }
  close(IN);

  return;
}


sub target_specific_genes{
  my $self = shift;
  my %args = @_;
  my $genes = $args{'-genes'};
  my $gmap = $args{'-gmap'};
  my $target_gene_list = $self->target_gene_list;

  my %target_genes;
  open (GENES, $target_gene_list) || die $self->error_message("Could not open target gene list: $target_gene_list");
  while(<GENES>){
    chomp($_);
    next unless ($_ =~ /ENSG/);
    my @line = split("\t", $_);
    my $gid = $line[0];
    $target_genes{$gid}=1;

    #Make sure this gene is recognized
    unless ($gmap->{$gid}){
      $self->warning_message("gene_id: $gid not recognized from target gene list file ($target_gene_list) - skipping");
      next;
    }

    #Make sure all these genes are added to the master list
    unless ($genes->{$gid}){
      $genes->{$gid} = ();
    }
  }
  close(GENES);

  #Remove any genes that are in the master list but do not match a target
  foreach my $gid (keys %{$genes}){
    delete $genes->{$gid} unless ($target_genes{$gid});
  }

  return;
}


sub remove_ignore_genes{
  my $self = shift;
  my %args = @_;
  my $genes = $args{'-genes'};
  my $ignore_gene_list = $self->ignore_gene_list;

  my %ignore_genes;
  open (GENES, $ignore_gene_list) || die $self->error_message("Could not open ignore gene list: $ignore_gene_list");
  while(<GENES>){
    chomp($_);
    if ($_ =~ /(ENSG\d+)/){
      $ignore_genes{$1}=1;    
    }
  }
  close(GENES);

  #Remove any genes that are in the master list that are also in the ignore list
  foreach my $gid (keys %{$genes}){
    delete $genes->{$gid} if ($ignore_genes{$gid});
  }

  return;
}


sub print_final_table{
  my $self = shift;
  my %args = @_;
  my $genes = $args{'-genes'};
  my $gmap = $args{'-gmap'};
  my $subject_labels = $args{'-subject_labels'};
  my $event_labels = $args{'-event_labels'};

  #Get a string of column names for the subjects
  my @subject_list;
  foreach my $bid (sort {$subject_labels->{$a}->{order} <=> $subject_labels->{$b}->{order}} keys %{$subject_labels}){
    push(@subject_list, $subject_labels->{$bid}->{name_abr});
  }
  my $subject_list_string = join("\t", @subject_list);

  #Store subject counts for each event type for this gene
  my %event_counter;
  foreach my $event (keys %{$event_labels}){
    $event_counter{$event} = 0;
  }

  my @event_counter_names;
  foreach my $event_counter_name (sort keys %event_counter){
    push (@event_counter_names, $event_counter_name . "_subject_count");
  }
  my $event_counter_names_string = join("\t", @event_counter_names);

  my $outfile = $self->outdir . "/events_final.tsv";
  my $outfile_categorical = $self->outdir . "/events_final_categorical.tsv";
  my $outfile_numerical = $self->outdir . "/events_final_numerical.tsv";
  my $outfile_labels = $self->outdir . "/events_final_labels.tsv";

  $self->status_message("\nPrinting final table to $outfile\n");
  open (OUT1, ">$outfile") || die $self->error_message("Could not open outfile: $outfile for writing");
  open (OUT2, ">$outfile_categorical") || die $self->error_message("Could not open outfile: $outfile_categorical for writing");
  open (OUT3, ">$outfile_numerical") || die $self->error_message("Could not open outfile: $outfile_numerical for writing");
  open (OUT4, ">$outfile_labels") || die $self->error_message("Could not open outfile: $outfile_labels for writing");

  #print "ensg_id\tensg_name\t$subject_list_string\n";
  print OUT1 "ensg_id\tensg_name\t$subject_list_string\t$event_counter_names_string\tgrand_subject_count\n";
  print OUT2 "ensg_id\tensg_name\tsubject\tevents\n";
  print OUT3 "ensg_id\tensg_name\t$subject_list_string\n";
  print OUT4 "ensg_id\tensg_name\t$subject_list_string\n";

  foreach my $gid (sort keys %{$genes}){
    my %grand_subject_list;
    my $gene_name = $gmap->{$gid}->{ensg_name};
    my @subject_events;
    my %gene_event_counter = %event_counter;
    foreach my $bid (sort {$subject_labels->{$a}->{order} <=> $subject_labels->{$b}->{order}} keys %{$subject_labels}){
      my $subject = $subject_labels->{$bid}->{name_abr};
      my @events;
      if (defined($genes->{$gid}->{$subject})){
        foreach my $event (keys %{$event_labels}){
          if ($genes->{$gid}->{$subject}->{$event}){
            push(@events, $event);
            $gene_event_counter{$event}++;
            $grand_subject_list{$subject}=1;
          }
        }
      }
      my @events_sorted = sort @events;
      my $events_string = join(",", @events_sorted);
      print OUT2 "$gid\t$gene_name\t$subject\t$events_string\n" if $events_string;
      $events_string = "-" unless $events_string;
      push (@subject_events, $events_string); 
    }

    my @gene_event_counts;
    foreach my $event (sort keys %gene_event_counter){
      push(@gene_event_counts, $gene_event_counter{$event});
    }
    my $gene_event_counts_string = join("\t", @gene_event_counts);

    my $grand_subject_count = keys %grand_subject_list;
    my $subject_events_string = join("\t", @subject_events);
    #print "$gid\t$gene_name\t$subject_events_string\n";
    print OUT1 "$gid\t$gene_name\t$subject_events_string\t$gene_event_counts_string\t$grand_subject_count\n";
    print OUT4 "$gid\t$gene_name\t$subject_events_string\n";

    #Generate a numerical classification of each event class and print as a table
    my @subject_events2;
    my $event_class_count = keys %{$event_labels};
    my $multi_event_class = $event_class_count + 1;
    foreach my $subject_event (@subject_events){
      my $class;
      if ($subject_event eq "-"){
        $class = 0;
      }elsif($event_labels->{$subject_event}){
        $class = $event_labels->{$subject_event}->{index};
      }else{
        $class = $multi_event_class;
      }
      push(@subject_events2, $class);
    }
    my $subject_events_string2 = join("\t", @subject_events2);
    print OUT3 "$gid\t$gene_name\t$subject_events_string2\n";

  }
  close(OUT1);
  close(OUT2);
  close(OUT3);

  return;
}


1;

