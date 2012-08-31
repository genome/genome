package Genome::Model::Tools::Capture::DefineModels;

use warnings;
use strict;
use IO::File;
use Genome;

class Genome::Model::Tools::Capture::DefineModels {
  is => 'Command',
  has_input => [
    ref_align_model_group_id => { is => 'Text', doc => "ID of model-group containing reference-alignments to pair up" },
    output_dir => { is => 'Text', doc => "Output directory to store various useful lists and statistics" },
    som_var_model_group_name => { is => 'Text', doc => "A model-group name for the SomaticVariation models to define", is_optional => 1 },
    processing_profile_name => { is => 'Text', doc => "Processing profile for the SomaticVariation models to define", is_optional => 1, default => "Dec 2011 Default Somatic Variation Exome" },
    annotation_build => { is => 'Text', doc => "Annotation build for the SomaticVariation models to define", is_optional => 1, default => "NCBI-human.combined-annotation/58_37c_v2" },
    dbsnp_build_id => { is => 'Text', doc => "The build that whitelists previously discovered variations", is_optional => 1, default => "110108854" },
    expected_refseq_name => { is => 'Text', doc => "Refseq that reference-alignments are expected to use", is_optional => 1, default => "GRCh37-lite-build37" },
    expected_roi_set_name => { is => 'Text', doc => "ROI set name that reference-alignments are expected to use", is_optional => 1, default => "agilent_sureselect_exome_version_2_broad_refseq_cds_only_hs37" },
  ],
  doc => "Defines SomaticVariation models, given a model-group of tumor-normal ref-alignment pairs",
};

sub help_detail {
  return <<HELP;
Given a model-group containing reference-alignments of tumor and normal samples (usually sequenced
in batches), this tool generates the files listed below, and optionally defines SomaticVariation
models for each tumor-normal pair of reference-alignments models. If som_var_model_group_name is
not specified, then SomaticVariation models are not defined.

IMPORTANT: It is recommended to run this tool first without specifying som_var_model_group_name,
to ensure that the ref_align_model_group contains models that use expected inputs and the latest
available instrument data. Fix/ignore the warnings reported, and then run this tool again with a
som_var_model_group_name specified, to start defining SomaticVariation models.

Files written to output_dir:

paired_model_ids - Lists tumor-normal pairs in the following format:
[Tumor_Sample Tumor_Model_ID Normal_Sample Normal_Model_ID]

build_start_commands - If som_var_model_group_name was specified and models defined, then this
file lists the commands to start builds. Ideally, you don't want to start too many builds at once,
but this tool assumes that you're a responsible grown-up, and won't tell you what not to do.
HELP
}

sub _doc_authors {
  return " Cyriac Kandoth, Ph.D.";
}

sub execute {
  my $self = shift;
  my $refaln_group_id = $self->ref_align_model_group_id;
  my $output_dir = $self->output_dir;
  my $som_var_model_group_name = $self->som_var_model_group_name;
  my $processing_profile_name = $self->processing_profile_name;
  my $annotation_build = $self->annotation_build;
  my $dbsnp_build_id = $self->dbsnp_build_id;
  my $expected_refseq_name = $self->expected_refseq_name;
  my $expected_roi_set_name = $self->expected_roi_set_name;

  $output_dir =~ s/(\/)+$//; # Remove trailing forward-slashes if any

  # A helpful hash, especially if you have more than just pairs of solid tumors and blood normals
  my %lookup_tcga_sample_type = (
    01 => "Primary Solid Tumor", 02 => "Recurrent Solid Tumor",
    03 => "Primary Blood Derived Cancer", 04 => "Recurrent Blood Derived Cancer",
    05 => "Additional New Primary", 06 => "Metastatic", 07 => "Additional Metastatic",
    10 => "Blood Derived Normal", 11 => "Solid Tissue Normal", 12 => "Buccal Cell Normal",
    13 => "EBV Immortalized Normal", 14 => "Bone Marrow Normal", 20 => "Cell Line Control",
  );

  # Check on all the input data before starting work
  my $refaln_group = Genome::ModelGroup->get( $refaln_group_id );
  print STDERR "ERROR: Could not find a model-group with ID: $refaln_group_id\n" unless( defined $refaln_group );
  print STDERR "ERROR: Output directory not found: $output_dir\n" unless( -e $output_dir );
  return undef unless( defined $refaln_group && -e $output_dir );

  my @refaln_models = $refaln_group->models;
  my %cases; # Hash to facilitate the pairing of tumors and normals for each case
  foreach my $refaln_model ( @refaln_models )
  {
    my $sample = $refaln_model->subject_name;
    unless( $sample =~ m/^\w+-\w{2}-\w{4}-\w{2}/ )
    {
      print STDERR "ERROR: Sample name $sample does not match the TCGA format. Unable to pair tumor-normals\n";
      return undef;
    }
    my ( $case, $type ) = $sample =~ m/^\w+-(\w{2}-\w{4})-(\w{2})/;
    unless( $type =~ m/^0\d$|^10$/ )
    {
      print STDERR "WARNING: Skipping ", $lookup_tcga_sample_type{$type}, " sample $sample\n" if defined( $lookup_tcga_sample_type{$type} );
      print STDERR "WARNING: Skipping sample $sample of unknown type\n" unless defined( $lookup_tcga_sample_type{$type} );
      next;
    }

    $sample = "$case-$type"; # Drop a load

    # Find out if this model is using the expected params
    unless( defined $refaln_model->reference_sequence_name and $refaln_model->reference_sequence_name eq $expected_refseq_name )
    {
      print STDERR "WARNING: Unexpected value for reference_sequence_name (", $refaln_model->reference_sequence_name;
      print STDERR ") in model ", $refaln_model->id, " for sample $sample\n";
    }
    unless( defined $refaln_model->region_of_interest_set_name and $refaln_model->region_of_interest_set_name eq $expected_roi_set_name )
    {
      print STDERR "WARNING: Unexpected value for region_of_interest_set_name (", $refaln_model->region_of_interest_set_name;
      print STDERR ") in model ", $refaln_model->id, " for sample $sample\n";
    }

    # Look through all available ref-alignment models to ensure that we have the latest one
    my @models = Genome::Model::ReferenceAlignment->get( subject_name => { operator => 'like' , value=> "\%$sample\%" } );
    my @models_we_can_use = (); # Shortlist those models that use the expected params
    foreach my $model ( @models )
    {
      # If everything about this model checks out, then we can use it for analysis
      if( defined $model->reference_sequence_name and defined $model->region_of_interest_set_name and
          $model->reference_sequence_name eq $expected_refseq_name and
          $model->region_of_interest_set_name eq $expected_roi_set_name )
      {
        push( @models_we_can_use, $model );
      }
    }

    # Look through all the available instrument-data to find out if we're using all of them in ref-alignments
    my @instrument_data = Genome::InstrumentData::Solexa->get( sample_name => { operator => 'like' , value=> "\%$sample\%" } );
    my @expected_inst_data_ids = (); # Exclude those that have been silently deleted
    foreach my $inst_data ( @instrument_data )
    {
      push( @expected_inst_data_ids, $inst_data->id ) if( -s $inst_data->bam_path );
    }
    @expected_inst_data_ids = sort @expected_inst_data_ids;

    my $latest_model;
    if( scalar( @models_we_can_use ) == 0 )
    {
      print STDERR "WARNING: No reference-alignment models found for sample $sample that use the expected parameters\n";
    }
    elsif( scalar( @models_we_can_use ) > 0 )
    {
      $latest_model = $models_we_can_use[0];
      # If there are multiple models, find the most recently defined one (the largest model_id)
      foreach my $model ( @models_we_can_use )
      {
        $latest_model = $model if( $model->id > $latest_model->id );
      }
      my @inst_data_ids = $latest_model->instrument_data_ids;
      @inst_data_ids = sort @inst_data_ids;
      if( $latest_model->id > $refaln_model->id && @inst_data_ids == @expected_inst_data_ids )
      {
        print STDERR "WARNING: Newer ref-aln model ", $latest_model->id, " which uses the same inst-data as ";
        print STDERR $refaln_model->id, " (mdl-grp $refaln_group_id), is available for sample $sample\n";
      }
      elsif( $latest_model->id > $refaln_model->id && @inst_data_ids != @expected_inst_data_ids )
      {
        print STDERR "WARNING: Newer ref-aln model ", $latest_model->id, " which uses different inst-data as ";
        print STDERR $refaln_model->id, " (mdl-grp $refaln_group_id), is available for sample $sample\n";
      }
      elsif( $latest_model->id eq $refaln_model->id and @inst_data_ids ne @expected_inst_data_ids )
      {
        print STDERR "WARNING: Ref-aln model ", $refaln_model->id, " (mdl-grp $refaln_group_id) ";
        print STDERR "for sample $sample does not use all available inst-data\n";
      }
    }

    # After printing any warnings for the user, add this sample to a hash for pairing up later
    $cases{$case}{$type}{wu_name} = $refaln_model->subject_name;
    $cases{$case}{$type}{model_id} = $refaln_model->id;
  }

  # For each patient, print the tumor-normal pairs into a file
  my $outFh = IO::File->new( "$output_dir/paired_model_ids", ">" ) or die "Failed to open file. $!\n";
  foreach my $case ( keys %cases )
  {
    my ( $t_sample, $t_id, $n_sample, $n_id );
    foreach my $type ( keys %{$cases{$case}} )
    {
      ( $t_sample, $t_id ) = ( $cases{$case}{$type}{wu_name}, $cases{$case}{$type}{model_id} ) if( $type =~ m/^0\d$/ );
      ( $n_sample, $n_id ) = ( $cases{$case}{$type}{wu_name}, $cases{$case}{$type}{model_id} ) if( $type eq "10" );
    }
    if( defined $t_sample and defined $t_id and defined $n_sample and defined $n_id )
    {
      $outFh->print( join( "\t", $t_sample, $t_id, $n_sample, $n_id ), "\n" );
    }
    else
    {
      print STDERR "WARNING: Skipping case ", $case, " which does not have a tumor-normal pair of ref-alignments available\n"
    }
  }
  $outFh->close;
  print "Wrote available tumor-normal pairs to $output_dir/paired_model_ids\n";

  if( defined $som_var_model_group_name )
  {
    my @defined_model_ids = ();
    print "\nDefining SomaticVariation models for each case in paired_model_ids...\n";
    my $pairsFh = IO::File->new( "$output_dir/paired_model_ids" ) or die "Failed to open file. $!\n";
    while( my $line = $pairsFh->getline )
    {
      chomp( $line );
      my ( $t_sample, $t_id, $n_sample, $n_id ) = split( /\t/, $line );
      $t_sample =~ s/^\w+-(\w{2}-\w{4}-\w{2}).*$/$1/;
      $n_sample =~ s/^\w+-(\w{2}-\w{4}-\w{2}).*$/$1/;
      my $model_name = "SomaticVariation-$expected_refseq_name-$t_sample\_vs_$n_sample";

      # If someone can write this in a nicer object-oriented way, that'd be great!
      my @stdoe = `genome model define somatic-variation --processing-profile "$processing_profile_name" --annotation-build "$annotation_build" --previously-discovered-variations-build $dbsnp_build_id --tumor-model $t_id --normal-model $n_id --model-name $model_name 2>&1`;
      if( $stdoe[0] =~ m/^Created model/ )
      {
        my ( $model_id ) = $stdoe[2] =~ m/\x1b\[31mid\x1b\[0m: \x1b\[36m(\d+)\x1b\[0m/; # Parse thru that unnecessarily colorful gunk
        print "Defined model $model_name with ID: $model_id\n";
        push( @defined_model_ids, $model_id );
      }
      else
      {
        print STDERR "ERROR: Failed to define model for $t_sample\_vs_$n_sample\n", @stdoe;
        return undef;
      }
    }
    $pairsFh->close;

    # Define a model-group for all the models that were defined above
    my $model_ids = join( " ", @defined_model_ids );
    my @stdoe = `genome model-group create "$som_var_model_group_name" $model_ids 2>&1`;
    if( scalar( grep( m/^Created model group/, @stdoe ) > 0 ))
    {
      my ( $line ) = grep( m/^ID: \d+, NAME: /, @stdoe );
      my ( $model_group_id ) = $line =~ m/^ID: (\d+), NAME: /;
      print "Defined model-group $som_var_model_group_name with ID: $model_group_id\n";
    }
    else
    {
      print STDERR "ERROR: Failed to define model-group $som_var_model_group_name\n", @stdoe;
    }

    # Write a file containing all the SomaticVariation model IDs, and the commands to start builds
    my $buildCmdFh = IO::File->new( "$output_dir/build_start_commands", ">" ) or die "Failed to open file. $!\n";
    $buildCmdFh->print( "# Start 5 builds at a time, around 6 hours apart, to ease load on LSF/disks (unless most software results are ready)\n" );
    for( my $i = 0; $i < scalar( @defined_model_ids ); ++$i )
    {
      $buildCmdFh->print( "\ngenome model build start" ) if( $i % 5 == 0 );
      $buildCmdFh->print( " ", $defined_model_ids[$i] );
    }
    $buildCmdFh->print( "\n" );
    $buildCmdFh->close;
    print "\nCommands to start builds can be found in $output_dir/build_start_commands\n";
  }

  return 1;
}
