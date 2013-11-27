package Genome::Model::ClinSeq::Command::DumpIgvXml;
#Written by Malachi Griffith

use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq::Util;

use Genome::Utility::List;

class Genome::Model::ClinSeq::Command::DumpIgvXml {
    is => 'Command::V2',
    has_input => [
        builds => { 
              is => 'Genome::Model::Build::ClinSeq',
              is_many => 1,
              shell_args_position => 1,
              require_user_verify => 0,
              doc => 'clinseq build to summarize',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
    ],
    doc => 'Based on the inputs of a clinseq build create a series of IGV session XML files',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq dump-igv-xml  --outdir=/tmp/  119971814

genome model clin-seq dump-igv-xml  --outdir=/tmp/  id=119971814

genome model clin-seq dump-igv-xml  --outdir=/tmp/  'id in [119971814,119971932]'

EOS
}

sub help_detail {
    return <<EOS

Create XML files to allow easy loading of IGV sessions with varying amounts of detail

For example: BAM files only, including junctions.bed, including SNV and INDEL calls, etc.

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
  return @errors;
}


#Get the builds that underlie the input clinseq build: wgs somatic variation (and their ref aligns), exome somatic variation (and their ref aligns), rna-seq

#Get the following file 'resources'
#WGS tumor/normal BAM, Exome tumor/normal BAM, RNA-seq tumor/normal BAM
# - Each BAM file will be used to create an IGV panel with two tracks: 1 coverage + 1 mapped reads
#WGS somatic variation SNV and INDEL files (tier1-4, 'novel' and 'previously discovered')
#Exome somatic variation SNV and INDEL files (tier1-4, 'novel' and 'previously discovered')
#RNA-seq junctions file

#To label each track we will need to retrieve or determine the following information for each sample
# - sample type: tumor/normal?
# - data type: wgs/exome/rna-seq?
# - sample name

sub get_levels {
  my $self = shift;
  #For different levels of output verbosity, generate an IGV XML session file

  #Possible 'BED' tracks to utilize from WGS and Exome somatic variation results:
  #SNVs (repeated for 'indels')
  #effects/snvs.hq.novel.tier1.v2.bed, effects/snvs.hq.novel.tier2.v2.bed, effects/snvs.hq.novel.tier3.v2.bed, effects/snvs.hq.novel.tier4.v2.bed
  #effects/snvs.hq.previously_detected.tier1.v2.bed, effects/snvs.hq.previously_detected.tier2.v2.bed, effects/snvs.hq.previously_detected.tier3.v2.bed, effects/snvs.hq.previously_detected.tier4.v2.bed
  #effects/snvs.lq.tier1.v2.bed, effects/snvs.lq.tier2.v2.bed, effects/snvs.lq.tier3.v2.bed, effects/snvs.lq.tier4.v2.bed
  #variants/snvs.hq.bed
  #variants/snvs.lq.bed

  #Level 1 list:
  my @level1 = qw (effects/snvs.hq.novel.tier1.v2.bed variants/snvs.hq.bed);

  #Level 2 list:
  my @level2 = qw (effects/snvs.hq.novel.tier1.v2.bed variants/snvs.hq.bed
                   effects/indels.hq.novel.tier1.v2.bed variants/indels.hq.bed);

  #Level 3 list
  my @level3 = qw (effects/snvs.hq.novel.tier1.v2.bed effects/snvs.hq.novel.tier2.v2.bed effects/snvs.hq.novel.tier3.v2.bed effects/snvs.hq.novel.tier4.v2.bed
                   effects/indels.hq.novel.tier1.v2.bed effects/indels.hq.novel.tier2.v2.bed effects/indels.hq.novel.tier3.v2.bed effects/indels.hq.novel.tier4.v2.bed);

  #Level 4 list:
  my @level4 = qw (effects/snvs.hq.novel.tier1.v2.bed variants/snvs.hq.bed variants/snvs.lq.bed);

  #Level 5 list:
  my @level5 = qw (effects/snvs.hq.novel.tier1.v2.bed variants/snvs.hq.bed variants/snvs.lq.bed effects/indels.hq.novel.tier1.v2.bed variants/indels.hq.bed variants/indels.lq.bed);

  #Level 6 list
  my @level6 = qw (effects/snvs.hq.novel.tier1.v2.bed effects/snvs.hq.novel.tier2.v2.bed effects/snvs.hq.novel.tier3.v2.bed effects/snvs.hq.novel.tier4.v2.bed 
                   effects/snvs.hq.previously_detected.tier1.v2.bed effects/snvs.hq.previously_detected.tier2.v2.bed effects/snvs.hq.previously_detected.tier3.v2.bed effects/snvs.hq.previously_detected.tier4.v2.bed
                   effects/indels.hq.novel.tier1.v2.bed effects/indels.hq.novel.tier2.v2.bed effects/indels.hq.novel.tier3.v2.bed effects/indels.hq.novel.tier4.v2.bed
                   effects/indels.hq.previously_detected.tier1.v2.bed effects/indels.hq.previously_detected.tier2.v2.bed effects/indels.hq.previously_detected.tier3.v2.bed effects/indels.hq.previously_detected.tier4.v2.bed);


  #Level 7 list
  my @level7 = qw (effects/snvs.hq.novel.tier1.v2.bed effects/snvs.hq.novel.tier2.v2.bed effects/snvs.hq.novel.tier3.v2.bed effects/snvs.hq.novel.tier4.v2.bed
                   effects/snvs.hq.previously_detected.tier1.v2.bed effects/snvs.hq.previously_detected.tier2.v2.bed effects/snvs.hq.previously_detected.tier3.v2.bed effects/snvs.hq.previously_detected.tier4.v2.bed
                   effects/indels.hq.novel.tier1.v2.bed effects/indels.hq.novel.tier2.v2.bed effects/indels.hq.novel.tier3.v2.bed effects/indels.hq.novel.tier4.v2.bed
                   effects/indels.hq.previously_detected.tier1.v2.bed effects/indels.hq.previously_detected.tier2.v2.bed effects/indels.hq.previously_detected.tier3.v2.bed effects/indels.hq.previously_detected.tier4.v2.bed
                   variants/snvs.lq.bed variants/indels.lq.bed);

  my %levels;
  $levels{'1'} = \@level1;
  $levels{'2'} = \@level2;
  $levels{'3'} = \@level3;
  $levels{'4'} = \@level4;
  $levels{'5'} = \@level5;
  $levels{'6'} = \@level6;
  $levels{'7'} = \@level7;

  return(\%levels);
}

sub execute {
  my $self = shift;
  my @clinseq_builds = $self->builds;
  my $outdir = $self->outdir;

  $self->queue_status_messages(1);

  my $levels = $self->get_levels;
  unless ($outdir =~ /\/$/){
    $outdir .= "/";
  }

  #Go through each clinseq build and create XML session files
  my $clinseq_build_count = scalar(@clinseq_builds);

  foreach my $clinseq_build (@clinseq_builds){
    my $clinseq_build_id = $clinseq_build->id;

    my $reference_build = Genome::Model::ClinSeq::Util::resolve_reference_sequence_build($clinseq_build);
    my $reference_genome_name = $reference_build->name;

    #Hardcoded list of allowed reference builds and their mappings to names used in IGV
    my $genome_build = "";
    my $gene_track_name = "";
    my $starting_locus = "";
    if ($reference_genome_name=~/GRCh37\-lite\-build37/){ #Assumes that if name is "like" this you want b37
      $genome_build = "b37";
      $starting_locus = "12:25398182-25398361";
    }elsif ($reference_genome_name=~/GRCh37\-lite\-human\-build37/){ #Assumes that if name is "like" this you want b37
      $genome_build = "b37";
      $starting_locus = "12:25398182-25398361";
    }else{
      die $self->error_message("Unrecognized reference genome name ($reference_genome_name) supplied to DumpIgvXml.pm");
    }
    $gene_track_name = $genome_build . "_genes";
    
    my $panel_count = 1;
        
    #If there is more than one clinseq build supplied.  Create sub-directories for each
    my $build_outdir;
    if ($clinseq_build_count > 1){
       $build_outdir = $outdir . $clinseq_build_id . "/";
    }else{
       $build_outdir = $outdir;
    }

    $self->status_message("\n\nWriting IGV XML session files for clinseq build $clinseq_build_id to: $build_outdir");

    #XML blocks to be created.  Each should be of the resource type 'bam', 'junction', or 'bed'
    my @resources;
    my $main_features_track_xml = "";
    my $extra_features_track_xml = "";

    if ($clinseq_build->wgs_build){
      #WGS somvar build - source of wgs SNV/Indel BED feature tracks
      my $wgs_somvar_build = $clinseq_build->wgs_build;

      #WGS normal refalign build - source of normal WGS BAMs track
      my $wgs_normal_refalign_build = $wgs_somvar_build->normal_build if ($wgs_somvar_build);
      my $wgs_normal_track = $self->generate_track_xml('-build'=>$wgs_normal_refalign_build, '-resource_type'=>'bam');
      my $wgs_normal_xml = $wgs_normal_track->{xml};
      $wgs_normal_xml = $self->generate_panel_xml('-track_xml'=>$wgs_normal_xml, '-panel_count'=>$panel_count++);
      push(@resources, $wgs_normal_track->{resource});
      $main_features_track_xml .= $wgs_normal_xml;

      #WGS tumor refalign build - source of tumor WGS BAMs track
      my $wgs_tumor_refalign_build = $wgs_somvar_build->tumor_build if ($wgs_somvar_build);
      my $wgs_tumor_track = $self->generate_track_xml('-build'=>$wgs_tumor_refalign_build, '-resource_type'=>'bam');
      my $wgs_tumor_xml = $wgs_tumor_track->{xml};
      $wgs_tumor_xml = $self->generate_panel_xml('-track_xml'=>$wgs_tumor_xml, '-panel_count'=>$panel_count++);
      push(@resources, $wgs_tumor_track->{resource});
      $main_features_track_xml .= $wgs_tumor_xml;

    }
    if ($clinseq_build->exome_build){
      #Exome somvar build - source of exome SNV/Indel BED feature tracks
      my $exome_somvar_build = $clinseq_build->exome_build;
      

      #Exome normal refalign build - source of normal Exome BAMs track
      my $exome_normal_refalign_build = $exome_somvar_build->normal_build if ($exome_somvar_build);
      my $exome_normal_track = $self->generate_track_xml('-build'=>$exome_normal_refalign_build, '-resource_type'=>'bam');
      my $exome_normal_xml = $exome_normal_track->{xml};
      $exome_normal_xml = $self->generate_panel_xml('-track_xml'=>$exome_normal_xml, '-panel_count'=>$panel_count++);
      push(@resources, $exome_normal_track->{resource});
      $main_features_track_xml .= $exome_normal_xml;

      #Exome tumor refalign build - source of tumor Exome BAMs track
      my $exome_tumor_refalign_build = $exome_somvar_build->tumor_build if ($exome_somvar_build);
      my $exome_tumor_track = $self->generate_track_xml('-build'=>$exome_tumor_refalign_build, '-resource_type'=>'bam');
      my $exome_tumor_xml = $exome_tumor_track->{xml};
      $exome_tumor_xml = $self->generate_panel_xml('-track_xml'=>$exome_tumor_xml, '-panel_count'=>$panel_count++);
      push(@resources, $exome_tumor_track->{resource});
      $main_features_track_xml .= $exome_tumor_xml;

    }
    if ($clinseq_build->normal_rnaseq_build){
      #RNAseq normal RNA build - source of normal RNAseq BAMs track
      my $rnaseq_normal_build = $clinseq_build->normal_rnaseq_build;
      my $rnaseq_normal_track = $self->generate_track_xml('-build'=>$rnaseq_normal_build, '-resource_type'=>'bam');
      my $rnaseq_normal_xml = $rnaseq_normal_track->{xml};
      $rnaseq_normal_xml = $self->generate_panel_xml('-track_xml'=>$rnaseq_normal_xml, '-panel_count'=>$panel_count++);
      push(@resources, $rnaseq_normal_track->{resource});
      $main_features_track_xml .= $rnaseq_normal_xml;

      #Also the source of the normal junctions.bed track
      my $rnaseq_normal_junctions_track = $self->generate_track_xml('-build'=>$rnaseq_normal_build, '-resource_type'=>'junctions');
      my $rnaseq_normal_track_xml = $rnaseq_normal_junctions_track->{xml};
      push(@resources, $rnaseq_normal_junctions_track->{resource});
      $extra_features_track_xml .= $rnaseq_normal_track_xml;
    }
    if ($clinseq_build->tumor_rnaseq_build){
      #RNAseq tumor RNA build - source of tumor RNAseq BAMs track
      my $rnaseq_tumor_build = $clinseq_build->tumor_rnaseq_build;
      my $rnaseq_tumor_track = $self->generate_track_xml('-build'=>$rnaseq_tumor_build, '-resource_type'=>'bam');
      my $rnaseq_tumor_xml = $rnaseq_tumor_track->{xml};
      $rnaseq_tumor_xml = $self->generate_panel_xml('-track_xml'=>$rnaseq_tumor_xml, '-panel_count'=>$panel_count++);
      push(@resources, $rnaseq_tumor_track->{resource});
      $main_features_track_xml .= $rnaseq_tumor_xml;

      #Also the source of the tumor junctions.bed track
      my $rnaseq_tumor_junctions_track = $self->generate_track_xml('-build'=>$rnaseq_tumor_build, '-resource_type'=>'junctions');
      my $rnaseq_tumor_track_xml = $rnaseq_tumor_junctions_track->{xml};
      push(@resources, $rnaseq_tumor_junctions_track->{resource});
      $extra_features_track_xml .= $rnaseq_tumor_track_xml;
    }


    #For each level of detail (increasing amounts of feature tracks) produce an output XML file
    foreach my $l (sort {$a <=> $b} keys %{$levels}){

      #Start with the existing extra features tracks and add more
      my $extra_features_track_xml_level = $extra_features_track_xml;
      my @resources_level = @resources;

      my @bed_files = @{$levels->{$l}};
      foreach my $bed_file (@bed_files){
        if ($clinseq_build->wgs_build){
          #WGS somvar build - source of wgs SNV/Indel BED feature tracks
          my $wgs_somvar_build = $clinseq_build->wgs_build;
          my $extra_track = $self->generate_track_xml('-build'=>$wgs_somvar_build, '-resource_type'=>'bed', '-bed_file'=>$bed_file, '-bed_data_type'=>'WGS');
          if ($extra_track){
            my $extra_track_xml = $extra_track->{xml};
            push(@resources_level, $extra_track->{resource});
            $extra_features_track_xml_level .= $extra_track_xml;
          }
        }
        if ($clinseq_build->exome_build){
          #Exome somvar build - source of exome SNV/Indel BED feature tracks
          my $exome_somvar_build = $clinseq_build->exome_build;
          my $extra_track = $self->generate_track_xml('-build'=>$exome_somvar_build, '-resource_type'=>'bed', '-bed_file'=>$bed_file, '-bed_data_type'=>'Exome');
          if ($extra_track){
            my $extra_track_xml = $extra_track->{xml};
            push(@resources_level, $extra_track->{resource});
            $extra_features_track_xml_level .= $extra_track_xml;
          }
        }
      }

      #Generate the resources XML block
      my $resource_track_xml = $self->generate_resource_xml('-resources'=>\@resources_level);

      #Generate the 'FeaturePanel' content
      my $feature_panel_xml = $self->generate_feature_panel_xml('-extra_features_track_xml'=>$extra_features_track_xml_level, '-gene_track_name'=>$gene_track_name);

      #Generate the 'PanelLayout' content based on the number of panels to be rendered (one per resource, plus the feature panel)
      my $panel_layout_xml = $self->generate_panel_layout_xml('-panel_count'=>$panel_count);

      #Put all the XML blocks together
      my $final_xml = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>";
      $final_xml .= "\n\n<Session genome=\"$genome_build\" locus=\"$starting_locus\" version=\"4\">";
      $final_xml .= "\n\n$resource_track_xml";
      $final_xml .= "\n$main_features_track_xml";
      $final_xml .= "\n\n$feature_panel_xml";
      $final_xml .= "\n\n$panel_layout_xml";
      $final_xml .= "\n\n</Session>";
      #print "\n\n\n$final_xml\n\n\n";
      my $outfile_name = $build_outdir . "level$l"."_igv_session.xml";
      my $outfile = IO::File->new(">$outfile_name");
      $outfile->print($final_xml);
    }
  }

  #Summarize the tracks used for each level
  foreach my $l (sort {$a <=> $b} keys %{$levels}){
    my @bed_files = @{$levels->{$l}};
    $self->status_message("\n\nLevel $l IGV session file will contain the following tracks (where they are available):");
    foreach my $bed_file (@bed_files){
      $self->status_message("$bed_file");
    }
  }
  $self->status_message("\n\n");

  my @output = $self->status_messages();
  my $log = IO::File->new(">$outdir/DumpIgvXml.log.txt");
  $log->print(join("\n", @output));

  return 1;
}

####################################################################################################
#Generate XML session track lines for a single 'file resource' for a single build                  #
####################################################################################################
sub generate_track_xml {
  my $self = shift;
  my %args = @_;
  my $build = $args{'-build'};
  my $resource_type = $args{'-resource_type'};
  my $bed_file = $args{'-bed_file'};
  my $bed_data_type = $args{'-bed_data_type'};

  #Input build could be a 'reference alignment', 'rna seq', or 'somatic variation'

  my $build_id = $build->id;
  my $model = $build->model;
  my $model_id = $model->id;
  my $pp_id = $model->processing_profile_id;
  my $pp = Genome::ProcessingProfile->get($pp_id);
  my $pp_name = $pp->name;
  my $pp_type = $pp->type_name;
  my $ref_name = "n/a";
  if ($model->can("reference_sequence_build")){
     my $rb = $model->reference_sequence_build;
     $ref_name = $rb->name;
  }

  my $subject = $build->subject;
  my $subject_name = $subject->name;
  my $build_dir = $build->data_directory;
  my $tissue_desc = "tissue_unknown";
  my $extraction_type = "extraction_unknown";
  my $common_name = "n/a";
  if ($subject->can("common_name")){
    $common_name = $build->subject->common_name || "n/a";
  }
  if ($subject->can("tissue_desc")){
    $tissue_desc = $build->subject->tissue_desc || "n/a";
  }
  if ($subject->can("extraction_type")){
    $extraction_type = $build->subject->extraction_type || "n/a";
  }

  #Infer whether the data is wgs or exome based on the presence of target region set names
  my $sequence_type = "n/a"; 
  if ($pp_type eq "reference alignment"){
    my @lanes = $build->instrument_data;
    my $lane_count = scalar(@lanes);
    my $trsn_count = 0;
    for my $id (@lanes) {
      if ($id->target_region_set_name){
        $trsn_count++;
      }
    }
    if ($trsn_count == $lane_count){
      $sequence_type = "Exome";
    }elsif($trsn_count == 0){
      $sequence_type = "WGS";
    }else{
      $sequence_type = "Mixed";
    }
  }

  #Determine the data type: ('WGS', 'Exome', 'WGS+Exome', 'RNA-Seq')
  my $data_type;
  if ($pp_type eq "reference alignment"){
    $data_type = $sequence_type;
  }elsif($pp_type eq "somatic variation"){
    $data_type = "Unknown";
    if ($pp_name =~ /exome/i){
      $data_type = "Exome";
    }
    if ($pp_name =~ /wgs/i){
      $data_type = "WGS";
    }
  }elsif($pp_type eq "rna seq"){
    $data_type = "RNA-Seq";
  }else{
    die $self->error_message("Processing profile type not recognized: $pp_type");
  }

  #Determine the tissue description: ('normal', 'tumor', 'somatic')
  #Should be defined for most samples but for somatic variation results, change it to 'somatic' to indicate a comparison of tissues
  if ($pp_type eq "somatic variation"){
    $tissue_desc = "somatic";
  }

  #Get the file resource according to the resource type specified
  my $resource_file;
  my $resource_file_url;
  if ($resource_type eq 'bam'){
    my $bam_file;
    if ($build->can("whole_rmdup_bam_file")){
      $bam_file = $build->whole_rmdup_bam_file;
    }elsif($build->can("alignment_result")){
      my $alignment_result = $build->alignment_result;
      if ($alignment_result->can("bam_file")){
        $bam_file = $alignment_result->bam_file;
      }
    }else{
      $self->error_message("Could not find BAM file for build: $build_id\tmodel: $model_id - did you specify the correct resource type?");
    }
    $resource_file = $bam_file;

  }elsif ($resource_type eq 'bed' && $pp_type eq 'somatic variation'){
    #User specifies actual file name - will be used to find file in a somatic variation result
    unless ($bed_file){
      die $self->error_message("A bed file name must be defined in dump-igv-xml to generate a track of type 'bed'");
    }
    $resource_file = $build_dir . "/$bed_file";

    #Unless this file is actually present.  Print a warning and return
    unless (-e $resource_file){
      $self->status_message("BED file not found for this somatic variation build:\n$bed_file\n");
      return(0);
    }
  }elsif ($resource_type eq 'junctions' && $pp_type eq 'rna seq'){
    if (-e $build_dir . "/alignments/junctions.bed"){
      $resource_file = $build_dir . "/alignments/junctions.bed";
    }elsif(-e $build_dir . "/junctions/junctions.bed"){
      $resource_file = $build_dir . "/junctions/junctions.bed";
    }
    unless (-e $resource_file){
      die $self->error_message("junctions.bed file not found for rna seq build: $build_dir");
    }
  }else{
    die $self->error_message("Resource type and PP type combination not recognized: $resource_type (must be one of 'bam', 'junctions' (rna-seq models only), or 'bed')");
  }


  #XML parameters
  my $min = "0.0";
  my $max = "100.0";
  my $font_size = 11;
  my $color_option = "UNEXPECTED_PAIR";

  #Set the resource file url path
  $resource_file_url = Genome::Utility::List::join_with_single_slash($ENV{GENOME_SYS_SERVICES_FILES_URL}, $resource_file);

  #Abbreviate some aspects of the base track name
  if ($tissue_desc =~ /skin\,\s+nos/i){
    $tissue_desc = "skin";
  }
  if ($extraction_type =~ /genomic\s+dna/i){
    $extraction_type = "gDNA";
  }
  if ($extraction_type =~ /rna/i){
    $extraction_type = "RNA";
  }
  if ($common_name =~ /normal/i){
    $common_name = "Normal";
  }
  if ($common_name =~ /tumor/i){
    $common_name = "Tumor";
  }

  my $base_track_name = "$data_type $common_name ($tissue_desc - $extraction_type)";

  #Create xml according to the type of resource file specified:
  my $xml;

  if ($resource_type eq 'bam'){
    #WGS, Exome or RNA-seq BAM files, SNV/InDel bed files, and junction bed files
    my $coverage_track_name = "$base_track_name Coverage";
    my $read_track_name = "$base_track_name Reads";
    my $resource_file_coverage = $resource_file . "_coverage";

    my $resource_file_coverage_url = Genome::Utility::List::join_with_single_slash($ENV{GENOME_SYS_SERVICES_FILES_URL}, $resource_file_coverage);

  $xml=<<XML;
    <Track altColor="0,0,178" autoScale="true" color="175,175,175" colorScale="ContinuousColorScale;0.0;9062.0;255,255,255;175,175,175" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="$font_size" id="$resource_file_coverage_url" name="$coverage_track_name" showDataRange="true" visible="true">
      <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="$max" minimum="$min" type="LINEAR"/>
    </Track>
    <Track altColor="0,0,178" color="0,0,178" colorOption="$color_option" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="$font_size" id="$resource_file_url" name="$read_track_name" showDataRange="true" sortByTag="" visible="true"/>
XML

  }elsif ($resource_type eq 'junctions'){
    my $read_track_name = "$base_track_name Junctions";

  $xml=<<XML;
    <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="$font_size" height="60" id="$resource_file_url" name="$read_track_name" showDataRange="true" visible="true" windowFunction="count"/>
XML

  }elsif ($resource_type eq 'bed'){
    my $bed_file_name;
    if ($bed_file =~ /^\w+\/(.*)/){
      $bed_file_name = $1;
    }else{
      die $self->error_message("Could not determine bed file name from bed file parameter: $bed_file");
    }

    my $bed_track_name = "$bed_data_type $bed_file_name";
  $xml=<<XML;
    <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="$font_size" height="45" id="$resource_file_url" name="$bed_track_name" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>
XML
  }

  #Summarize values found from genome API
  #$self->status_message("\n\nbid: $build_id\tmid: $model_id\tpp_type: $pp_type\tref_name: $ref_name");
  #$self->status_message("subject: $subject_name\tcommon_name: $common_name\ttissue_desc: $tissue_desc\textraction_type: $extraction_type");
  #$self->status_message("build_dir: $build_dir");
  #$self->status_message("data_type: $data_type");
  #$self->status_message("resource_type: $resource_type\tresource_file_url: $resource_file_url");

  #Return the xml block and the resource file
  my %result;
  $result{xml} = $xml;
  $result{resource} = $resource_file_url;

  return (\%result);
}

####################################################################################################
#Generate an XML track for the resources                                                           #
####################################################################################################
sub generate_resource_xml {
  my $self = shift;
  my %args = @_;
  my @resource_list = @{$args{'-resources'}};

  unless (scalar(@resource_list) > 0){
    die $self->error_message("\n\nNo resources found for creation of a resource XML section for IGV in DumpIgvXml.pm\n\n");
  }

  #Create an XML block like this:
  # <Resources>
  #   <Resource path="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc7001/info/model_data/2879615516/build114445127/alignments/114469152.bam"/>
  #   <Resource path="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc7001/info/model_data/2879616958/build114445247/alignments/114445014.bam"/>
  #   <Resource path="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc2016/info/model_data/2880794613/build115909698/alignments/accepted_hits.bam"/>
  #   <Resource path="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc2016/info/model_data/2880794613/build115909698/alignments/junctions.bed"/>
  # </Resources>

  my $xml = "  <Resources>";
  foreach my $resource (sort @resource_list){
    $xml .= "\n    <Resource path=\"$resource\"/>";
  }
  $xml .= "\n  </Resources>";

  return($xml);
}


####################################################################################################
#Wrap XML tracks into 'Panel' XML                                                                  #
####################################################################################################
sub generate_panel_xml {
  my $self = shift;
  my %args = @_;
  my $track_xml = $args{'-track_xml'};
  my $panel_count = $args{'-panel_count'};

  my $panel_name = "Panel" . $panel_count++;
  my $panel_xml = '';

  $panel_xml .= "\n  <Panel name=\"$panel_name\">";
  $panel_xml .= "\n$track_xml";
  $panel_xml .= "  </Panel>";

  return($panel_xml);
}


####################################################################################################
#Create the 'FeaturePanel' XML                                                                     #
####################################################################################################
sub generate_feature_panel_xml {
  my $self = shift;
  my %args = @_;
  my $extra_features_track_xml = $args{'-extra_features_track_xml'};
  my $gene_track_name = $args{'-gene_track_name'};

  #XML parameters
  my $min = "0.0";
  my $max = "100.0";
  my $font_size = 11;

  my $xml = '';
  $xml .= "  <Panel name=\"FeaturePanel\">";
  $xml .= "\n     <Track altColor=\"0,0,178\" color=\"0,0,178\" displayMode=\"COLLAPSED\" featureVisibilityWindow=\"-1\" fontSize=\"$font_size\" id=\"Reference sequence\" name=\"Reference sequence\" showDataRange=\"true\" sortable=\"false\" visible=\"true\"/>";
  $xml .= "\n     <Track altColor=\"0,0,178\" color=\"0,0,178\" colorScale=\"ContinuousColorScale;0.0;160.0;255,255,255;0,0,178\" displayMode=\"COLLAPSED\" featureVisibilityWindow=\"-1\" fontSize=\"$font_size\" height=\"35\" id=\"$gene_track_name\" name=\"Gene\" renderer=\"BASIC_FEATURE\" showDataRange=\"true\" sortable=\"false\" visible=\"true\" windowFunction=\"count\">";
  $xml .= "\n       <DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"$max\" minimum=\"$min\" type=\"LINEAR\"/>";
  $xml .= "\n     </Track>";
  $xml .= "\n$extra_features_track_xml";
  $xml .= "\n  </Panel>";
        
  return($xml);
}


####################################################################################################
#Generate the 'PanelLayout' content based on the number of panels to be rendered                   #
####################################################################################################
sub generate_panel_layout_xml{
  my $self = shift;
  my %args = @_;
  my $panel_count = $args{'-panel_count'};
 
  my $xml;

  #Every session should have at least 2 panels (the feature panel + at least one data panel) otherwise there may be an input problem and/or the igv session would 
  #If for some reason only the feature panel was present, simply return and empty string for the panel layout as a layout would be uneccessary
  if ($panel_count == 0){
    die $self->error_message("\n\nNo panels defined for creation panel layout XML in DumpIgvXml.pm\n\n");
  }elsif($panel_count == 1){
    $xml = '';
  }else{
    #If an IGV session has say 6 panels (including the feature panel), it will have 6 entries in this list
    #The first entry always seems to be near 0.  The last entry is 1 - (the desired relative space used by the feature panel)
    # <PanelLayout dividerFractions="0.007352941176470588,0.32107843137254904,0.6017156862745098,0.8762254901960784"/>
    my @divs;
    my $current_divider = 0;
    $current_divider = sprintf("%.3f", $current_divider);
    push(@divs, $current_divider);
    
    my $extra_feature_space = 0.1;
    my $space_left = 1 - $extra_feature_space;
    my $spacer = $space_left/$panel_count;
    for (my $i = 1; $i < $panel_count; $i++){
      $current_divider += $spacer;
      $current_divider = sprintf("%.3f", $current_divider);
      push(@divs, $current_divider);
    }
    my $div_string = join(",", @divs);
    $xml = "  <PanelLayout dividerFractions=\"$div_string\"/>";
  }

  return($xml);
}

#NOTES on IGV XML.

#Header line
#<?xml version="1.0" encoding="UTF-8" standalone="no"?>

#Level 1 tag is a 'Session':
#<Session genome="b37" locus="3:124447213-124447214" version="4">
#$content
#</Session>

#Level 2 contains 4 tags: 'Resources', 'Panel', 'PanelLayout', and the optional 'Regions' 

#Example IGV XML for 'Resources' section
# <Resources>
#   <Resource path="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc7001/info/model_data/2879615516/build114445127/alignments/114469152.bam"/>
#   <Resource path="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc7001/info/model_data/2879616958/build114445247/alignments/114445014.bam"/>
#   <Resource path="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc2016/info/model_data/2880794613/build115909698/alignments/accepted_hits.bam"/>
#   <Resource path="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc2016/info/model_data/2880794613/build115909698/alignments/junctions.bed"/>
# </Resources>

#Example IGV XML for 'Panel' section.  Data panels have these long random seeming numeric names. Each file should also have a Panel with the name 'FeaturePanel'
# <Panel height="3561" name="Panel1341709698396" width="1613">
#   $track_content
# </Panel>

#Panels seem to only contain 'Track' records.  Some 'Track' records also contain a 'DataRange' record.  
#It seems that a 'Track' marked as 'displayMode="COLLAPSED"' must have the 'DataRange' record
#It also seems that a '<Track /> record only has a closing </Track> tag if it contains a 'DataRange' record

#Example IGV XML for 'PanelLayout' - divider fractions determine the relative position of panel boundaries (each with a value from 0.0 - 1.0)
#If an IGV session has say 6 panels (including the feature panel), it will have 6 entries in this list
#The first entry always seems to be near 0.  The last entry is 1 - the desired relative space used by the feature panel
# <PanelLayout dividerFractions="0.007352941176470588,0.32107843137254904,0.6017156862745098,0.8762254901960784"/>

#Example IGV XML for a FeaturePanel 'Panel' - bed tracks will be added to this panel
# <Panel height="96" name="FeaturePanel" width="1613">
#   <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" showDataRange="true" sortable="false" visible="true"/>
#   <Track altColor="0,0,178" color="0,0,178" colorScale="ContinuousColorScale;0.0;160.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="35" id="b37_genes" name="Gene" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count">
#     <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="160.0" minimum="0.0" type="LINEAR"/>
#   </Track>
# </Panel>

#Example IGV XML for a reference alignment or RNA-seq BAM 'Track'
# <Track altColor="0,0,178" autoScale="true" color="175,175,175" colorScale="ContinuousColorScale;0.0;9062.0;255,255,255;175,175,175" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc2016/info/model_data/2880794613/build115909698/alignments/accepted_hits.bam_coverage" name="RNA-seq Tumor Coverage" showDataRange="true" visible="true">
#   <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="273.0" minimum="0.0" type="LINEAR"/>
# </Track>
# <Track altColor="0,0,178" color="0,0,178" colorOption="UNEXPECTED_PAIR" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" id="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc2016/info/model_data/2880794613/build115909698/alignments/accepted_hits.bam" name="RNA-seq Tumor Reads" showDataRange="true" sortByTag="" visible="true"/>

#Example IGV XML for a junctions.bed 'Track'.  These tracks will be placed in the 'FeaturePanel'
# <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="60" id="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc2016/info/model_data/2880794613/build115909698/alignments/junctions.bed" name="RNA-seq Tumor Junctions" showDataRange="true" visible="true" windowFunction="count"/>

#Example IGV XML for a snvs.bed 'Track'
# <Track altColor="0,0,178" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="12" height="45" id="$ENV{GENOME_SYS_SERVICES_FILES_URL}gscmnt/gc8002/info/model_data/2882504846/build119390903/effects/snvs.hq.novel.tier1.v2.bed" name="WGS snvs.hq.novel.tier1" renderer="BASIC_FEATURE" showDataRange="true" sortable="false" visible="true" windowFunction="count"/>

#Example IGV XML for a 'Regions' of interest section
# <Regions>
#   <Region chromosome="3" description="UMPS" end="124470119" start="124447213"/>
#   <Region chromosome="7" description="EGFR" end="55277031" start="55084725"/>
# </Regions>


1;

