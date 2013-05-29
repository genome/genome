package Genome::Model::Tools::Capture::SomaticVariationGroup;

use warnings;
use strict;
use IO::File;
use Genome;

class Genome::Model::Tools::Capture::SomaticVariationGroup {
  is => 'Command',

  has => [
    group_id => { is => 'Text', doc => "Specify the SomaticVariation model-group to work with", is_optional => 0 },
    output_build_dirs => { is => 'Text', doc => "Specify a filename for a list of last succeeded build directories" , is_optional => 1 },
    output_bam_list => { is => 'Text', doc => "Specify a filename for a Tumor-normal BAM list" , is_optional => 1 },
    output_germline_calls => { is => 'Text', doc => "Specify a directory to output germline SNV/indel calls" , is_optional => 1 },
    germline_roi_file => { is => 'Text', doc => "A BED file within which to restrict the reported germline calls" , is_optional => 1 },
    reference => { is => 'Text', doc => "Reference FASTA for use by FP filters", is_optional => 0, example_values => ['build37   /gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa'] },
    reference_transcripts => { is => 'Text', doc => "Annotation build to use to annotate variants" , is_optional => 1, example_values => ['NCBI-human.combined-annotation/58_37c_v2'] },
  ],
  doc => "Perform various operations on a given SomaticVariation model-group",
};

sub help_synopsis {
  return <<HELP;
 gmt capture somatic-variation-group --group-id 24350 \\
    --output-germline-calls germline_calls \\
    --germline-roi-file all_coding_regions.bed
HELP
}

sub help_detail {
  return <<HELP;
Given a model-group of SomaticVariation models, this tool can perform the following on each model:
  - Print data directories of the last succeeded builds
  - Output lists of annotated germline variants that lie within user-specified regions of interest
  - ::TODO:: Other useful things
HELP
}

sub _doc_authors {
  return <<AUTHS;
 Cyriac Kandoth
 Daniel C Koboldt
AUTHS
}

sub execute {
  my $self = shift;
  my $group_id = $self->group_id;
  my %stats = ();

  # If the user wants a lost of data dirs, then open up a file for output
  my $build_dir_fh;
  if( $self->output_build_dirs ) {
    $build_dir_fh = IO::File->new( $self->output_build_dirs, ">" ) or die "Cannot open $self->output_build_dirs for output. $!";
  }

  # Fetch all the models in the model group
  my $model_group = Genome::ModelGroup->get( $group_id );
  my @models = $model_group->models;


	if($self->output_bam_list)
	{
	  open(BAMLIST, ">" . $self->output_bam_list) or die "Can't open BAM list: $!\n";
	}

  foreach my $model ( @models ) {

    unless ( $model->class =~ /Genome::Model::SomaticVariation/ ) {
        die( 'This tool can only be used with Somatic Variation Models! ' . $model->name .
            ' is of type: ' . $model->class . ".\n" );
    };

    $stats{models_in_group}++;

    # Pull information about this SomVar model and all of it's builds
    my $model_id = $model->genome_model_id;
    my $model_name = $model->name;
    my $subject_name = $model->subject_name;
    $subject_name = "Model" . $model_id unless( defined $model->subject_name );
    my @builds = $model->builds;

    # Get normal and tumor models, and look all through their builds
    my $normal_model = $model->normal_model;
    my $tumor_model = $model->tumor_model;
    my $normal_model_id = $normal_model->id;
    my $tumor_model_id = $tumor_model->id;
    my @normal_builds = $normal_model->builds;
    my @tumor_builds = $tumor_model->builds;

    # Get the latest normal and tumor builds
    my ( $normal_build_dir, $tumor_build_dir ) = ( "", "" );

    foreach my $build ( @normal_builds ) {
      my $build_status = $build->status;
      my $build_dir = $build->data_directory;
      $normal_build_dir = $build_dir if( $build_status eq "Succeeded" );
    }

    foreach my $build ( @tumor_builds ) {
      my $build_status = $build->status;
      my $build_dir = $build->data_directory;
      $tumor_build_dir = $build_dir if( $build_status eq "Succeeded" );
    }

	## Output the BAM list ##
      if($self->output_bam_list)
      {
        my $tumor_model = $model->tumor_model;
        my $normal_model = $model->normal_model;
	if($tumor_model->last_succeeded_build && $normal_model->last_succeeded_build)
	{
		my $tumor_bam = $tumor_model->last_succeeded_build->whole_rmdup_bam_file;
		my $normal_bam = $normal_model->last_succeeded_build->whole_rmdup_bam_file;
		print BAMLIST join("\t", $tumor_model_id, $tumor_model->last_succeeded_build->id, "TUMOR_NOBAM") . "\n" if(!$tumor_bam);
		print BAMLIST join("\t", $normal_model_id, $normal_model->last_succeeded_build->id, "NORMAL_NOBAM") . "\n" if(!$normal_bam);

		if($tumor_bam && $normal_bam)
		{
			my $tumor_sample = $tumor_model->subject_name;
			my $normal_sample = $normal_model->subject_name;
			print BAMLIST join("\t", $tumor_sample, $normal_bam, $tumor_bam, $normal_sample) . "\n";
#			print join("\t", $tumor_sample, $normal_bam, $tumor_bam, $normal_sample) . "\n";
		}
	      }		
	}
	elsif($normal_model->last_succeeded_build)
	{
		warn "No Succeeded build for tumor model $tumor_model_id\n";	
	}
	elsif($tumor_model->last_succeeded_build)
	{
		warn "No Succeeded build for normal model $normal_model_id\n";	
	}
	else
	{
		warn "No Succeeded build for tumor model $tumor_model_id\n";
		warn "No Succeeded build for normal model $normal_model_id\n";	
	}


    my $last_build_dir = "";
    my $model_status = "New";
    my $final_build_result = "";
    my $last_build_id = 0;

    if(( $self->output_build_dirs || $self->output_germline_calls ) && scalar( @builds ) > 0 ) {
      my ( $build_ids, $build_statuses ) = ( "", "" );
      $model_status = "Building";

      foreach my $build ( @builds ) {
        my $build_id = $build->id;
        my $build_status = $build->status;
        my $build_dir = $build->data_directory;

        $build_ids .= "," if( $build_ids );
        $build_statuses .= "," if( $build_statuses );
        $build_ids .= $build_id;
        $build_statuses .= $build_status;

        if( $model_status eq "New" || $build_status eq "Succeeded" || $build_status eq "Running" ) {
          $model_status = $build_status;
          $last_build_dir = $build_dir;
        }
      }
      



      if( $model->last_succeeded_build_directory ) {
        $model_status = "Succeeded";  # Override if we have successful build dir
        $last_build_dir = $model->last_succeeded_build_directory;
        my $tumor_model = $model->tumor_model;
        my $normal_model = $model->normal_model;
        my $tumor_sample = $tumor_model->subject_name;
        my $normal_sample = $normal_model->subject_name;

        if( $self->output_build_dirs ) {
          $build_dir_fh->print( "$tumor_sample\t$normal_sample\t$last_build_dir\n" );
        }

        if( $self->output_germline_calls ) {
          my $dir = $self->output_germline_calls;
          unless( -d $dir ) {
            mkdir( $dir ) or die "Can't create output directory: $!\n";
          }

          # Make output dir for this model
	  my $this_model_name = $model_name;
	  $this_model_name =~ s/[^0-9A-Za-z\.\-\_]/\_/g;
          my $germline_dir = $self->output_germline_calls . "/" . $this_model_name;
	  if(-d $germline_dir)
	  {
		
	  }
	  else
	  {
#		unless( -d $germline_dir ) {
		  mkdir( $germline_dir ) or die "Can't create output directory: $!\n";
#		}		

		my $normal_model = $model->normal_model;
		my $tumor_model = $model->tumor_model;
      
		output_germline_files( $self, $last_build_dir, $germline_dir, $normal_model, $tumor_model );

	  }


        }
      }
    }
    print join( "\t", $model_id, $model_name, $normal_model_id, $tumor_model_id, $final_build_result ) . "\n";
  }
    
    close(BAMLIST) if($self->output_bam_list);

  # Close any open file handles
  $build_dir_fh->close if( $self->output_build_dirs );

  print $stats{models_in_group} . " models in group\n" if( defined $stats{models_in_group} );
  print $stats{models_running} . " models running\n" if( defined $stats{models_running} );
  print $stats{models_finished} . " models finished\n" if( defined $stats{models_finished} );

  return 1;
}

################################################################################################
# Output Germline Files - output the files for germline analysis
################################################################################################

sub output_germline_files
{
  my ( $self, $build_dir, $germline_dir, $normal_model, $tumor_model ) = @_;

  # Get File of Somatic SNVs
  my $somatic_snvs = "$build_dir/variants/snvs.hq.bed";

  # Get Normal/Tumor Models and BAMs
  my $normal_model_dir = $normal_model->last_succeeded_build_directory;
  my $normal_bam = $normal_model->last_succeeded_build->whole_rmdup_bam_file;
  warn "Could not locate Normal Ref-alignment BAM within build dir $normal_model_dir\n" unless( -e $normal_bam );

  my $tumor_model_dir = $tumor_model->last_succeeded_build_directory;
  my $tumor_bam = $tumor_model->last_succeeded_build->whole_rmdup_bam_file;
  warn "Could not locate Tumor Ref-alignment BAM within build dir $tumor_model_dir\n" unless( -e $tumor_bam );

  ###################################################################
  ### SAMTOOLS REFALIGN VARIANTS
  ###################################################################

  # Get variant file to parse
  my ( $normal_variant_file ) = glob("$normal_model_dir/*/filtered.variants.post_annotation");
  my ( $tumor_variant_file ) = glob("$tumor_model_dir/*/filtered.variants.post_annotation");

  my $outfile_tumor_snp = $germline_dir . "/samtools.tumor.snp";
  my $outfile_normal_snp = $germline_dir . "/samtools.normal.snp";

  # Parse out Normal SNPs
  warn "Parsing Normal SAMtools SNPs file...\n";
  my $num_normal_snps = get_samtools_snps($normal_variant_file, $outfile_normal_snp);
  print "$num_normal_snps SNPs obtained for normal\n";

   if($self->germline_roi_file){
    warn "Limiting SNPs to ROI...\n";
    # Limit SNPs to Germline ROI BED
    system("java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.jar limit $outfile_normal_snp --regions-file " . $self->germline_roi_file . " --output-file $outfile_normal_snp.roi");
    $outfile_normal_snp .= ".roi";
  }

  # Parse out Tumor SNPs
  warn "Parsing Tumor SAMtools SNPs file...\n";
  my $num_tumor_snps = get_samtools_snps($tumor_variant_file, $outfile_tumor_snp);
  print "$num_tumor_snps SNPs obtained for tumor\n";

  if($self->germline_roi_file) {
    warn "Limiting SNPs to ROI...\n";
    # Limit SNPs to Germline ROI BED
    system("java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.jar limit $outfile_tumor_snp --regions-file " . $self->germline_roi_file . " --output-file $outfile_tumor_snp.roi");
    $outfile_tumor_snp .= ".roi";
  }

  # Remove somatic SNPs
  warn "Removing Somatic SNVs from normal list...\n";
  system("java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.jar limit $outfile_normal_snp --regions-file " . $somatic_snvs . " --not-file $outfile_normal_snp.germline");
  $outfile_normal_snp .= ".germline";

  warn "Removing Somatic SNVs from tumor list...\n";
  system("java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.jar limit $outfile_tumor_snp --regions-file " . $somatic_snvs . " --not-file $outfile_tumor_snp.germline");
  $outfile_tumor_snp .= ".germline";

  # Apply FP-filter to SNPs, and also annotate them
  warn "Submitting a job to apply FP-filter on normal SNPs...\n";
  my $cmd = "bsub -R 'select[mem>8000] rusage[mem=8000]' -M 8000000 'gmt somatic filter-false-positives ";
  $cmd .= "--reference " . $self->reference . " ";
  $cmd .= "--variant-file $outfile_normal_snp --bam-file $normal_bam --output-file $outfile_normal_snp.fpfilter ";
  $cmd .= "--filtered-file $outfile_normal_snp.fpfilter.removed --max-mm-qualsum-diff 100; ";
  $cmd .= "gmt annotate transcript-variants --annotation-filter top --reference-transcripts " . $self->reference_transcripts . " ";
  $cmd .= "--variant-file $outfile_normal_snp.fpfilter --output-file $outfile_normal_snp.fpfilter.anno'";
  system($cmd);

  warn "Submitting a job to apply FP-filter on tumor SNPs...\n";
  $cmd = "bsub -R 'select[mem>8000] rusage[mem=8000]' -M 8000000 'gmt somatic filter-false-positives ";
  $cmd .= "--reference " . $self->reference . " ";
  $cmd .= "--variant-file $outfile_tumor_snp --bam-file $tumor_bam --output-file $outfile_tumor_snp.fpfilter ";
  $cmd .= "--filtered-file $outfile_tumor_snp.fpfilter.removed --max-mm-qualsum-diff 100; ";
  $cmd .= "gmt annotate transcript-variants --annotation-filter top --reference-transcripts " . $self->reference_transcripts . " ";
  $cmd .= "--variant-file $outfile_tumor_snp.fpfilter --output-file $outfile_tumor_snp.fpfilter.anno'";
  system($cmd);

  ###################################################################
  ### VARSCAN GERMLINE/LOH VARIANTS
  ###################################################################

  warn "Parsing VarScan Germline/LOH SNP files...\n";
  my ( $varscan_germline_snps ) = glob("$build_dir/variants/snv/varscan-somatic-*/varscan-high-confidence-*/snvs.Germline.hc");
  my ( $varscan_loh_snps ) = glob("$build_dir/variants/snv/varscan-somatic-*/varscan-high-confidence-*/snvs.LOH.hc");

  if(-e $varscan_germline_snps && -e $varscan_loh_snps) {
    my $output_snp_varscan = $germline_dir . "/varScan.germline.snp";
    system("cat $varscan_germline_snps $varscan_loh_snps >$output_snp_varscan");
    system("gmt capture sort-by-chr-pos --input $output_snp_varscan --output $output_snp_varscan");

    if($self->germline_roi_file) {
      warn "Limiting SNPs to ROI...\n";
      system("java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.jar limit $output_snp_varscan --regions-file " . $self->germline_roi_file . " --output-file $output_snp_varscan.roi");
      $output_snp_varscan .= ".roi";
    }

    # Apply FP-filter using Normal BAM, and annotate the output
    warn "Submitting a job to apply FP-filter on SNPs, and annotate the output...\n";
    my $cmd = "bsub -R 'select[mem>8000] rusage[mem=8000]' -M 8000000 'gmt somatic filter-false-positives ";
    $cmd .= "--reference " . $self->reference . " --variant-file $output_snp_varscan --bam-file $normal_bam ";
    $cmd .= "--output-file $output_snp_varscan.fpfilter --filtered-file $output_snp_varscan.fpfilter.removed --max-mm-qualsum-diff 100; ";
    $cmd .= "gmt annotate transcript-variants --annotation-filter top --reference-transcripts " . $self->reference_transcripts . " ";
    $cmd .= "--variant-file $output_snp_varscan.fpfilter --output-file $output_snp_varscan.fpfilter.anno'";
    system($cmd);
  }
  else {
    warn "Warning: Could not locate VarScan Germline/LOH SNPs within $build_dir\n$varscan_germline_snps\n$varscan_loh_snps\n";
  }

  ###################################################################
  ### GATK GERMLINE INDELS
  ###################################################################

  warn "Parsing GATK Germline Indel files...\n";
  my ( $gatk_germline_indels ) = glob("$build_dir/variants/indel/gatk-somatic-indel-*/gatk_output_file");

  if(-e $gatk_germline_indels) {
    my $output_indel_gatk = $germline_dir . "/gatk.germline.indel";
    system("grep GERMLINE $gatk_germline_indels >$output_indel_gatk");

    # If ROIs are defined, select only the indels within
    if($self->germline_roi_file) {
      warn "Limiting Indels to ROI...\n";
      system("java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.jar limit $output_indel_gatk --regions-file " . $self->germline_roi_file . " --output-file $output_indel_gatk.roi");
      $output_indel_gatk .= ".roi";
    }

    # Convert the GATK output format to something that the annotator can parse
    my $cmd_obj = Genome::Model::Tools::Gatk::FormatIndels->create(
      variants_file => $output_indel_gatk,
      output_file => "$output_indel_gatk.bed",
    );
    $cmd_obj->execute;

    # Apply FP-filter using Normal BAM, and annotate the output
    warn "Submitting a job to apply FP-filter on GATK Indels, and annotate the output...\n";
    my $cmd = "bsub -R 'select[mem>8000] rusage[mem=8000]' -M 8000000 'gmt somatic filter-false-indels ";
    $cmd .= "--reference " . $self->reference . " --variant-file $output_indel_gatk.bed --bam-file $normal_bam ";
    $cmd .= "--output-file $output_indel_gatk.fpfilter --filtered-file $output_indel_gatk.fpfilter.removed --max-mm-qualsum-diff 100; ";
    $cmd .= "gmt annotate transcript-variants --annotation-filter top --reference-transcripts " . $self->reference_transcripts . " ";
    $cmd .= "--variant-file $output_indel_gatk.fpfilter --output-file $output_indel_gatk.fpfilter.anno'";
    system($cmd);
  }
  else {
    warn "Warning: Could not locate GATK indels within $build_dir\n$gatk_germline_indels\n";
  }

  ###################################################################
  ### VARSCAN GERMLINE/LOH INDELS
  ###################################################################

  warn "Parsing VarScan Germline/LOH Indel files...\n";
  my ( $varscan_germline_indels ) = glob("$build_dir/variants/indel/varscan-somatic-*/varscan-high-confidence-indel-*/indels.Germline.hc");
  my ( $varscan_loh_indels ) = glob("$build_dir/variants/indel/varscan-somatic-*/varscan-high-confidence-indel-*/indels.LOH.hc");

  if(-e $varscan_germline_indels && -e $varscan_loh_indels) {
    my $output_indel_varscan = $germline_dir . "/varScan.germline.indel";
    system("cat $varscan_germline_indels $varscan_loh_indels >$output_indel_varscan");
    system("gmt capture sort-by-chr-pos --input $output_indel_varscan --output $output_indel_varscan");

    if($self->germline_roi_file)
    {
      warn "Limiting Indels to ROI...\n";
      system("java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.jar limit $output_indel_varscan --regions-file " . $self->germline_roi_file . " --output-file $output_indel_varscan.roi");
      $output_indel_varscan .= ".roi";
    }

    # Convert the VarScan Indel output format to something that the annotator can parse
    my $cmd_obj = Genome::Model::Tools::Capture::FormatIndels->create(
      variants_file => $output_indel_varscan,
      output_file => "$output_indel_varscan.bed",
    );
    $cmd_obj->execute;

    # Apply FP-filter using Normal BAM, and annotate the output
    warn "Submitting a job to apply FP-filter on Indels, and annotate the output...\n";
    my $cmd = "bsub -R 'select[mem>8000] rusage[mem=8000]' -M 8000000 'gmt somatic filter-false-indels ";
    $cmd .= "--reference " . $self->reference . " --variant-file $output_indel_varscan.bed --bam-file $normal_bam ";
    $cmd .= "--output-file $output_indel_varscan.fpfilter --filtered-file $output_indel_varscan.fpfilter.removed --max-mm-qualsum-diff 100; ";
    $cmd .= "gmt annotate transcript-variants --annotation-filter top --reference-transcripts " . $self->reference_transcripts . " ";
    $cmd .= "--variant-file $output_indel_varscan.fpfilter --output-file $output_indel_varscan.fpfilter.anno'";
    system($cmd);
  }
  else {
    warn "Warning: Could not locate VarScan germline Indels within $build_dir\n$varscan_germline_indels\n$varscan_loh_indels\n";
  }
}

################################################################################################
# Get SAMtools SNPs - get the SNPs from SAMtools filtered post-annotation file
################################################################################################

sub get_samtools_snps
{
  my ( $variant_file, $output_file ) = @_;
  my $num_snps = 0;

  if(-e $variant_file) {
    my $output = IO::File->new( $output_file, ">" ) or die "Cannot open $output_file. $!";
    my $input = IO::File->new( $variant_file ) or die "Cannot open $variant_file. $!";

    while ( <$input> ) {
      chomp;
      my $line = $_;
      my ( $chrom, $chr_start, $chr_stop, $ref, $var, $var_type ) = split( /\t/, $line );

      if( $var_type eq "SNP" ) {
        $num_snps++;
        $output->print( join( "\t", $chrom, $chr_start, $chr_stop, $ref, $var ), "\n" );
      }
    }

    $input->close;
    $output->close;
  }
  else {
    warn "Warning: Could not locate Samtools variant file $variant_file\n";
  }

  return( $num_snps );
}

1;
