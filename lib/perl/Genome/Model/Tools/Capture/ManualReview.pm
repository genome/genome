package Genome::Model::Tools::Capture::ManualReview;

use warnings;
use strict;
use IO::File;
use Genome;
use Sort::Naturally qw(nsort);

class Genome::Model::Tools::Capture::ManualReview {
  is => 'Command',
  has_input => [
    som_var_model_group_id => { is => 'Text', doc => "ID of model-group containing SomaticVariation models with variants to manually review" },
    output_dir => { is => 'Text', doc => "Directory where manual review files will be written to or read from" },
    exclude_pindel => { is => 'Boolean', doc => "Keep aside indels unique to Pindel, since many cannot be reviewed on a BWA aligned BAM", is_optional => 1, default => 1 },
    exclude_uhf_snvs => { is => 'Boolean', doc => "Keep aside SNV calls that pass the ultra-high-confidence filter (UHF)", is_optional => 1, default => 1 },
    max_variants => { is => 'Number', doc => "Maximum allowed SNVs+Indels (after pindel/UHF exclusions) to consider a case for review", is_optional => 1, default => 250 },
    refseq_version => { is => 'Text', doc => "Reference sequence to use with the UHF (GRCh37-lite-build37 or NCBI-human-build36)", is_optional => 1, example_values => ["GRCh37-lite-build37"] },
    read_review => { is => 'Boolean', doc => "Read existing manual review files and create WU annotation files per case", is_optional => 1, default => 0 },
  ],
  doc => "Reads/creates manual review files, given a model-group of SomaticVariation models",
};

sub help_detail {
  return <<HELP;
Given a model-group of SomaticVariation models, this tool will gather resulting tier1 variants and
prepare them for manual review. Existing review files will not be overwritten. Review files for
cases with more than --max-variants will be written to a subdirectory named too_many_to_review.

After manual review, this tool can be used again with the argument --read-review to parse the
reviewed files (expects file extensions *.reviewed.csv), and print somatic SNVs and Indels in
WashU annotation format for each case (variants with review codes 'S' or 'V').

--exclude-pindel
  If set, calls unique to Pindel are stored in a separate review file. Pindel calls that are also
  found by GATK are stored in the Indel review file (along with calls unique to GATK).

--exclude-uhf-snvs
  If set, the Ultra High-confidence Filter (UHF) is invoked. Calls that pass the filter are stored
  separately in annotation format, and the remaining SNVs are written to the SNV review file.
HELP
}

sub _doc_authors {
  return <<AUTHS;
 Cyriac Kandoth, Ph.D.
 David Larson, Ph.D.
AUTHS
}

sub execute {
  my $self = shift;
  my $som_var_model_group_id = $self->som_var_model_group_id;
  my $output_dir = $self->output_dir;
  my $exclude_pindel = $self->exclude_pindel;
  my $exclude_uhf_snvs = $self->exclude_uhf_snvs;
  my $max_variants = $self->max_variants;
  my $refseq_version = $self->refseq_version;
  my $read_review = $self->read_review;
  $output_dir =~ s/(\/)+$//; # Remove trailing forward-slashes if any

  # Check on all the input data before starting work
  my $somvar_group = Genome::ModelGroup->get( $som_var_model_group_id );
  print STDERR "ERROR: Could not find a model-group with ID: $som_var_model_group_id\n" unless( defined $somvar_group );
  print STDERR "ERROR: Output directory not found: $output_dir\n" unless( $output_dir and -e $output_dir );
  return undef unless( defined $somvar_group && -e $output_dir );

  my @somvar_models = $somvar_group->models;
  my %bams; # Hash to store the tumor-normal BAM pairs
  print "Finding latest succeeded builds in model-group ", $somvar_group->name, "...\n";
  foreach my $model ( @somvar_models )
  {
    my $build = $model->last_succeeded_build;
    unless( defined $build )
    {
      print STDERR "WARNING: Skipping model ", $model->id, " that has no succeeded builds\n";
      next;
    }
    my $tcga_patient_id = $build->tumor_build->model->subject->extraction_label;
    $tcga_patient_id = $build->tumor_build->model->subject_name unless( $tcga_patient_id =~ /^TCGA/ );
    my $tumor_bam = $build->tumor_build->whole_rmdup_bam_file;
    my $normal_bam = $build->normal_build->whole_rmdup_bam_file;

    if( exists( $bams{$tcga_patient_id} ))
    {
      print STDERR "ERROR: Multiple models in model-group $som_var_model_group_id for sample $tcga_patient_id";
      return undef;
    }
    else
    {
      $bams{$tcga_patient_id}{tumor} = $tumor_bam;
      $bams{$tcga_patient_id}{normal} = $normal_bam;
    }
    my $build_dir = $build->data_directory;

    # Check if the necessary SNV and Indel files were created by this build
    my ( $snv_anno ) = glob( "$build_dir/effects/snvs.hq.tier1.v1.annotated.top" );
    print STDERR "ERROR: Tier1 SNV annotations for $tcga_patient_id not found\n" unless( $snv_anno and -e $snv_anno );
    my ( $indel_anno ) = glob( "$build_dir/effects/indels.hq.tier1.v1.annotated.top" );
    print STDERR "ERROR: Tier1 Indel annotations for $tcga_patient_id not found\n" unless( $indel_anno and -e $indel_anno );
    my ( $gatk_indels ) = glob( "$build_dir/variants/indel/gatk-somatic-indel-*/indels.hq.bed" );
    print STDERR "ERROR: GATK calls for $tcga_patient_id not found\n" unless( $gatk_indels and -e $gatk_indels );
    my ( $varscan_indels ) = glob( "$build_dir/variants/indel/varscan-somatic-*/varscan-high-confidence-indel-*/indels.hq.bed" );
    unless( $varscan_indels and -e $varscan_indels ) {
      print STDERR "WARNING: VarScan indel calls for $tcga_patient_id not found\n";
      $varscan_indels = "";
    }
    my ( $pindel_indels ) = glob( "$build_dir/variants/indel/pindel-*/pindel-somatic-calls-*/pindel-vaf-filter-*/pindel-read-support-*/indels.hq.bed" );
    unless( $pindel_indels and -e $pindel_indels ) {
      print STDERR "WARNING: Pindel calls for $tcga_patient_id not found\n";
      $pindel_indels = "";
    }
    return undef unless( -e $snv_anno && -e $indel_anno && -e $gatk_indels );

    $bams{$tcga_patient_id}{snvs} = $snv_anno;
    $bams{$tcga_patient_id}{indels} = $indel_anno;
    if( $exclude_pindel )
    {
      $bams{$tcga_patient_id}{gatk} = $gatk_indels;
      $bams{$tcga_patient_id}{varscan_indel} = $varscan_indels;
      $bams{$tcga_patient_id}{pindel} = $pindel_indels;
    }
  }

  # Unless it already exists, create a subdirectory to keep aside cases with too many variants
  mkdir "$output_dir/too_many_to_review" unless( -e "$output_dir/too_many_to_review" );

  foreach my $case ( keys %bams )
  {
    # If we're parsing reviewed files, then don't create any new ones
    if( $read_review )
    {
      print "Reading review files for case $case... ";
      my ( $snv_uhf_anno, $snv_review, $indel_review, $out_anno_file ) = ( "", "", "", "" );
      if( -e "$output_dir/$case.snv.reviewed.csv" and -e "$output_dir/$case.indel.reviewed.csv" )
      {
        $snv_uhf_anno = "$output_dir/$case.snv.uhf.anno";
        $snv_review = "$output_dir/$case.snv.reviewed.csv";
        $indel_review = "$output_dir/$case.indel.reviewed.csv";
        $out_anno_file = "$output_dir/$case.reviewed.anno";
      }
      elsif( -e "$output_dir/too_many_to_review/$case.snv.reviewed.csv" and -e "$output_dir/too_many_to_review/$case.indel.reviewed.csv" )
      {
        $snv_uhf_anno = "$output_dir/too_many_to_review/$case.snv.uhf.anno";
        $snv_review = "$output_dir/too_many_to_review/$case.snv.reviewed.csv";
        $indel_review = "$output_dir/too_many_to_review/$case.indel.reviewed.csv";
        $out_anno_file = "$output_dir/too_many_to_review/$case.reviewed.anno";
      }
      elsif( -e "$output_dir/too_many_to_review/$case.snv.review.csv" and -e "$output_dir/too_many_to_review/$case.indel.review.csv" )
      {
        print "Case is in folder 'too_many_to_review' with no reviewed files found for SNVs and/or Indels. Skipping.\n";
        next;
      }
      else
      {
        print "No reviewed files found for SNVs and/or Indels. Skipping.\n";
        next;
      }

      # Grab the high confidence tier1 SNV and Indel annotations from their respective files
      my ( $snv_anno, $indel_anno ) = ( $bams{$case}{snvs}, $bams{$case}{indels} );
      my @snv_lines = `cat $snv_anno`;
      my @indel_lines = `cat $indel_anno`;
      chomp( @snv_lines, @indel_lines );

      # Store the variants into a hash to help search through them quickly
      my %anno_lines = ();
      for my $line ( @snv_lines )
      {
        my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
        $anno_lines{snvs}{$chr}{$start}{$stop} = $line;
      }
      for my $line ( @indel_lines )
      {
        my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
        $anno_lines{indels}{$chr}{$start}{$stop} = $line;
      }

      my $out_anno_fh = IO::File->new( $out_anno_file, ">" ) or die "Cannot open $out_anno_file. $!";

      # Read in the UHF SNVs, if found, and write them straight to output
      if ( $snv_uhf_anno and -e $snv_uhf_anno )
      {
        my @uhf_snv_lines = `cat $snv_uhf_anno`;
        chomp( @uhf_snv_lines ); # Chomping and reprinting newlines cuz we shouldn't make assumptions
        $out_anno_fh->print( join( "\n", @uhf_snv_lines ), "\n" );
      }

      # Load the manual review files and print WU annotation lines for call codes 'S' or 'V'. If it has no code, let it through anyway.
      ( @snv_lines, @indel_lines ) = ((), ());
      @snv_lines = `cat $snv_review` if( -e $snv_review );
      @indel_lines = `cat $indel_review` if( -e $indel_review );
      chomp( @snv_lines, @indel_lines );
      for my $line ( @snv_lines )
      {
        my ( $chr, $start, $stop, undef, undef, $call ) = split( /\t/, $line );
        next if( $chr eq 'Chr' ); # Skip header
        if( !defined $call or $call =~ m/^\s*$/ or $call =~ m/[SVsv]/ )
        {
          $out_anno_fh->print( $anno_lines{snvs}{$chr}{$start}{$stop}, "\n" ) if( defined $anno_lines{snvs}{$chr}{$start}{$stop} );
        }
      }
      for my $line ( @indel_lines )
      {
        my ( $chr, $start, $stop, undef, undef, $call ) = split( /\t/, $line );
        next if( $chr eq 'Chr' ); # Skip header
        if( !defined $call or $call =~ m/^\s*$/ or $call =~ m/[SVsv]/ )
        {
          $out_anno_fh->print( $anno_lines{indels}{$chr}{$start}{$stop}, "\n" ) if( defined $anno_lines{indels}{$chr}{$start}{$stop} );
        }
      }
      $out_anno_fh->close;
			$out_anno_file =~ s/$output_dir/output-dir/;
      print "wrote annotation files to $out_anno_file\n";

      next;
    }

    # If we're not parsing reviewed files, then go ahead and create new ones unless they exist
    print "\nPreparing review files for $case... ";

    # Check if any review files exist. We don't want to overwrite reviewed variants
    if( -e "$output_dir/$case.snv.review.csv" or -e "$output_dir/$case.indel.review.csv" or
        -e "$output_dir/$case.snv.reviewed.csv" or -e "$output_dir/$case.indel.reviewed.csv" or
        -e "$output_dir/too_many_to_review/$case.snv.review.csv" or -e "$output_dir/too_many_to_review/$case.indel.review.csv" or
        -e "$output_dir/too_many_to_review/$case.snv.reviewed.csv" or -e "$output_dir/too_many_to_review/$case.indel.reviewed.csv" )
    {
      print "files exist. Will not overwrite.\n";
      next;
    }

    # Find the indels unique to Pindel, if the user wants to keep them aside
    my %uniq_to_pindel = ();
    if( $exclude_pindel )
    {
      my ( $gatk_indels, $varscan_indels, $pindel_indels ) = ( $bams{$case}{gatk}, $bams{$case}{varscan_indel}, $bams{$case}{pindel} );
      # I know I shouldn't use backticks like this, but think of all the precious lines of code we're saving!
      my %gatk_varscan_lines = map {chomp; s/\*/0/; $_=>1} `cut -f 1-4 $varscan_indels $gatk_indels`;
      %uniq_to_pindel = map {chomp; $_=>1} `cut -f 1-4 $pindel_indels` if( $pindel_indels and -e $pindel_indels );
      foreach my $indel ( keys %gatk_varscan_lines )
      {
        delete $uniq_to_pindel{$indel} if( defined $uniq_to_pindel{$indel} );
      }
    }

    # Grab the high confidence calls from their respective files
    my ( $snv_anno, $indel_anno ) = ( $bams{$case}{snvs}, $bams{$case}{indels} );
    my @snv_lines = `cat $snv_anno`;
    my @indel_lines = `cat $indel_anno`;
    chomp( @snv_lines, @indel_lines );

    # Store the variants into a hash to help sort variants by loci, and to remove duplicates
    my %review_lines = ();
    my ( $snv_cnt, $indel_cnt ) = ( 0, 0 );
    for my $line ( @snv_lines )
    {
      my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
      ++$snv_cnt unless( defined $review_lines{snvs}{$chr}{$start}{$stop} );
      $review_lines{snvs}{$chr}{$start}{$stop} = $line; # Save annotation for use with the UHF
    }
    for my $line ( @indel_lines )
    {
      my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
      my ( $base0_start, $base0_stop ) = ( $start - 1, $stop - 1 );
      my $refvar = ( $ref eq '-' ? "0/$var" : "$ref/0" );
      if( $exclude_pindel && ( defined $uniq_to_pindel{"$chr\t$base0_start\t$stop\t$refvar"} ||
                               defined $uniq_to_pindel{"$chr\t$start\t$base0_stop\t$refvar"} ))
      {
        $review_lines{pindels}{$chr}{$start}{$stop} = join( "\t", $chr, $start, $stop, $ref, $var );
      }
      else
      {
        ++$indel_cnt unless( defined $review_lines{indels}{$chr}{$start}{$stop} );
        $review_lines{indels}{$chr}{$start}{$stop} = join( "\t", $chr, $start, $stop, $ref, $var );
      }
    }

    my ( $tumor_bam, $normal_bam ) = ( $bams{$case}{tumor}, $bams{$case}{normal} );

    # If user wants, filter out the ultra-high-confidence SNVs
    my $uhf_snv_file = "$output_dir/$case.snv.uhf.anno";
    if( $exclude_uhf_snvs )
    {
      print "Running ultra-high-confidence filter. This could take a while...\n";
      # Print out a de-duplicated snv annotation file for use with the UHF
      my $snv_anno_file = Genome::Sys->create_temp_file_path();
      my $snv_anno_fh = IO::File->new( $snv_anno_file, ">" ) or die "Cannot open $snv_anno_file. $!";
      for my $chr ( nsort keys %{$review_lines{snvs}} )
      {
        for my $start ( sort {$a <=> $b} keys %{$review_lines{snvs}{$chr}} )
        {
          for my $stop ( sort {$a <=> $b} keys %{$review_lines{snvs}{$chr}{$start}} )
          {
            $snv_anno_fh->print( $review_lines{snvs}{$chr}{$start}{$stop}, "\n" );
          }
        }
      }
      $snv_anno_fh->close;
      my $snv_filtered_file = Genome::Sys->create_temp_file_path();

      # Fetch the path to the reference sequence FASTA file to use with the UHF
      my $refseq_build = Genome::Model::Build::ReferenceSequence->get( name => $refseq_version );
      my $refseq_fasta = $refseq_build->data_directory . "/all_sequences.fa";

      # Run the UHF
      `gmt somatic ultra-high-confidence --variant-file $snv_anno_file --normal-bam-file $normal_bam --tumor-bam-file $tumor_bam --reference $refseq_fasta --output-file $uhf_snv_file --filtered-file $snv_filtered_file`;
      `rm -f $uhf_snv_file.readcounts.normal $uhf_snv_file.readcounts.tumor`; # Remove intermediate files

      # Only the variants that didn't pass the UHF will need to be reviewed
      undef %{$review_lines{snvs}}; # Reset the review line hash
      my @snv_filtered_lines = `cat $snv_filtered_file`;
      chomp( @snv_filtered_lines );
      $snv_cnt = 0;
      for my $line ( @snv_filtered_lines )
      {
        my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
        ++$snv_cnt unless( defined $review_lines{snvs}{$chr}{$start}{$stop} );
        $review_lines{snvs}{$chr}{$start}{$stop} = $line;
      }
    }

    # If there are more variants than our max-variants threshold, store review files separately
    my $output_dir_new = $output_dir;
    if(( $snv_cnt + $indel_cnt ) > $max_variants )
    {
      $output_dir_new = "$output_dir/too_many_to_review";
      # Move over files we've already created for this case
      if( $exclude_uhf_snvs )
      {
        my $uhf_snv_file_new = "$output_dir_new/$case.snv.uhf.anno";
        `mv $uhf_snv_file $uhf_snv_file_new`;
      }
    }

    # Create review files for manual reviewers to record comments, and bed files for use with IGV
    my $snv_review_file = "$output_dir_new/$case.snv.review.csv";
    my $indel_review_file = "$output_dir_new/$case.indel.review.csv";
    my $snv_bed_file = "$output_dir_new/$case.snv.bed";
    my $indel_bed_file = "$output_dir_new/$case.indel.bed";
    # If user wants to exclude calls unique to pindel, create a separate file (no need for a .bed)
    my $pindel_review_file = "$output_dir_new/$case.pindel.review.csv";

    # Write SNV review files
    my $snv_review_fh = IO::File->new( $snv_review_file, ">" ) or die "Cannot open $snv_review_file. $!";
    my $snv_bed_fh = IO::File->new( $snv_bed_file, ">" ) or die "Cannot open $snv_bed_file. $!";
    $snv_review_fh->print( "Chr\tStart\tStop\tRef\tVar\tCall\tNotes\n" );
    for my $chr ( nsort keys %{$review_lines{snvs}} )
    {
      for my $start ( sort {$a <=> $b} keys %{$review_lines{snvs}{$chr}} )
      {
        for my $stop ( sort {$a <=> $b} keys %{$review_lines{snvs}{$chr}{$start}} )
        {
          my ( undef, undef, undef, $ref, $var ) = split( /\t/, $review_lines{snvs}{$chr}{$start}{$stop} );
          $snv_review_fh->print( join( "\t", $chr, $start, $stop, $ref, $var ), "\n" );
          $snv_bed_fh->printf( "%s\t%d\t%d\t%s\t%s\n", $chr, $start-1, $stop, $ref, $var );
        }
      }
    }
    $snv_bed_fh->close;
    $snv_review_fh->close;

    # Write indel review files
    my $indel_review_fh = IO::File->new( $indel_review_file, ">" ) or die "Cannot open $indel_review_file. $!";
    my $indel_bed_fh = IO::File->new( $indel_bed_file, ">" ) or die "Cannot open $indel_bed_file. $!";
    $indel_review_fh->print( "Chr\tStart\tStop\tRef\tVar\tCall\tNotes\n" );
    for my $chr ( nsort keys %{$review_lines{indels}} )
    {
      for my $start ( sort {$a <=> $b} keys %{$review_lines{indels}{$chr}} )
      {
        for my $stop ( sort {$a <=> $b} keys %{$review_lines{indels}{$chr}{$start}} )
        {
          my ( undef, undef, undef, $ref, $var ) = split( /\t/, $review_lines{indels}{$chr}{$start}{$stop} );
          $indel_review_fh->print( $review_lines{indels}{$chr}{$start}{$stop}, "\n" );
          $indel_bed_fh->printf( "%s\t%d\t%d\t%s\t%s\n", $chr, $start-1, $stop, $ref, $var );
        }
      }
    }
    $indel_bed_fh->close;
    $indel_review_fh->close;

    # If the user doesn't want to review Pindel calls, store them in a separate file
    if( $exclude_pindel )
    {
      my $pindel_review_fh = IO::File->new( $pindel_review_file, ">" ) or die "Cannot open $pindel_review_file. $!";
      $pindel_review_fh->print( "Chr\tStart\tStop\tRef\tVar\tCall\tNotes\n" );
      for my $chr ( nsort keys %{$review_lines{pindels}} )
      {
        for my $start ( sort {$a <=> $b} keys %{$review_lines{pindels}{$chr}} )
        {
          for my $stop ( sort {$a <=> $b} keys %{$review_lines{pindels}{$chr}{$start}} )
          {
            $pindel_review_fh->print( $review_lines{pindels}{$chr}{$start}{$stop}, "\n" );
          }
        }
      }
      $pindel_review_fh->close if( $exclude_pindel );
    }

    # Dump IGV XML files to make it easy on the manual reviewers
    my $reference = (( $refseq_version eq "GRCh37-lite-build37" ) ? "b37" : "reference" );
    unless( Genome::Model::Tools::Analysis::DumpIgvXml->execute( tumor_bam => $tumor_bam, normal_bam => $normal_bam, reference_name => $reference,
            review_bed_file => $indel_bed_file, review_description => "High-confidence Tier1 Indels", genome_name => "$case.indel", output_dir => $output_dir_new )->result and
            Genome::Model::Tools::Analysis::DumpIgvXml->execute( tumor_bam => $tumor_bam, normal_bam => $normal_bam, reference_name => $reference,
            review_bed_file => $snv_bed_file, review_description => "High-confidence Tier1 SNVs", genome_name => "$case.snv", output_dir => $output_dir_new )->result)
    {
      print STDERR "WARNING: Unable to generate IGV XMLs for $case\n";
    }
    print "done.\n";
  }

  return 1;
}

1;
