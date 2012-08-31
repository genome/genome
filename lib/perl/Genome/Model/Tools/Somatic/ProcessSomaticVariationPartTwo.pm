package Genome::Model::Tools::Somatic::ProcessSomaticVariationPartTwo;

use warnings;
use strict;
use IO::File;
use Genome;
use Sort::Naturally qw(nsort);
use Genome::Info::IUB;

class Genome::Model::Tools::Somatic::ProcessSomaticVariationPartTwo {
  is => 'Command',
  has_input => [
      somatic_variation_model_id => {
          is => 'Text',
          doc => "ID of SomaticVariation model",
      },
      
      output_dir => {
          is => 'Text',
          doc => "Directory where output will be stored (under a subdirectory with the sample name)",
      },

      ],
  
  has_optional_input => [
      add_readcounts =>{
          is => 'Boolean',
          is_optional => 1,
          default => 1,
          doc => "append readcounts to the output table",
      },

      add_tiers =>{
          is => 'Boolean',
          is_optional => 1,
          default => 0,
          doc => "add a tier column to the output table",
      },

      review_file =>{
          is => 'String',
          is_optional => 1,
          doc => "the script will try to grab this automatically, but if the review files are found in a non-standard place or have a non-standard name, specify the path here",
      },

      ],
};


sub help_detail {
  return <<HELP;
After gmt somatic process-somatic-variation is run and the variants are reviewed, this tool can be run to
process the review files and produce output (annotation format, mafs, etc)
HELP
}

sub _doc_authors {
  return <<AUTHS;
 Chris Miller
AUTHS
}


sub bedToAnno{
    my ($chr,$start,$stop,$ref,$var,@rest) = split("\t",$_[0]);
    if ($ref =~ /0|\-|\*/){ #indel INS
        $stop = $stop+1;
    } else { #indel DEL or SNV
        $start = $start+1;
    }
    return(join("\t",($chr,$start,$stop,$ref,$var,@rest)));
}


sub execute {
  my $self = shift;
  my $somatic_variation_model_id = $self->somatic_variation_model_id;
  my $output_dir = $self->output_dir;
  $output_dir =~ s/(\/)+$//; # Remove trailing forward-slashes if any

  # Check on the input data before starting work
  my $model = Genome::Model->get( $somatic_variation_model_id );
  print STDERR "ERROR: Could not find a model with ID: $somatic_variation_model_id\n" unless( defined $model );
  print STDERR "ERROR: Output directory not found: $output_dir\n" unless( -e $output_dir );
  return undef unless( defined $model && -e $output_dir );


  #grab the info from the model
  my $build = $model->last_succeeded_build;
  unless( defined($build) ){
      print STDERR "WARNING: Model ", $model->id, "has no succeeded builds\n";
      return undef;
  }

  my $tumor_model = $model->tumor_model;
  my $normal_model = $model->normal_model;

  my $ref_seq_build_id = $tumor_model->reference_sequence_build->build_id;
  my $ref_seq_build = Genome::Model::Build->get($ref_seq_build_id);
  my $ref_seq_fasta = $ref_seq_build->full_consensus_path('fa');
  my $sample_name = $tumor_model->subject->name;
  my $normal_sample_name = $normal_model->subject->name;
  print STDERR "processing model with sample_name: " . $sample_name . "\n";
  my $tumor_bam = $tumor_model->last_succeeded_build->whole_rmdup_bam_file;
  my $normal_bam = $normal_model->last_succeeded_build->whole_rmdup_bam_file;
  my $build_dir = $build->data_directory;

  my $annotation_build_id = $tumor_model->annotation_reference_build->id;

  
  #this will need to be fixed if applying to mouse 
  my $genome_build_number;
  if($ref_seq_build->name =~/36/){
      $genome_build_number = 36;
  }  elsif($ref_seq_build->name =~/37/){
      $genome_build_number = 37;
  }



  # Check if the necessary files exist
  my $review_file;
  if(defined($self->review_file)){
      $review_file = $self->review_file;
  } else {
      $review_file = glob("$output_dir/review/" . $sample_name . "_reviewed.*");
  }
  unless( -e $review_file ){
      die "ERROR: review file for $sample_name not found at $review_file. Specify the correct location using --review-file\n";
  }

  die "snv directory not found at $output_dir/$sample_name/snvs - did you specify the right output directory?" unless( -e "$output_dir/$sample_name/snvs" );
  die "indel directory not found at $output_dir/$sample_name/indels - did you specify the right output directory?" unless( -e "$output_dir/$sample_name/indels" );




  #clean up the review file, convert to one-based
  #--------------------------------------------------------------
  open(SNVREVFILE,">$output_dir/$sample_name/snvs/snvs.reviewed") || die ("couldn't open output file");
  open(INDELREVFILE,">$output_dir/$sample_name/indels/indels.reviewed") || die ("couldn't open output file");
  open(SNVSOMFILE,">$output_dir/$sample_name/snvs/snvs.reviewed.somatic") || die ("couldn't open output file");
  open(INDELSOMFILE,">$output_dir/$sample_name/indels/indels.reviewed.somatic") || die ("couldn't open output file");
  my $inFh = IO::File->new( $review_file ) || die "can't open review file $review_file\n";
  while( my $line = $inFh->getline )
  {
      chomp($line);
      $line =~ s/"//g; #spreadsheet often adds quotes
    
      #convert to one-based coords
      $line = bedToAnno($line);

      my ( $chr, $start, $stop, $ref, $var, $status, $notes ) = split( /\t/, $line );

      #notes are optional
      $notes = "" unless defined($notes);
      #status is not
      die "no review status on $line\n" unless defined($status);

#      next if $line =~/ambiguous but looks like normal contamination - restored/;

      if(($ref =~ /0|\-|\*/) || ($var =~ /0|\-|\*/)){ #indel
          print INDELREVFILE join("\t",( $chr, $start, $stop, $ref, $var, $status, $notes)) . "\n";
          if($status eq "S"){
              print INDELSOMFILE join("\t",( $chr, $start, $stop, $ref, $var)) . "\n";
          }
      } else { #snv
          print SNVREVFILE join("\t",( $chr, $start, $stop, $ref, $var, $status, $notes)) . "\n";
          if($status eq "S"){
              print SNVSOMFILE join("\t",( $chr, $start, $stop, $ref, $var)) . "\n";
          }
      }
  }
  close(SNVREVFILE);
  close(INDELREVFILE);
  close(INDELSOMFILE);

  #--------------
  #add in the sites that passed the uhc filter, if present
  my $uhc_file = glob("$output_dir/$sample_name/snvs/*.passuhc");
  if(defined($uhc_file)){
      if( -e $uhc_file){
          print STDERR "re-adding pass uhc file $uhc_file\n";
          $inFh = IO::File->new( $uhc_file ) || die "can't open file $uhc_file\n";
          while( my $line = $inFh->getline )
          {
              chomp($line);
              my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
              print SNVSOMFILE join("\t",( $chr, $start, $stop, $ref, $var)) . "\n";
          }
      }
  }

  close(SNVSOMFILE);
          
  #sort these files
  `joinx sort $output_dir/$sample_name/snvs/snvs.reviewed.somatic >$output_dir/$sample_name/zz && mv -f $output_dir/$sample_name/zz $output_dir/$sample_name/snvs/snvs.reviewed.somatic`;
  `joinx sort $output_dir/$sample_name/indels/indels.reviewed.somatic >$output_dir/$sample_name/zz && mv -f $output_dir/$sample_name/zz $output_dir/$sample_name/indels/indels.reviewed.somatic`;

  

  #------------
  #create annotation files
  print STDERR "annotating snvs and indels\n";
  if( -s "$output_dir/$sample_name/snvs/snvs.reviewed.somatic"){
      my $annoResult = Genome::Model::Tools::Annotate::TranscriptVariants->create(
          variant_file => "$output_dir/$sample_name/snvs/snvs.reviewed.somatic",
          output_file => "$output_dir/$sample_name/snvs/snvs.reviewed.somatic.anno",
          build_id => "$annotation_build_id",
          annotation_filter => "top",
          );
      unless ($annoResult->execute) {
          die "Failed to annotate somatic snvs";
      }
  } else {
      `touch $output_dir/$sample_name/snvs/snvs.reviewed.somatic.anno`;
  }

  if( -s "$output_dir/$sample_name/indels/indels.reviewed.somatic"){
      my $annoResult = Genome::Model::Tools::Annotate::TranscriptVariants->create(
          variant_file => "$output_dir/$sample_name/indels/indels.reviewed.somatic",
          output_file => "$output_dir/$sample_name/indels/indels.reviewed.somatic.anno",
          build_id => "$annotation_build_id",
          annotation_filter => "top",
          );
      unless ($annoResult->execute) {
          die "Failed to annotate somatic indels";
      }
  } else {
      `touch $output_dir/$sample_name/indels/indels.reviewed.somatic.anno`;
  }
  
 
  #--------------
  #create table, add tiers and readcounts if specified 
  my $tmpfile = "$output_dir/$sample_name/tmp";
  `cat $output_dir/$sample_name/snvs/snvs.reviewed.somatic.anno >$tmpfile`;
  `tail -n +2 $output_dir/$sample_name/indels/indels.reviewed.somatic.anno >>$tmpfile`;

  #skip if the tmpfile has zero size
  if( -s $tmpfile ){
      #add tiers
      if($self->add_tiers){
          print STDERR "adding tiers\n";
          my $tier_cmd  = Genome::Model::Tools::FastTier::AddTiers->create(
              input_file => $tmpfile,
              output_file => $tmpfile . ".tiered",
              build => $genome_build_number,
              );
          unless ($tier_cmd->execute) {
              die "Failed to tier variants.\n";
          }
          $tmpfile = $tmpfile . ".tiered";
      }
  
      #--------------
      if($self->add_readcounts){
          print STDERR "adding readcounts\n";
          #get readcounts from the normal bam
          my $normal_rc_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
              bam_file => $normal_bam,
              variant_file => $tmpfile,
              output_file => $tmpfile . ".nrm",
              genome_build => $ref_seq_fasta,
              header_prefix => "Normal",
              );
          unless ($normal_rc_cmd->execute) {
              die "Failed to obtain normal readcounts.\n";
          }
          
          #get readcounts from the tumor bam
          my $tumor_rc_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
              bam_file => $tumor_bam,
              variant_file => $tmpfile . ".nrm",
              output_file => $tmpfile . ".nrm.tum",
              genome_build => $ref_seq_fasta,
              header_prefix => "Tumor",
              );
          unless ($tumor_rc_cmd->execute) {
              die "Failed to obtain tumor readcounts.\n";
          }
          $tmpfile = $tmpfile . ".nrm.tum";
      }
  
      `mv $tmpfile $output_dir/$sample_name/somatic.variants.anno`;
      `rm -f $output_dir/$sample_name/tmp*`;
  } else {
      `touch $output_dir/$sample_name/somatic.variants.anno`;
  }


  #--------------
  print STDERR "creating maf file\n";
  #create maf file for this sample
  
  #have to tier things first so that we only include tier1 variants in the maf. grar.
  my $sfile = "$output_dir/$sample_name/snvs/snvs.reviewed.somatic";
  my $sannofile = "$output_dir/$sample_name/snvs/snvs.reviewed.somatic.anno";
  my $ifile = "$output_dir/$sample_name/indels/indels.reviewed.somatic";
  my $iannofile = "$output_dir/$sample_name/indels/indels.reviewed.somatic.anno";
  my @files = ($sfile, $sannofile, $ifile, $iannofile);

  foreach my $file (@files){
      if( -s $file){
          my $tier_cmd  = Genome::Model::Tools::FastTier::AddTiers->create(
              input_file => $file,
              output_file => $file . ".tiered",
              build => $genome_build_number,
              );
          unless ($tier_cmd->execute) {
              die "Failed to tier variants.\n";
          }  
          #grep out tier 1 variants from these files
          my $tfile = $file . ".tiered";
          `grep "tier1" $tfile >$file.tier1`;
      } else {
          `touch $file.tiered`;
          `touch $file.tier1`;
      }
  }


  
  my $maf_cmd = Genome::Model::Tools::Capture::CreateMafFile->create(
      output_file => "$output_dir/$sample_name/somatic.maf",
      genome_build => $genome_build_number,
      indel_file => "$output_dir/$sample_name/indels/indels.reviewed.somatic.tier1",
      indel_annotation_file => "$output_dir/$sample_name/indels/indels.reviewed.somatic.anno.tier1",
      snv_file => "$output_dir/$sample_name/snvs/snvs.reviewed.somatic.tier1",
      snv_annotation_file => "$output_dir/$sample_name/snvs/snvs.reviewed.somatic.anno.tier1",
      normal_sample => "$normal_sample_name",
      tumor_sample => "$sample_name",
      );
  unless ($maf_cmd->execute) {
      die "Failed to create maf file";
  }

  
  #----------------------
  # output some bam info we need for music
  open(OUTFILE,">$output_dir/$sample_name/bam_list");
  print OUTFILE join("\t",($sample_name,$normal_bam,$tumor_bam)) . "\n";
  close(OUTFILE);
      
  
return 1;
}

1;
