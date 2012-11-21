package Genome::Model::Tools::Validation::ProcessSomaticValidation;

use warnings;
use strict;
use IO::File;
use Genome;
use Sort::Naturally qw(nsort);
use Genome::Info::IUB;
use File::Basename;

class Genome::Model::Tools::Validation::ProcessSomaticValidation {
  is => 'Command',
  has_input => [
      somatic_validation_model_id => {
          is => 'Text',
          doc => "ID of SomaticValidation model",
      },

      output_dir => {
          is => 'Text',
          doc => "Directory where output will be stored (under a subdirectory with the sample name)",
      },

  ],

  has_optional_input => [

      variant_list => {
          is => 'Text',
          is_optional => 1,
          doc => "List of variants that we're trying to validate, in bed format. If not specified, it will try to grab this list from the model. If no list is found, it will assume that this is an extension experiment (and that there are no variants to validate)",
      },

      somatic_variation_model_id => {
          is => 'Text',
          is_optional => 1,
          doc => "Somatic variation model that was used to call the variants in the first place. Only used if --review-non-validated-calls is true. If not specified, we'll try to grab this from the som-val model link",
      },

      igv_reference_name => {
          is => 'Text',
          is_optional => 1,
          doc => "Name of the igv reference to use",
          default => "reference_build36",
      },

      tumor_only => {
          is => 'Boolean',
          is_optional => 1,
          default => 0,
          doc => "Model is an extension with only tumor data",
      },

      filter_sites => {
          is => 'Text',
          is_optional => 1,
          doc => "Comma separated list of bed files containing sites to be removed (example - removing cell-line sites from tumors grown on them). Only sites with these exact coordinates and that match the ref and var alleles are removed.",
      },

      filter_regions => {
          is => 'Text',
          is_optional => 1,
          doc => "Comma separated list of bed files containing regions to be removed (example - removing paralog regions) Any sites intersecting these regions are removed.",
      },

      restrict_to_target_regions => {
          is => 'Boolean',
          is_optional => 1,
          default => 1,
          doc => "Only keep snv calls within the target regions. These are pulled from the build",
      },

      tier1_only => {
          is => 'Boolean',
          is_optional => 1,
          default => 0,
          doc => "Only keep and review calls that are tier 1",
      },

      tier_file_location => {
          is => 'String',
          is_optional => 1,
          doc => "If tier1-only is specified, this needs to be a path to the appropriate tiering files",
      },

      recover_non_validated_calls => {
          is => 'Boolean',
          is_optional => 1,
          doc => "Do readcounting on calls that are non-validated, find sites that have high coverage (>20x tumor and normal) in the original, but poor (<8x tumor/ <6x norm) coverage in the validation. Send these for review",
          default => 0,
      },

      variant_list_is_one_based => {
          is => 'Boolean',
          is_optional => 1,
          doc => "The variant list you provide is in annotation (one-based) format, instead of bed. This flag fixes that.",
          default => 0,
      },

      dbsnp_filter => {
          is => 'String',
          is_optional => 1,
          doc => "Path to a dbsnp bed file of sites that should be removed from consideration",
      },

      ##restrict to targeted region - grab from build

      # read_review => {
      #     is => 'Boolean',
      #     doc => "Read existing manual review files and create WU annotation files per case",
      #     is_optional => 1,
      #     default => 0
      #   },

  ],
};

sub help_detail {
  return <<HELP;
Given a SomaticValidation model, this tool will gather the resulting variants, remove
off-target sites, tier the variants, optionally filter them, and match them up with the
initial predictions sent for validation. It will then divide them into categories (validated,
non-validated, and new calls). New calls are prepped for manual review in the review/ directory.
HELP
}

sub _doc_authors {
  return <<AUTHS;
 Chris Miller
AUTHS
}

sub bedToAnno {
    my ($chr,$start,$stop,$ref,$var) = split("\t",$_[0]);
    #print STDERR join("|",($chr,$start,$stop,$ref,$var)) . "\n";
    if ($ref =~ /^\-/) { #indel INS
        $stop = $stop+1;
    }
    else { #indel DEL or SNV
        $start = $start+1;
    }
    return(join("\t",($chr,$start,$stop,$ref,$var)));
}

sub annoToBed {
    my ($chr,$start,$stop,$ref,$var) = split("\t",$_[0]);
    if ($ref =~ /^\-/) { #indel INS
        $stop = $stop-1;
    }
    else { #indel DEL or SNV
        $start = $start-1;
    }
    return(join("\t",($chr,$start,$stop,$ref,$var)));
}

sub annoFileToSlashedBedFile {
    my $fh = shift;
    my $input_file = shift;

    my $inFh = IO::File->new( $input_file ) or die "Can't open file!";
    while( my $line = $inFh->getline ) {
        chomp($line);
        my $tmp = annoToBed($line);
        my @tmp2 = split("\t",$tmp);
        print $fh join("\t",(@tmp2[0..2], ($tmp2[3] . "/" . $tmp2[4]))) . "\n";
    }
    close($fh);
    close($inFh);
}

sub slashedBedFileToAnnoFile {
    my $fh = shift;
    my $input_file = shift;

    my $inFh = IO::File->new( $input_file ) or die "Can't open file!";
    while( my $line = $inFh->getline ) {
        chomp($line);
        my ( $chr, $start, $stop, $refvar, @rest ) = split( /\t/, $line );
        my @alleles = split("/",$refvar);
        my $line2 = bedToAnno(join("\t",($chr, $start, $stop, $alleles[0], $alleles[1]))) . "\n";
        print $fh $line2 . "\n";
    }
    close($fh);
    close($inFh);
}

sub intersects {
    my ($st,$sp,$st2,$sp2) = @_;
    if ((($sp2 >= $st) && ($sp2 <= $sp)) || (($sp >= $st2) && ($sp <= $sp2))) {
        return 1;
    }
    return 0;
}

sub fixIUB {
    my ($ref,$var) = @_;
    my @vars = Genome::Info::IUB->variant_alleles_for_iub($ref,$var);
    return @vars;
}

sub execute {
  my $self = shift;
  my $tumor_only = $self->tumor_only;
  my $somatic_validation_model_id = $self->somatic_validation_model_id;
  my $output_dir = $self->output_dir;
  $output_dir =~ s/(\/)+$//; # Remove trailing forward-slashes if any

  # Check on the input data before starting work
  my $model = Genome::Model->get( $somatic_validation_model_id );
  print STDERR "ERROR: Could not find a model with ID: $somatic_validation_model_id\n" unless( defined $model );
  print STDERR "ERROR: Output directory not found: $output_dir\n" unless( -e $output_dir );
  return undef unless( defined $model && -e $output_dir );

  #grab the info from all of the models
  my %bams; # Hash to store the model info
  my $build = $model->last_succeeded_build;
  unless( defined($build) ){
      print STDERR "WARNING: Model ", $model->id, "has no succeeded builds\n";
      return undef;
  }

  my $ref_seq_build_id = $model->reference_sequence_build->build_id;
  my $ref_seq_build = Genome::Model::Build->get($ref_seq_build_id);
  my $ref_seq_fasta = $ref_seq_build->full_consensus_path('fa');
  my $sample_name = $model->tumor_sample->name;
  $self->status_message( "Processing model with sample_name: " . $sample_name . "\n" );
  my $tumor_bam = $build->tumor_bam;
  my $build_dir = $build->data_directory;

  my $normal_bam;
  unless($tumor_only){
      $normal_bam = $build->normal_bam;
  }

  # Check if the necessary files exist in this build
  my $snv_file = "$build_dir/variants/snvs.hq.bed";
  unless( -e $snv_file ){
      die "ERROR: SNV results annotations for $sample_name not found at $snv_file\n";
  }
  my $indel_file = "$build_dir/variants/indels.hq.bed";
  unless( -e $indel_file ){
      die "ERROR: INDEL results annotations for $sample_name not found at $indel_file\n";
  }

  # create subdirectories, get files in place
  mkdir "$output_dir/$sample_name" unless( -e "$output_dir/$sample_name" );
  mkdir "$output_dir/$sample_name/snvs" unless( -e "$output_dir/$sample_name/snvs" );
  mkdir "$output_dir/$sample_name/indels" unless( -e "$output_dir/$sample_name/indels" );
  mkdir "$output_dir/review" unless( -e "$output_dir/review" );
  `ln -s $build_dir $output_dir/$sample_name/build_directory`;
  `ln -s $snv_file $output_dir/$sample_name/snvs/` unless( -e "$output_dir/$sample_name/snvs/$snv_file");
  `ln -s $indel_file $output_dir/$sample_name/indels/` unless( -e "$output_dir/$sample_name/indels/$indel_file");
  $snv_file = "$output_dir/$sample_name/snvs/snvs.hq.bed";
  $indel_file = "$output_dir/$sample_name/indels/indels.hq.bed";

  my %dups;
  #munge through SNV file to remove duplicates and fix IUB codes
  #--------------------------------------------------------------
  my $filenm = $snv_file;
  $filenm =~ s/.bed/.clean.bed/g;

  open(SNVFILE,">$filenm") or die "Can't open filter file!";
  my $inFh = IO::File->new( $snv_file ) or die "Can't open file!";
  while( my $line = $inFh->getline ) {
      chomp($line);
      my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
      if($ref =~ /\//){
          ( $ref, $var ) = split(/\//, $ref);
      }
      my @vars = fixIUB($ref, $var);
      foreach my $v (@vars){
          unless (exists($dups{join("\t",($chr, $start, $stop, $ref, $v ))})){
              print SNVFILE join("\t",($chr, $start, $stop, $ref, $v )) . "\n";
          }
          $dups{join("\t",($chr, $start, $stop, $ref, $v ))} = 1;
      }
  }
  close(SNVFILE);
  $snv_file = $filenm;

  #clean up the indel file too
  $filenm = $indel_file;
  $filenm =~ s/.bed/.clean.bed/g;

  open(INDELFILE,">$filenm") or die "Can't open filter file!";
  $inFh = IO::File->new( $indel_file ) or die "Can't open file!";
  while( my $line = $inFh->getline ) {
      chomp($line);
      my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
      if($ref =~ /\//){
          ( $ref, $var ) = split(/\//, $ref);
      }
      $ref =~ s/0/-/g;
      $var =~ s/0/-/g;
      $ref =~ s/\*/-/g;
      $var =~ s/\*/-/g;
      unless (exists($dups{join("\t",($chr, $start, $stop, $ref, $var ))})){
          print INDELFILE join("\t",($chr, $start, $stop, $ref, $var )) . "\n";
      }
      $dups{join("\t",($chr, $start, $stop, $ref, $var ))} = 1;
  }
  close(INDELFILE);
  $indel_file = $filenm;

  #-------------------------------------------------
  #store the previously called variants into a hash
  my %prevCalls;

  my $variant_list = "";
  if(defined($self->variant_list)){
      $variant_list = $self->variant_list;
  }
  else { #look in model
      if(defined($model->snv_variant_list)){
          $variant_list = $model->snv_variant_list->original_file_path
      }
  }
  if(defined($variant_list)){
      print "Using variant list: $variant_list\n";
  }

  if( -e $variant_list){
      my $inFh = IO::File->new( $variant_list ) or die "Can't open file!";
      while( my $line = $inFh->getline )
      {
          chomp($line);
          #handle either 5 col (Ref\tVar) or 4 col (Ref/Var) bed
          my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
          if($ref =~ /\//){
              ( $ref, $var ) = split(/\//, $ref);
          }
          $ref =~ s/0/-/g;
          $var =~ s/0/-/g;
          $ref =~ s/\*/-/g;
          $var =~ s/\*/-/g;

          if(($ref eq "-") || ($var eq "-")){
              my $key = join("\t",($chr, $start, $stop, $ref, $var ));
              if($self->variant_list_is_one_based){
                  $prevCalls{annoToBed($key)} = 0;
              }
              else {
                  $prevCalls{$key} = 0;
              }
          }
          else {
              my @vars = fixIUB($ref, $var);
              foreach my $v (@vars){
                  my $key = join("\t",($chr, $start, $stop, $ref, $v ));
                  if($self->variant_list_is_one_based){
                      $prevCalls{annoToBed($key)} = 0;
                  }
                  else {
                      $prevCalls{$key} = 0;
                  }
              }
          }
      }
      close($inFh);

  }
  else {
      $self->status_message( "WARNING: BED file of targeted variants not found for $sample_name.\nAssuming that there were no calls for this model (and it was an extension experiment)\n" );
  }

  #--------------------------------------------------------
  #output the (non)validated calls before doing any filtering
   # Grab the high confidence calls from their respective files
   #---first snvs---
  open(VALFILE,">$output_dir/$sample_name/snvs/snvs.validated");
  open(NEWFILE,">$output_dir/$sample_name/snvs/snvs.newcalls");
  open(MISFILE,">$output_dir/$sample_name/snvs/snvs.notvalidated");

  $inFh = IO::File->new( $snv_file ) or die "can't open file!";
  while( my $line = $inFh->getline )
  {
      chomp($line);
      my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
      if($ref =~ /\//){
          ( $ref, $var ) = split(/\//, $ref);
      }

      if (exists($prevCalls{join("\t",($chr, $start, $stop, $ref, $var ))})){
          print VALFILE bedToAnno(join("\t",($chr, $start, $stop, $ref, $var ))) . "\n";
          $prevCalls{join("\t",($chr, $start, $stop, $ref, $var ))} = 1;
      }
      else {
          print NEWFILE bedToAnno(join("\t",($chr, $start, $stop, $ref, $var ))) . "\n";
      }
  }

  close($inFh);
  #case 3: called in original, not in validation
  foreach my $k (keys(%prevCalls)){
      next if $prevCalls{$k} == 1;
      my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $k );
      #skip indels
      unless (($ref =~ /-|0/) || ($var =~ /-|0/)){
          print MISFILE bedToAnno(join("\t",($chr, $start, $stop, $ref, $var ))) . "\n";
      }
  }

  close(VALFILE);
  close(NEWFILE);
  close(MISFILE);

  #---now indels---
  open(VALFILE,">$output_dir/$sample_name/indels/indels.validated");
  open(NEWFILE,">$output_dir/$sample_name/indels/indels.newcalls");
  open(MISFILE,">$output_dir/$sample_name/indels/indels.notvalidated");
  $inFh = IO::File->new( $indel_file ) or die "Can't open file!";
  while( my $line = $inFh->getline )
  {
      chomp($line);
      my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
      if($ref =~ /\//){
          ( $ref, $var ) = split(/\//, $ref);
      }
      $ref =~ s/0/-/g;
      $var =~ s/0/-/g;
      #case 1 - previously found variant, now validated
      if (exists($prevCalls{join("\t",($chr, $start, $stop, $ref, $var ))})){
          print VALFILE bedToAnno(join("\t",($chr, $start, $stop, $ref, $var ))) . "\n";
          #mark as found
          $prevCalls{join("\t",($chr, $start, $stop, $ref, $var ))} = 1;
      }
      else { #case 2: new indel not found in original targets
          print NEWFILE bedToAnno(join("\t",($chr, $start, $stop, $ref, $var ))) . "\n";
      }
  }
  #case 3: called in original, not in validation
  foreach my $k (keys(%prevCalls)){
      next if $prevCalls{$k} == 1;

      my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $k );
      #only indels
      if (($ref =~ /-/) || ($var =~ /-/)){
          print MISFILE bedToAnno(join("\t",($chr, $start, $stop, $ref, $var ))) . "\n";
      }
  }
  close(VALFILE);
  close(NEWFILE);
  close(MISFILE);

  $indel_file = "$output_dir/$sample_name/indels/indels.newcalls";
  $snv_file = "$output_dir/$sample_name/snvs/snvs.newcalls";

  #-------------------------------------------------
  #filter out the off-target regions, if target regions are available
  if ($self->restrict_to_target_regions) {
      $self->status_message( "Filtering out off-target regions...\n" );

      #Locate the BED file of targeted regions used in the SomVal model
      my ( $featurelist_name, $featurelist_id ) = ( $model->target_region_set->name, $model->target_region_set->id );
      $self->status_message( "Using target regions in feature-list $featurelist_name (ID: $featurelist_id)\n" );
      my $featurelist_file = Genome::FeatureList->get($featurelist_id)->file_path;

      if ( -e $featurelist_file ) {
          my %targetRegions;
          my $inFh = IO::File->new( $featurelist_file ) or die "Can't open file!";
          while( my $line = $inFh->getline ) {
              chomp($line);
              next if $line =~ /^track name/;
              my ($chr, $start, $stop, @rest) = split( /\t/, $line );
              #remove chr if present
              $chr =~ s/^chr//g;
              $targetRegions{$chr}{join("\t",($start, $stop))} = 0;
          }
          close($inFh);

          #compare the snvs to the targets
          open(TARFILE,">$snv_file.ontarget") or die "Can't open target file!";
          $inFh = IO::File->new( $snv_file ) or die "Can't open file!";
          while( my $line = $inFh->getline ) {
              chomp($line);
              my ( $chr, $start, $stop, @rest ) = split( /\t/, $line );

              #if we run into huge lists, this will be slow - refactor to use joinx - TODO
              my $found = 0;
              foreach my $pos (keys(%{$targetRegions{$chr}})) {
                  my ($tst, $tsp) = split("\t",$pos);
                  if(intersects($start,$stop,$tst,$tsp)) {
                      $found = 1;
                  }
              }
              if($found) {
                  print TARFILE $line . "\n";
              }
          }
          close($inFh);
          close(TARFILE);
          $snv_file = "$snv_file.ontarget";

          #compare the indels to the targets
          open(TARFILE,">$indel_file.ontarget") or die "Can't open target file!";
          $inFh = IO::File->new( $indel_file ) or die "Can't open file!";
          while( my $line = $inFh->getline ) {
              chomp($line);
              my ( $chr, $start, $stop, @rest ) = split( /\t/, $line );
              foreach my $pos (keys(%{$targetRegions{$chr}})) {
                  my ($tst, $tsp) = split("\t",$pos);
                  if(intersects($start,$stop,$tst,$tsp)) {
                      print TARFILE $line . "\n";
                  }
              }
          }
          close($inFh);
          close(TARFILE);
          $indel_file = "$indel_file.ontarget";
      }
      else {
          $self->status_message( "WARNING: Feature list not found at $featurelist_file\nNo target region filtering being done\n" );
      }
  }

  ##------------------------------------------------------
  #remove all but tier 1 sites, if that option is specified
  if($self->tier1_only) {
  $self->status_message( "Doing Tiering...\n" );
      my $tier_cmd = Genome::Model::Tools::FastTier::FastTier->create(
          tier_file_location => $self->tier_file_location,
          variant_bed_file => $snv_file,
      );
      unless ($tier_cmd->execute) {
          die "Failed to tier variants successfully.\n";
      }
      $snv_file = "$snv_file.tier1";

      $tier_cmd = Genome::Model::Tools::FastTier::FastTier->create(
          tier_file_location => $self->tier_file_location,
          variant_bed_file => $indel_file,
      );
      unless ($tier_cmd->execute) {
          die "Failed to tier variants successfully.\n";
      }
      $indel_file = "$indel_file.tier1";
  }

  ##-------------------------------------------------
  #use joinx to remove dbsnp sites
  if(defined($self->dbsnp_filter)) {
      $self->status_message( "Applying dbsnp filter...\n" );

      #snvs - convert back to slashed bed file before joinx intersect
      my ($fh,$temp_bed_file) = Genome::Sys->create_temp_file;
      annoFileToSlashedBedFile($fh, $snv_file);
      close($fh);

      my ($fh2,$temp_novel_file) = Genome::Sys->create_temp_file;
      my $cmd = "joinx intersect --dbsnp-match --miss-a $temp_novel_file -a $temp_bed_file -b " . $self->dbsnp_filter . ">/dev/null";
      my $result = Genome::Sys->shellcmd( cmd => "$cmd" );
      unless($result) {
          $self->error_message("Failed to execute joinx: Returned $result");
          die $self->error_message;
      }

      my $ofh;
      open($ofh,">$snv_file.novel");
      slashedBedFileToAnnoFile($ofh, $temp_novel_file);
      close($ofh);
      $snv_file = "$snv_file.novel";

      #indels - convert back to slashed bed file before joinx intersect
      my ($fh3,$temp_bed_file2) = Genome::Sys->create_temp_file;
      annoFileToSlashedBedFile($fh3, $indel_file);
      close($fh3);

      my ($fh4,$temp_novel_file2) = Genome::Sys->create_temp_file;
      $cmd = "joinx intersect --dbsnp-match --miss-a $temp_novel_file2 -a $temp_bed_file2 -b " . $self->dbsnp_filter . ">/dev/null";
      $result = Genome::Sys->shellcmd( cmd => "$cmd" );
      unless($result) {
          $self->error_message("Failed to execute joinx: Returned $result");
          die $self->error_message;
      }

      open($ofh,">$indel_file.novel");
      slashedBedFileToAnnoFile($ofh, $temp_novel_file);
      close($ofh);
      # $inFh = IO::File->new( "$indel_file.novel" );
      # slashedBedFileToAnnoFile($inFh, $temp_novel_file2);
      # close($inFh);
      $indel_file = "$indel_file.novel";

  }

  #-------------------------------------------------
  #remove filter sites specified by the user
  if(defined($self->filter_sites)) {
      $self->status_message( "Applying user-supplied filter...\n" );
      my @filters = split(",",$self->filter_sites);
      my %filterSites;

      foreach my $filterfile (@filters) {
          if( -e $filterfile) {
              #store sites to filter out in a hash
              my $inFh = IO::File->new( $filterfile ) or die "Can't open file!";
              while( my $line = $inFh->getline ) {
                  chomp($line);
                  #handle either 5 col (Ref\tVar) or 4 col (Ref/Var) bed
                  my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
                  if($ref =~ /\//){
                      ( $ref, $var ) = split(/\//, $ref);
                  }
                  $ref =~ s/0/-/g;
                  $var =~ s/0/-/g;

                  my @vars = fixIUB($ref, $var);
                  foreach my $v (@vars){
                      $filterSites{join("\t",($chr, $start, $stop, $ref, $v ))} = 0;
                  }
              }
              close($inFh);
          }
          else {
              die("Filter sites file does not exist: " . $filterfile);
          }

          #remove snvs
          open(FILFILE,">$snv_file.filtered") or die "Can't open filter file!";
          $inFh = IO::File->new( $snv_file ) or die "Can't open file!";
          while( my $line = $inFh->getline ) {
              chomp($line);
              my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
              if($ref =~ /\//){
                  ( $ref, $var ) = split(/\//, $ref);
              }
              unless (exists($filterSites{join("\t",($chr, $start, $stop, $ref, $var ))})){
                  print FILFILE $line . "\n";
              }
          }
          close(FILFILE);
          $snv_file = "$snv_file.filtered";

          #remove indels
          open(FILFILE,">$indel_file.filtered") or die "Can't open filter file!";
          $inFh = IO::File->new( $indel_file ) or die "Can't open file!";
          while( my $line = $inFh->getline ) {
              chomp($line);
              my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
              if($ref =~ /\//){
                  ( $ref, $var ) = split(/\//, $ref);
              }

              unless (exists($filterSites{join("\t",($chr, $start, $stop, $ref, $var ))})){
                  print FILFILE $line . "\n";
              }
          }
          close(FILFILE);
          $indel_file = "$indel_file.filtered";
      }
  }

  #-------------------------------------------------
  #remove filter regions specified by the user
  if(defined($self->filter_regions)){
      #cat all the filter files together into one bed file
      my @filters = split(",",$self->filter_regions);

      my ($fh,$temp_file) = Genome::Sys->create_temp_file;

      foreach my $filterfile (@filters){
          my $inFh2 = IO::File->new( $filterfile ) or die "Can't open file!";
          while( my $line = $inFh2->getline ) {
              print $fh $line;
          }
          close $inFh2;
      }
      close($fh);

      $self->status_message( "Removing user-specified filter...\n" );
      my $cmd = "joinx intersect --miss-a $snv_file.filteredReg -a $snv_file -b $temp_file >/dev/null";
      my $result = Genome::Sys->shellcmd( cmd => "$cmd" );
      unless($result) {
          $self->error_message("Failed to execute joinx: Returned $result");
          die $self->error_message;
      }
      $snv_file = "$snv_file.filteredReg";

      $cmd = "joinx intersect --miss-a $indel_file.filteredReg -a $indel_file -b $temp_file >/dev/null";
      $result = Genome::Sys->shellcmd( cmd => "$cmd" );
      unless($result) {
          $self->error_message("Failed to execute joinx: Returned $result");
          die $self->error_message;
      }
      $indel_file = "$indel_file.filteredReg";
  }

  #add snv readcounts
  $self->status_message( "Getting snv readcounts...\n" );
  mkdir "$output_dir/$sample_name/snvs/readcounts";
  foreach my $file ("snvs.validated","snvs.notvalidated",basename($snv_file)){

      my $dir = "$output_dir/$sample_name/snvs/";
      if( -s "$dir/$file" ){
          unless($tumor_only){
              #get readcounts from the normal bam
              my $normal_rc_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
                  bam_file => $normal_bam,
                  output_file => "$dir/readcounts/$file.nrm.cnt",
                  variant_file => "$dir/$file",
                  genome_build => $ref_seq_fasta,
                  );
              unless ($normal_rc_cmd->execute) {
                  die "Failed to obtain normal readcounts for file $file.\n";
              }
          }

          #get readcounts from the tumor bam
          my $tumor_rc_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
              bam_file => $tumor_bam,
              output_file => "$dir/readcounts/$file.tum.cnt",
              variant_file => "$dir/$file",
              genome_build => $ref_seq_fasta,
              );
          unless ($tumor_rc_cmd->execute) {
              die "Failed to obtain tumor readcounts for file $file.\n";
          }
      }
  }

  #add indel readcounts
  $self->status_message( "Getting indel readcounts...\n" );
  mkdir "$output_dir/$sample_name/indels/readcounts";
  foreach my $file ("indels.validated","indels.notvalidated",basename($indel_file)){

      my $dir = "$output_dir/$sample_name/indels/";
      if( -s "$dir/$file" ){
          unless($tumor_only){
              #get readcounts from the normal bam
              my $normal_rc_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
                  bam_file => $normal_bam,
                  output_file => "$dir/readcounts/$file.nrm.cnt",
                  variant_file => "$dir/$file",
                  genome_build => $ref_seq_fasta,
                  indel_size_limit => 4,
                  );
              unless ($normal_rc_cmd->execute) {
                  die "Failed to obtain normal readcounts for file $file.\n";
              }
          }

          #get readcounts from the tumor bam
          my $tumor_rc_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
              bam_file => $tumor_bam,
              output_file => "$dir/readcounts/$file.tum.cnt",
              variant_file => "$dir/$file",
              genome_build => $ref_seq_fasta,
              indel_size_limit => 4,
              );
          unless ($tumor_rc_cmd->execute) {
              die "Failed to obtain tumor readcounts for file $file.\n";
          }
      }
  }

  ## we're not going to do this in most cases - if not called in validation, we just ditch the call,
  ## but it may be useful in some cases.

  #-------------------------------------------------
  #look at the calls that were missed (case 3 - called in original, failed validation)
  #to determine whether they were missed due to poor coverage.
  #if coverage is fine, dump them (most), but if coverage was poor in validation (and good in wgs), send for review
  #we can only really do this for snvs at the moment, until indel readcounting is tweaked
  my $somvarid;
  if(defined($self->somatic_variation_model_id)){
      $somvarid = $self->somatic_variation_model_id;
  }
  else {
      if(defined($model->snv_variant_list)){
          if(defined($model->snv_variant_list->source_build)){
              $somvarid = $model->snv_variant_list->source_build->model->id;
          }
      }
  }

  my $som_var_model;
  if(defined($somvarid)){
      my $som_var_model = Genome::Model->get($somvarid);
  }

  if( -s "$output_dir/$sample_name/snvs/snvs.notvalidated" ){
      if($self->recover_non_validated_calls){

          if(!(defined($somvarid))){
              print STDERR "ERROR: No somatic variation build found - can't recover non-validated calls.\n";
          }
          else {
              #get original bams
              my $tumor_bam_var = $som_var_model->tumor_model->last_succeeded_build->whole_rmdup_bam_file;
              print STDERR "ERROR: Somatic Variation tumor bam not found\n" unless( -e $tumor_bam_var );

              my $normal_bam_var = $som_var_model->normal_model->last_succeeded_build->whole_rmdup_bam_file;
              print STDERR "ERROR: Somatic Variation normal bam not found\n" unless( -e $normal_bam_var );

              if(defined($normal_bam_var)){
                  #get readcounts from the original normal sample
                  my $normal_rc2_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
                      bam_file => $normal_bam_var,
                      output_file => "$output_dir/$sample_name/snvs/readcounts/snvs.notvalidated.orig.nrm.cnt",
                      variant_file => "$output_dir/$sample_name/snvs/snvs.notvalidated",
                      genome_build => $ref_seq_fasta,
                      );
                  unless ($normal_rc2_cmd->execute) {
                      die "Failed to obtain normal readcounts.\n";
                  }
              }

              if(defined($tumor_bam_var)){
                  #get readcounts from the original tumor sample
                  my $tumor_rc2_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
                      bam_file => $tumor_bam_var,
                      output_file => "$output_dir/$sample_name/snvs/readcounts/snvs.notvalidated.orig.tum.cnt",
                      variant_file => "$output_dir/$sample_name/snvs/snvs.notvalidated",
                      genome_build => $ref_seq_fasta,
                      );
                  unless ($tumor_rc2_cmd->execute) {
                      die "Failed to obtain tumor readcounts.\n";
                  }
              }

              #read in all the validation readcounts, keep only those with poor coverage
              #require 8 reads in tumor, 6 in normal, per varscan cutoffs

              my %poorly_covered_snvs;

              unless($tumor_only){
                  my $inFh = IO::File->new( "$output_dir/$sample_name/snvs/readcounts/snvs.notvalidated.nrm.cnt" ) or die "Can't open file!";
                  while( my $line = $inFh->getline ) {
                      chomp($line);
                      my ( $chr, $start, $ref, $var, $refcnt, $varcnt, $vaf) = split("\t",$line);
                      if(($refcnt+$varcnt) < 6){
                          $poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))} = 0;
                      }
                  }
                  close($inFh);
              }

              $inFh = IO::File->new( "$output_dir/$sample_name/snvs/readcounts/snvs.notvalidated.tum.cnt" ) or die "Can't open file!";
              while( my $line = $inFh->getline ) {
                  chomp($line);
                  my ( $chr, $start, $ref, $var, $refcnt, $varcnt, $vaf) = split("\t",$line);
                  if(($refcnt+$varcnt) < 8){
                      $poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))} = 0;
                  }
              }
              close($inFh);

              #now, go through the original readcounts, and flag any that do have good coverage for manual review
              $inFh = IO::File->new( "$output_dir/$sample_name/snvs/readcounts/snvs.notvalidated.orig.nrm.cnt" ) or die "Can't open file!";
              while( my $line = $inFh->getline ) {
                  chomp($line);
                  my ( $chr, $start, $ref, $var, $refcnt, $varcnt, $vaf) = split("\t",$line);

                  if(defined($poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))})){
                      if(($refcnt+$varcnt) >= 20){
                          $poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))} = 1;
                      }
                  }
              }
              close($inFh);

              my $poorCount = 0;
              my @poorOutput;
              $inFh = IO::File->new( "$output_dir/$sample_name/snvs/readcounts/snvs.notvalidated.orig.tum.cnt" ) or die "Can't open file!";
              while( my $line = $inFh->getline ) {
                  chomp($line);
                  my ( $chr, $start, $ref, $var, $refcnt, $varcnt, $vaf) = split("\t",$line);
                  if(($refcnt+$varcnt) >= 20){
                      if(defined($poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))})){
                          if ($poorly_covered_snvs{join(":",( $chr, $start, $ref, $var ))} == 1){
                              push(@poorOutput,join("\t",( $chr, $start-1, $start, $ref, $var )));
                              $poorCount++;
                          }
                      }
                  }
              }
              close($inFh);

              if(@poorOutput){
                  open(OUTFILE,">$output_dir/review/$sample_name.poorValCoverage.bed");
                  foreach my $site (@poorOutput){
                      print OUTFILE $site . "\n";
                  }
                  close(OUTFILE)
              }

              if($poorCount > 0){
                  #create the xml file for this 4-way review
                  my $dumpCovXML = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
                      bams => join(",",($normal_bam,$tumor_bam,$normal_bam_var,$tumor_bam_var)),
                      labels => join(",",("validation normal $sample_name","validation tumor $sample_name","original normal $sample_name","original tumor $sample_name")),
                      output_file => "$output_dir/review/$sample_name.poorValCoverage.xml",
                      genome_name => $sample_name,
                      review_bed_file => "$output_dir/review/$sample_name.poorValCoverage.bed",
                      reference_name => $self->igv_reference_name,
                      );
                  unless ($dumpCovXML->execute) {
                      die "Failed to dump IGV xml for poorly covered sites.\n";
                  }

                  $self->status_message( "--------------------------------------------------------------------------------\n" .
                      "Sites for review with poor coverage in validation but good coverage in original are here:\n" .
                      "$output_dir/review/$sample_name.poorValCoverage.bed\n" .
                      "IGV XML file is here:\n" .
                      "$output_dir/$sample_name/review/poorValCoverage.xml\n\n"
                  );
              }
          }
      }
  }

  #-------------------------------------------------
  # look at the new calls that were found in validation, but not in the first build (Case 2 above)
  # in the case of extension experiments, with no previous build, this will be all variants
  print "Gathering new sites...\n";

  if ( -s "$snv_file"){
      #can't run UHC if tumor-only:
      unless($tumor_only){
          print "Running UHC filter...\n";
          #run the uhc filter to remove solid calls
          my $uhc_cmd = Genome::Model::Tools::Somatic::UltraHighConfidence->create(
              normal_bam_file => $normal_bam,
              tumor_bam_file => $tumor_bam,
              output_file => "$snv_file.passuhc",
              variant_file => "$snv_file",
              reference => $ref_seq_fasta,
              filtered_file => "$snv_file.failuhc",
              );
          unless ($uhc_cmd->execute) {
              die "Failed to run UHC filter.\n";
          }
      }


      #now get the files together for review
      print "Generating Review files...\n";
      my $revfile;
      if ( -s "$snv_file.failuhc"){
          $revfile = "$snv_file.failuhc";
      }
      else {
          $revfile = "$snv_file";
      }

      open(OUTFILE2,">$output_dir/review/$sample_name.newcalls.bed") or die "Can't open outfile!";

      $inFh = IO::File->new( $revfile ) or die "Can't open file!";
      while( my $line = $inFh->getline ) {
          chomp($line);
          my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
          print OUTFILE2 annoToBed(join("\t",($chr, $start, $stop, $ref, $var ))) . "\n";
      }
      close(OUTFILE2);

      my $bam_files;
      my $labels;
      if($tumor_only){
          if(defined($somvarid)){
              #add tumor and normal from somatic-variation model
              my $som_var_model = Genome::Model->get( $somvarid );
              my $tbam = $som_var_model->tumor_model->last_succeeded_build->whole_rmdup_bam_file;
              my $nbam = $som_var_model->normal_model->last_succeeded_build->whole_rmdup_bam_file;
              $bam_files = join(",",($tumor_bam,$tbam,$nbam));
              $labels = join(",",("validation tumor $sample_name","original tumor $sample_name","original normal $sample_name"));
          }
          else {
              $bam_files = join(",",($tumor_bam));
              $labels = join(",",("validation tumor $sample_name"));
          }
      }
      else {
          $bam_files = join(",",($normal_bam,$tumor_bam));
          $labels = join(",",("validation normal $sample_name","validation tumor $sample_name"));
      }

      #create the xml file for review
      my $dumpXML = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
          bams => "$bam_files",
          labels => "$labels",
          output_file => "$output_dir/review/$sample_name.newcalls.xml",
          genome_name => $sample_name,
          review_bed_file => "$output_dir/review/$sample_name.newcalls.bed",
          reference_name => $self->igv_reference_name,
          );
      unless ($dumpXML->execute) {
          die "Failed to dump IGV xml for poorly covered sites.\n";
      }

      $self->status_message( "\n--------------------------------------------------------------------------------\n" .
          "Sites to review that were not found original genomes, but were found in validation are here:\n" .
          "$output_dir/review/$sample_name.newcalls.bed\n" .
          "IGV XML file is here:\n" .
          "$output_dir/review/$sample_name.newcalls.xml\n\n"
      );
  }
  else {
      print STDERR "No variants found that were called in the validation, but not found in original genomes\n";
  }

  return 1;
}

1;
