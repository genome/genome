package Genome::Model::Tools::Somatic::ProcessSomaticVariation;

use warnings;
use strict;
use IO::File;
use Genome;
use Sort::Naturally qw(nsort);
use Genome::Info::IUB;

class Genome::Model::Tools::Somatic::ProcessSomaticVariation {
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

      igv_reference_name =>{
          is => 'Text',
          is_optional => 1,
          doc => "name of the igv reference to use",
          example_values => ["reference_build36"],
      },

      filter_sites =>{
          is => 'Text',
          is_optional => 1,
          doc => "comma separated list of bed files containing sites to be removed (example - removing cell-line sites from tumors grown on them). Only sites with these exact coordinates and that match the ref and var alleles are removed.",
      },

      filter_regions =>{
          is => 'Text',
          is_optional => 1,
          doc => "comma separated list of bed files containing regions to be removed (example - removing paralog regions) Any sites intersecting these regions are removed.",
      },

      get_readcounts =>{
          is => 'Boolean',
          is_optional => 1,
          default => 1,
          doc => "create readcount files for the final variant lists",
      },

      restrict_to_target_regions =>{
          is => 'Boolean',
          is_optional => 1,
          default => 1,
          doc => "only keep snv calls within the target regions. These are pulled from the build if possible",
      },

      tier1_only =>{
          is => 'Boolean',
          is_optional => 1,
          default => 0,
          doc => "only keep and review calls that are tier 1",
      },

      tier123_only =>{
          is => 'Boolean',
          is_optional => 1,
          default => 0,
          doc => "only keep and review calls that are tier 1, 2, or 3",
      },

      tier_file_location =>{
          is => 'String',
          is_optional => 1,
          doc => "if tier1-only is specified, this needs to be a path to the appropriate tiering files",
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
          doc => "path to a dbsnp bed file of sites that should be removed from consideration",
      },

      process_svs => {
          is => 'Boolean',
          is_optional => 1,
          doc => "Build has sv calls (probably WGS data). Most exomes won't have this",
          default => 0,
      },
      
      sites_to_pass => {
          is => 'String',
          is_optional => 1,
          doc => "an annotation file (5 col, 1 based) containing sites that will automatically be passed. This is useful when sequencing a relapse - the sites already found in the tumor don't need to be manually reviewed",
      },
      
      ],
};


sub help_detail {
  return <<HELP;
Given a SomaticVariation model, this tool will gather the resulting variants, remove
off-target sites, tier the variants, optionally filter them, etc. Calls are prepped for 
manual review in the review/ directory.
HELP
}

sub _doc_authors {
  return <<AUTHS;
 Chris Miller
AUTHS
}

sub bedToAnno{
    my ($chr,$start,$stop,$ref,$var) = split("\t",$_[0]);
    #print STDERR join("|",($chr,$start,$stop,$ref,$var)) . "\n";
    if ($ref =~ /^\-/){ #indel INS
        $stop = $stop+1;
    } else { #indel DEL or SNV
        $start = $start+1;
    }
    return(join("\t",($chr,$start,$stop,$ref,$var)));
}

sub annoToBed{
    my ($chr,$start,$stop,$ref,$var) = split("\t",$_[0]);
    if ($ref =~ /^\-/){ #indel INS
        $stop = $stop-1;
    } else { #indel DEL or SNV
        $start = $start-1;
    }
    #handle 5 col or 4 col ref/var
    if(defined($var)){
        return(join("\t",($chr,$start,$stop,$ref,$var)));
    } else {
        return(join("\t",($chr,$start,$stop,$ref)));
    }

}

sub annoFileToSlashedBedFile{
    my $fh = shift;
    my $input_file = shift;
    
    my $inFh = IO::File->new( $input_file ) || die "can't open file1\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my $tmp = annoToBed($line);
        my @tmp2 = split("\t",$tmp);
        print $fh join("\t",(@tmp2[0..2], ($tmp2[3] . "/" . $tmp2[4]))) . "\n"; 
    }
    close($fh);
    close($inFh);
}

sub slashedBedFileToAnnoFile{
    my $fh = shift;
    my $input_file = shift;
    
    my $inFh = IO::File->new( $input_file ) || die "can't open file2\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my ( $chr, $start, $stop, $refvar, @rest ) = split( /\t/, $line );
        my @alleles = split("/",$refvar);
        my $line2 = bedToAnno(join("\t",($chr, $start, $stop, $alleles[0], $alleles[1]))) . "\n";
        print $fh $line2 . "\n";
            
    }
    close($fh);
    close($inFh);
}


sub intersects{
    my ($st,$sp,$st2,$sp2) = @_;
    if((($sp2 >= $st) && ($sp2 <= $sp)) ||
       (($sp >= $st2) && ($sp <= $sp2))){
        return 1;
    }
    return 0;
}

sub fixIUB{
    my ($ref,$var) = @_;
    my @vars = Genome::Info::IUB->variant_alleles_for_iub($ref,$var);
    return @vars;
}


sub execute {
  my $self = shift;
  my $somatic_variation_model_id = $self->somatic_variation_model_id;
  my $output_dir = $self->output_dir;
  $output_dir =~ s/(\/)+$//; # Remove trailing forward-slashes if any

  if(!(-e $output_dir)){
      mkdir($output_dir);
  }

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
  print STDERR "processing model with sample_name: " . $sample_name . "\n";
  my $tumor_bam = $tumor_model->last_succeeded_build->whole_rmdup_bam_file;
  my $normal_bam = $normal_model->last_succeeded_build->whole_rmdup_bam_file;
  my $build_dir = $build->data_directory;



  # Check if the necessary files exist in this build

  my $snv_file = "$build_dir/effects/snvs.hq.novel.tier1.v2.bed";
  unless( -e $snv_file ){
      die "ERROR: SNV results for $sample_name not found at $snv_file\n";
  }
  my $indel_file = "$build_dir/effects/indels.hq.novel.tier1.v2.bed";
  unless( -e $indel_file ){
      die "ERROR: INDEL results for $sample_name not found at $indel_file\n";
  }  
  my $sv_file;
  if($self->process_svs){
      $sv_file = "$build_dir/variants/svs.hq";
      unless( -e $sv_file ){
          die "ERROR: SV results for $sample_name not found at $sv_file\n";
      }
  }


  # create subdirectories, get files in place
  
  #if multiple models with the same name, add a suffix
  if( -e "$output_dir/$sample_name" ){
      my $suffix = 1;
      my $newname = $sample_name . "-" . $suffix;
      while ( -e "$output_dir/$newname" ){
          $suffix++;
          $newname = $sample_name . "-" . $suffix;
      }
      $sample_name = $newname;
  }
  #make the directory structure
  mkdir "$output_dir/$sample_name";
  mkdir "$output_dir/$sample_name/snvs" unless( -e "$output_dir/$sample_name/snvs" );
  mkdir "$output_dir/$sample_name/indels" unless( -e "$output_dir/$sample_name/indels" );  
  mkdir "$output_dir/review" unless( -e "$output_dir/review" );
  `ln -s $build_dir $output_dir/$sample_name/build_directory`;

  #cat all the filtered snvs together (same for indels)
  `cat $build_dir/effects/snvs.hq.novel.tier*.v2.bed $build_dir/effects/snvs.hq.previously_detected.tier*.v2.bed | joinx sort >$output_dir/$sample_name/snvs/snvs.hq.bed`;
  `cat $build_dir/effects/indels.hq.novel.tier*.v2.bed $build_dir/effects/indels.hq.previously_detected.tier*.v2.bed | joinx sort >$output_dir/$sample_name/indels/indels.hq.bed`;

#  `ln -s $snv_file $output_dir/$sample_name/snvs/` unless( -e "$output_dir/$sample_name/snvs/$snv_file");
#  `ln -s $indel_file $output_dir/$sample_name/indels/` unless( -e "$output_dir/$sample_name/indels/$indel_file");
  if($self->process_svs){
      `mkdir $output_dir/$sample_name/svs`;
      `ln -s $sv_file $output_dir/$sample_name/svs/svs.hq` unless( -e "$output_dir/$sample_name/svs/$sv_file");
  }
  $snv_file = "$output_dir/$sample_name/snvs/snvs.hq.bed";
  $indel_file = "$output_dir/$sample_name/indels/indels.hq.bed";



  my %dups;
  #munge through SNV file to remove duplicates and fix IUB codes
  #--------------------------------------------------------------
  my $filenm = $snv_file;
  $filenm =~ s/.bed/.clean.bed/g;

  open(SNVFILE,">$filenm") || die ("couldn't open filter file");
  my $inFh = IO::File->new( $snv_file ) || die "can't open file3\n";
  while( my $line = $inFh->getline )
  {
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

  open(INDELFILE,">$filenm") || die ("couldn't open filter file");
  $inFh = IO::File->new( $indel_file ) || die "can't open file4\n";
  while( my $line = $inFh->getline )
  {
      chomp($line);
      my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
      if($ref =~ /\//){
          ( $ref, $var ) = split(/\//, $ref);
      }
      $ref =~ s/0/-/g;
      $var =~ s/0/-/g;
      $ref =~ s/\*/-/g;
      $var =~ s/\*/-/g;
      my $key = join("\t",($chr, $start, $stop, $ref, $var ));
      unless (exists($dups{$key})){
          print INDELFILE $key . "\n";
      }
      $dups{$key} = 1;
  }
  close(INDELFILE);
  $indel_file = $filenm;


  #-------------------------------------------------
  #filter out the off-target regions, if target regions are available
  if($self->restrict_to_target_regions){
      print STDERR "Filtering out off-target regions...\n";
      my %targetRegions;
      
      my $featurelist_name;
      if(defined($model->tumor_model->target_region_set_name)){
          $featurelist_name = $model->tumor_model->target_region_set_name;
          my $featurelist = Genome::FeatureList->get(name=>$featurelist_name)->file_path;
          if ( -s $featurelist ){
              #clean up feature list
              open(FEATFILE,">$output_dir/$sample_name/featurelist.tmp");
              my $inFh = IO::File->new( $featurelist ) || die "can't open file feature file\n";
              while( my $line = $inFh->getline )
              {
                  chomp($line);
                  next if $line =~ /^track/;
                  my ( $chr, $start, $stop, @rest) = split( /\t/, $line );
                  #remove chr if present
                  $chr =~ s/^chr//g;
                  print FEATFILE join("\t",( $chr, $start, $stop, @rest)) . "\n";
              }
              close($inFh);
              close(FEATFILE);
              `joinx sort $output_dir/$sample_name/featurelist.tmp >$output_dir/$sample_name/featurelist`;
              `rm -f $output_dir/$sample_name/featurelist.tmp`;
              `joinx intersect -a $snv_file -b $output_dir/$sample_name/featurelist >$snv_file.ontarget`;
              $snv_file = "$snv_file.ontarget";
              `joinx intersect -a $indel_file -b $output_dir/$sample_name/featurelist >$indel_file.ontarget`;
              $indel_file = "$indel_file.ontarget";

              
              
              # #compare the snvs to the targets
              # open(TARFILE,">$snv_file.ontarget") || die ("couldn't open target file");
              # $inFh = IO::File->new( $snv_file ) || die "can't open file\n";
              # while( my $line = $inFh->getline )
              # {
              #     chomp($line);
              #     my ( $chr, $start, $stop, @rest ) = split( /\t/, $line );
                  
              #     #if we run into huge lists, this will be slow - refactor to use joinx - TODO
              #     my $found = 0;
              #     foreach my $pos (keys(%{$targetRegions{$chr}})){
              #         my ($tst, $tsp) = split("\t",$pos);
              #         if(intersects($start,$stop,$tst,$tsp)){
              #             $found = 1;
              #         }
              #     }
              #     if($found){
              #         print TARFILE $line . "\n";
              #     }
              # }
              # close($inFh);
              # close(TARFILE);
              # $snv_file = "$snv_file.ontarget";
              
              
              # #compare the indels to the targets
              # open(TARFILE,">$indel_file.ontarget") || die ("couldn't open target file");
              # $inFh = IO::File->new( $indel_file ) || die "can't open file\n";
              # while( my $line = $inFh->getline )
              # {
              #     chomp($line);
              #     my ( $chr, $start, $stop, @rest ) = split( /\t/, $line );
              #     foreach my $pos (keys(%{$targetRegions{$chr}})){
              #         my ($tst, $tsp) = split("\t",$pos);
              #         if(intersects($start,$stop,$tst,$tsp)){
              #             print TARFILE $line . "\n";
              #         }
              #     }
              # }
              # close($inFh);
              # close(TARFILE);
              # $indel_file = "$indel_file.ontarget";
          } else {
              print STDERR "WARNING: feature list not found, No target region filtering being done\n";
          }      
      } else {
          print STDERR "No target region filtering being done (expected if this is WGS)\n";
      }
  }

  ##------------------------------------------------------
  #remove all but tier 1 sites, if that option is specified
  if($self->tier1_only){
      print STDERR "Doing Tiering...\n";
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
  
  ##------------------------------------------------------
  #remove all but tier 1 2 and 3 sites, if that option is specified
  if($self->tier123_only){
      print STDERR "Doing Tiering...\n";
      my $tier_cmd = Genome::Model::Tools::FastTier::FastTier->create(
          tier_file_location => $self->tier_file_location,
          variant_bed_file => $snv_file,
          );
      unless ($tier_cmd->execute) {
          die "Failed to tier variants successfully.\n";
      }
      `cat $snv_file.tier1 $snv_file.tier2 $snv_file.tier3 | joinx sort >$snv_file.tier1-3`;
      $snv_file = "$snv_file.tier1-3";
      `rm -f $snv_file.tier1 $snv_file.tier2 $snv_file.tier3 $snv_file.tier4`;

      $tier_cmd = Genome::Model::Tools::FastTier::FastTier->create(
          tier_file_location => $self->tier_file_location,
          variant_bed_file => $indel_file,
      );
      unless ($tier_cmd->execute) {
          die "Failed to tier variants successfully.\n";
      }
      `cat $indel_file.tier1 $indel_file.tier2 $indel_file.tier3 | joinx sort >$indel_file.tier1-3`;
      $indel_file = "$indel_file.tier1-3";
      `rm -f $indel_file.tier1 $indel_file.tier2 $indel_file.tier3 $indel_file.tier4`;
  }

  ##-------------------------------------------------
  #use joinx to remove dbsnp sites
  if(defined($self->dbsnp_filter)){
      print STDERR "Applying dbsnp filter...\n";

      #snvs - convert back to slashed bed file before joinx intersect
      my ($fh,$temp_bed_file) = Genome::Sys->create_temp_file;
      annoFileToSlashedBedFile($fh, $snv_file);
      close($fh);

      my ($fh2,$temp_novel_file) = Genome::Sys->create_temp_file;
      my $cmd = "joinx intersect --dbsnp-match --miss-a $temp_novel_file -a $temp_bed_file -b " . $self->dbsnp_filter . ">/dev/null";
      my $result = Genome::Sys->shellcmd(
        cmd => "$cmd",
          );
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
      $result = Genome::Sys->shellcmd(
        cmd => "$cmd",
          );
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
  if(defined($self->filter_sites)){
      print STDERR "Applying user-supplied filter...\n";
      my @filters = split(",",$self->filter_sites);
      my %filterSites;

      foreach my $filterfile (@filters){
          if( -e $filterfile){
              #store sites to filter out in a hash
              my $inFh = IO::File->new( $filterfile ) || die "can't open file5\n";
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

                  my @vars = fixIUB($ref, $var);
                  foreach my $v (@vars){
                      $filterSites{join("\t",($chr, $start, $stop, $ref, $v ))} = 0;
                  }
              }
              close($inFh);
          } else {
              die("filter sites file does not exist: " . $filterfile);
          }

          #remove snvs
          open(FILFILE,">$snv_file.filtered") || die ("couldn't open filter file");
          $inFh = IO::File->new( $snv_file ) || die "can't open file6\n";
          while( my $line = $inFh->getline )
          {
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
          open(FILFILE,">$indel_file.filtered") || die ("couldn't open filter file");
          $inFh = IO::File->new( $indel_file ) || die "can't open file7\n";
          while( my $line = $inFh->getline )
          {
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
          my $inFh2 = IO::File->new( $filterfile ) || die "can't open file8\n";
          while( my $line = $inFh2->getline )
          {
              print $fh $line;
          }
          close $inFh2;
      }
      close($fh);

      print STDERR "Removing user-specified filter...\n";
      my $cmd = "joinx intersect --miss-a $snv_file.filteredReg -a $snv_file -b $temp_file >/dev/null";
      my $result = Genome::Sys->shellcmd(
          cmd => "$cmd",
          );
      unless($result) {
          $self->error_message("Failed to execute joinx: Returned $result");
          die $self->error_message;
      }
      $snv_file = "$snv_file.filteredReg";

      $cmd = "joinx intersect --miss-a $indel_file.filteredReg -a $indel_file -b $temp_file >/dev/null";
      $result = Genome::Sys->shellcmd(
          cmd => "$cmd",
          );
      unless($result) {
          $self->error_message("Failed to execute joinx: Returned $result");
          die $self->error_message;
      }
      $indel_file = "$indel_file.filteredReg";
  }

  #-------------------------------------------------
  #auto-pass sites specified by the user
  if(defined($self->sites_to_pass)){
      print STDERR "Automatically-passing sites...\n";

      my %passSites;
      if( -e  $self->sites_to_pass){
          #store sites to pass in a hash
          my $inFh = IO::File->new(  $self->sites_to_pass) || die "can't open file\n";
          while( my $line = $inFh->getline )
          {
              chomp($line);
              my $bedline = annoToBed($line);
              my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $bedline );
              $ref =~ s/0/-/g;
              $var =~ s/0/-/g;
              
              $passSites{join("\t",($chr, $start, $stop, $ref, $var ))} = 0;
          }
          close($inFh);
      } else {
          die("sites-to-pass file does not exist: " . $self->sites_to_pass);
      }

      
      #remove snvs
      open(FILFILE,">$snv_file.autopassed") || die ("couldn't open autopass snvs outfile");
      open(FILFILE2,">$snv_file.notautopassed") || die ("couldn't open not autopass snvs outfile");
      $inFh = IO::File->new( $snv_file ) || die "can't open file6\n";
      while( my $line = $inFh->getline )
      {
          chomp($line);
          my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
          if($ref =~ /\//){
              ( $ref, $var ) = split(/\//, $ref);
          }
          if (exists($passSites{join("\t",($chr, $start, $stop, $ref, $var ))})){
              print FILFILE $line . "\n";
          } else {
              print FILFILE2 $line . "\n";
          }
      }
      close(FILFILE);
      $snv_file = "$snv_file.notautopassed";
          
      #remove indels
      open(FILFILE,">$indel_file.autopassed") || die ("couldn't open autopass indels outfile");
      open(FILFILE2,">$indel_file.notautopassed") || die ("couldn't open autopass indels outfile");
      $inFh = IO::File->new( $indel_file ) || die "can't open file7\n";
      while( my $line = $inFh->getline )
      {
          chomp($line);
          my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
          if($ref =~ /\//){
              ( $ref, $var ) = split(/\//, $ref);
          }
          
          if(exists($passSites{join("\t",($chr, $start, $stop, $ref, $var ))})){
              print FILFILE $line . "\n";
          } else {
              print FILFILE2 $line . "\n";
          }
      }
      close(FILFILE);
      $indel_file = "$indel_file.notautopassed";
  }

  ##
  ## This isn't needed - pindel vaf filter solved most of this problem
  ##
  # #-------------------------------------------------
  # #remove pindel-called indels, if they exist and aren't called by the other indel callers
  # my $pindel_calls = glob( "$build_dir/variants/indel/pindel-*/pindel-somatic-calls-*/pindel-vaf-filter-*/pindel-read-support-*/indels.hq.bed" );
  
  # if($pindel_calls && -s $pindel_calls){
  #     my $tmpdir = "$output_dir/$sample_name/indels/tmp";
  #     mkdir $tmpdir unless ( -e $tmpdir );
      
  #     my @indel_callers = glob("$build_dir/variants/indel/*");
      
  #     open(PFILE,">>$tmpdir/pindelcalls");
  #     open(OFILE,">>$tmpdir/othercalls");
      
  #     foreach my $dir (@indel_callers){
  #         next if $dir =~ /union|intersect/;
  #         next if $dir =~ /pindel/;
          
  #         #grab the non-pindel variants
  #         my $calls;
  #         if($dir =~ /gatk/){
  #             $calls = "$dir/indels.hq.bed";
  #         } elsif ($dir =~ /varscan/){
  #             $calls = glob("$dir/varscan-high-confidence-indel*/false-indel*/indels.hq.bed");
  #         } else {
  #             $self->error_message("can't parse indel caller directory $dir\nEdit the code to provide the path to final filtered indels.hq.bed for this caller");
  #             return 0;
  #         }

  #         my $infile = IO::File->new( $calls ) || die "can't open file $calls\n";
  #         while( my $line = $infile->getline )
  #         {
  #             chomp($line);
  #             my ( $chr, $start, $stop, $refvar ) = split( /\t/, $line );
  #             $refvar =~ s/\*/-/g;
  #             $refvar =~ s/0/-/g;
  #             print OFILE join("\t",( $chr, $start, $stop, $refvar )) . "\n";
  #         } 
  #         close($infile);
  #     }

  #     my $infile = IO::File->new( $pindel_calls ) || die "can't open file10\n";
  #     while( my $line = $infile->getline )
  #     {
  #         chomp($line);
  #         my ( $chr, $start, $stop, $refvar ) = split( /\t/, $line );
  #         $refvar =~ s/\*/-/g;
  #         $refvar =~ s/0/-/g;
  #         print PFILE join("\t",( $chr, $start, $stop, $refvar )) . "\n";
  #     } 
  #     close($infile);
      
  #     ##TODO -wrap these commands properly,
  #     `joinx sort $tmpdir/othercalls | uniq >$tmpdir/othercalls.sorted`;
  #     `joinx sort $tmpdir/pindelcalls | uniq >$tmpdir/pindelcalls.sorted`;

  #     print STDERR "Splitting out pindel indels, since they can't be reviewed...\n";
  #     my $cmd = "joinx intersect --miss-a $indel_file.pindel -a $indel_file -b $tmpdir/othercalls.sorted >$indel_file.non_pindel";
  #     my $result = Genome::Sys->shellcmd(
  #         cmd => "$cmd",
  #         );
  #     unless($result) {
  #         $self->error_message("Failed to execute joinx: Returned $result");
  #         die $self->error_message;
  #     }
  #     $indel_file = "$indel_file.non_pindel";
  # }



  #-------------------------------------------------------
  #get readcounts
  if($self->get_readcounts){
      print STDERR "Getting readcounts...\n";
      mkdir "$output_dir/$sample_name/snvs/readcounts";
      my $dir = "$output_dir/$sample_name/snvs/";
      if( -s "$dir/$snv_file" ){
          #get readcounts from the normal bam
          my $normal_rc_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
              bam_file => $normal_bam,
              output_file => "$dir/readcounts/$snv_file.nrm.cnt",
              variant_file => "$dir/$snv_file",
              genome_build => $ref_seq_fasta,
              );
          unless ($normal_rc_cmd->execute) {
              die "Failed to obtain normal readcounts for file $snv_file.\n";
          }      
          
          #get readcounts from the tumor bam
          my $tumor_rc_cmd = Genome::Model::Tools::Analysis::Coverage::BamReadcount->create(
              bam_file => $tumor_bam,
              output_file => "$dir/readcounts/$snv_file.tum.cnt",
              variant_file => "$dir/$snv_file",
              genome_build => $ref_seq_fasta,
              );
          unless ($tumor_rc_cmd->execute) {
              die "Failed to obtain tumor readcounts for file $snv_file.\n";
          }
      }
  }

  #-------------------------------------------------
  # run UHC on sites

  print "Gathering new sites...\n";

  if ( -s "$snv_file"){
      print "Running UHC filter...\n";

      #convert to annotation format:
      `perl /gscuser/cmiller/oneoffs/bedToAnnotation.pl $snv_file >$snv_file.var`;

      #run the uhc filter to remove solid calls
      my $uhc_cmd = Genome::Model::Tools::Somatic::UltraHighConfidence->create(
          normal_bam_file => $normal_bam,
          tumor_bam_file => $tumor_bam,
          output_file => "$snv_file.var.passuhc",
          variant_file => "$snv_file.var",
          reference => $ref_seq_fasta,
          filtered_file => "$snv_file.var.failuhc",
          );
      unless ($uhc_cmd->execute) {
          die "Failed to run UHC filter.\n";
      }

      #now get the files together for review
      print "Generating Review files...\n";
      my $revfile;
      if ( -s "$snv_file.var.failuhc"){
          $revfile = "$snv_file.var.failuhc";
      } else {
          $revfile = "$snv_file";
      }

      open(OUTFILE2,">$output_dir/review/$sample_name.tmp") || die "couldn't open outfile";
      #add the snvs
      $inFh = IO::File->new( $revfile ) || die "can't open file11\n";
      while( my $line = $inFh->getline )
      {
          chomp($line);
          my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
          print OUTFILE2 annoToBed(join("\t",($chr, $start, $stop, $ref . "/" . $var ))) . "\n";
      }

      #and the indels
      $inFh = IO::File->new( $indel_file ) || die "can't open file11\n";
      while( my $line = $inFh->getline )
      {
          my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
          print OUTFILE2 join("\t",($chr, $start, $stop, $ref . "/" . $var )) . "\n";
      }
      close(OUTFILE2);
      `joinx sort $output_dir/review/$sample_name.tmp > $output_dir/review/$sample_name.bed`;
      `rm -f $output_dir/review/$sample_name.tmp`;


      my $bam_files;
      my $labels;
      $bam_files = join(",",($normal_bam,$tumor_bam));
      $labels = join(",",("normal $sample_name","tumor $sample_name"));

      #create the xml file for review
      my $dumpXML = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
          bams => "$bam_files",
          labels => "$labels",
          output_file => "$output_dir/review/$sample_name.xml",
          genome_name => $sample_name,
          review_bed_file => "$output_dir/review/$sample_name.bed",
          reference_name => $self->igv_reference_name,
          );
      unless ($dumpXML->execute) {
          die "Failed to dump IGV xml for poorly covered sites.\n";
      }

      print STDERR "\n--------------------------------------------------------------------------------\n";
      print STDERR "Sites to review are here:\n";
      print STDERR "$output_dir/review/$sample_name.bed\n";
      print STDERR "IGV XML file is here:";
      print STDERR "$output_dir/review/$sample_name.xml\n\n";
      
      

  } else {
      print STDERR "No variants found\n";
  }

  return 1;
}

1;
