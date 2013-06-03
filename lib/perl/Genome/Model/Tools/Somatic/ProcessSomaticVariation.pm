package Genome::Model::Tools::Somatic::ProcessSomaticVariation;

use warnings;
use strict;
use IO::File;
use Genome;
use Sort::Naturally qw(nsort);
use Genome::Info::IUB;
use Spreadsheet::WriteExcel;

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
          example_values => ["reference_build36","b37","mm9"],
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
          doc => "add readcounts to the final variant list",
      },

      restrict_to_target_regions =>{
          is => 'Boolean',
          is_optional => 1,
          default => 1,
          doc => "only keep snv calls within the target regions. These are pulled from the build if possible",
      },

      tier_file_location =>{
          is => 'String',
          is_optional => 1,
          doc => "if provided, will add tiering information to the output.",
      },

      variant_list_is_one_based => {
          is => 'Boolean',
          is_optional => 1,
          doc => "The variant list you provide is in annotation (one-based) format, instead of bed. This flag fixes that.",
          default => 0,
      },

      add_dbsnp_and_gmaf => {
          is => 'Boolean',
          is_optional => 1,
          default => 1,
          doc => "if this is a recent build with vcf files (Jan 2013 or later), will add the rsids and GMAF information for all SNVs",
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

      create_review_files => {
          is => 'Boolean',
          is_optional => 1,
          doc => "create xml and bed files for manual review",
          default => 0,
      },

      create_archive => {
          is => 'Boolean',
          is_optional => 1,
          doc => "create an archive suitable for passing to collaborators",
          default => 0,
      },

      include_vcfs_in_archive => {
          is => 'Boolean',
          is_optional => 1,
          doc => "include full vcf files in archive (very large files)",
          default => 0,
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
  my $annotation_build_name = $model->annotation_build->name;
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
      my @sv_files = glob("$build_dir/variants/sv/union-union-sv_breakdancer_*sv_squaredancer*/svs.merge.file.somatic");
      $sv_file = $sv_files[0];
      unless( -e $sv_file ){
          print STDERR "ERROR: SV results for $sample_name not found, skipping SVs\n";
          $self->process_svs = 0;
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
  mkdir "$output_dir/review" unless( -e "$output_dir/review" || !($self->create_review_files));
  `ln -s $build_dir $output_dir/$sample_name/build_directory`;

  #cat all the filtered snvs together (same for indels)
  `cat $build_dir/effects/snvs.hq.novel.tier*.v2.bed $build_dir/effects/snvs.hq.previously_detected.tier*.v2.bed | joinx sort >$output_dir/$sample_name/snvs/snvs.hq.bed`;
  `cat $build_dir/effects/indels.hq.novel.tier*.v2.bed $build_dir/effects/indels.hq.previously_detected.tier*.v2.bed | joinx sort >$output_dir/$sample_name/indels/indels.hq.bed`;

#  `ln -s $snv_file $output_dir/$sample_name/snvs/` unless( -e "$output_dir/$sample_name/snvs/$snv_file");
#  `ln -s $indel_file $output_dir/$sample_name/indels/` unless( -e "$output_dir/$sample_name/indels/$indel_file");
  if($self->process_svs){
      `mkdir $output_dir/$sample_name/svs`;
      `ln -s $sv_file $output_dir/$sample_name/svs/svs.hq` unless( -e "$output_dir/$sample_name/svs/$sv_file");

#       #annotate the svs
#       my $anno_cmd = Genome::Model::Tools::Annotate::Sv::Combine->create(
#           input-file => $sv_file,
#           output-file => "$output_dir/$sample_name/svs/svs.hq.annotated",
#           annotation-build
#           dbsnp-annotation-file /gsc/scripts/opt/genome/db/genome-db-dbsnp/human/build/37/132/dbsnp.csv
#           dbvar-annotation-file /gsc/scripts/opt/genome/db/dbvar/human/build37/dbvar.tsv
#           fusion-transcripts-fusion-output-file /tmp/zout.fusions
#           repeat-masker-annotation-file /gscuser/aregier/git/genome/vep/lib/perl/Genome/repeat_masker.tsv
#           annotator-list=Transcripts,FusionTranscripts,Dbsnp,Segdup,Dbvar
#           segdup-annotation-file /gsc/scripts/opt/genome/db/ucsc/human/build37/segdup.tsv
#           chrA-column 1
#           bpA-column 2
#           chrB-column 4
#           bpB-column 5
#           event-type-column 7
#           score-column 12
#           orient-column 8

#           );
#       unless ($anno_cmd->execute) {
#           die "Failed to annotate sv file\n";
#       }
# #gmt annotate sv combine --input-file /tmp/zin --output-file /tmp/zout1 --annotation-build 124434505 --dbsnp-annotation-file /gsc/scripts/opt/genome/db/genome-db-dbsnp/human/build37/132/dbsnp.csv --dbvar-annotation-file /gsc/scripts/opt/genome/db/dbvar/human/build37/dbvar.tsv --fusion-transcripts-fusion-output-file /tmp/zout.fusions --repeat-masker-annotation-file /gscuser/aregier/git/genome/vep/lib/perl/Genome/repeat_masker.tsv --annotator-list=Transcripts,FusionTranscripts,Dbsnp,Segdup,Dbvar --segdup-annotation-file /gsc/scripts/opt/genome/db/ucsc/human/build37/segdup.tsv --chrA-column 1 --bpA-column 2 --chrB-column 4 --bpB-column 5 --event-type-column 7 --score-column 12 --orient-column 8


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

          } else {
              print STDERR "WARNING: feature list not found, No target region filtering being done\n";
          }
      } else {
          print STDERR "No target region filtering being done (expected if this is WGS)\n";
      }
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


  ##------------------------------------------------------
  #do tiering and annotation

  print STDERR "Tiering...\n";
  my $tier_cmd = Genome::Model::Tools::FastTier::FastTier->create(
      tier_file_location => $self->tier_file_location,
      variant_bed_file => $snv_file,
      );
  unless ($tier_cmd->execute) {
      die "Failed to tier variants successfully.\n";
  }
  $tier_cmd = Genome::Model::Tools::FastTier::FastTier->create(
      tier_file_location => $self->tier_file_location,
      variant_bed_file => $indel_file,
      );
  unless ($tier_cmd->execute) {
      die "Failed to tier variants successfully.\n";
  }

  my $snv_file2 = $snv_file;
  $snv_file2 =~ s/\.bed//g;
  my $indel_file2 = $indel_file;
  $indel_file2 =~ s/\.bed//g;

  #keep tier 1 through 3 sites
  `cat $snv_file.tier1 $snv_file.tier2 $snv_file.tier3 | joinx sort >$snv_file2.tier1-3.bed`;
  my $tsnvfile=$snv_file;
  $snv_file = "$snv_file2.tier1-3";

  `cat $indel_file.tier1 $indel_file.tier2 $indel_file.tier3 | joinx sort >$indel_file2.tier1-3.bed`;
  my $tindelfile=$indel_file;
  $indel_file = "$indel_file2.tier1-3";

  #convert to 1-based annotation format
  open(OUTFILE,">$snv_file2.tier1-3");
  $inFh = IO::File->new( "$snv_file2.tier1-3.bed" ) || die "can't open file2b\n";
  while( my $line = $inFh->getline )
  {
      chomp($line);
      my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
      print OUTFILE bedToAnno(join("\t",($chr, $start, $stop, $ref, $var))) . "\n";
  }
  close(OUTFILE);

  open(OUTFILE,">$indel_file2.tier1-3");
  $inFh = IO::File->new( "$indel_file2.tier1-3.bed" ) || die "can't open file2c\n";
  while( my $line = $inFh->getline )
  {
      chomp($line);
      my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
      print OUTFILE bedToAnno(join("\t",($chr, $start, $stop, $ref, $var))) . "\n";
  }
  close(OUTFILE);
  $snv_file = "$snv_file2.tier1-3";
  $indel_file = "$indel_file2.tier1-3";


  #-------------------------------
  # annotation
  print STDERR "Doing Annotation...\n";
  my $anno_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
      variant_file => $snv_file,
      output_file => $snv_file . ".anno",
      reference_transcripts => $annotation_build_name,
      annotation_filter => "top",
          );
  unless ($anno_cmd->execute) {
      die "Failed to annotate variants successfully.\n";
  }
  $anno_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
      variant_file => $indel_file,
      output_file => $indel_file . ".anno",
      reference_transcripts => $annotation_build_name,
      annotation_filter => "top",
      );
  unless ($anno_cmd->execute) {
      die "Failed to annotate variants successfully.\n";
  }

  $snv_file = "$snv_file.anno";
  $indel_file = "$indel_file.anno";

  #----------------------------------------------------
  # add dbsnp/gmaf

  if ($self->add_dbsnp_and_gmaf){
      print STDERR "==== adding dbsnp ids ====\n";
      print STDERR "$build_dir/variants/snvs.annotated.vcf.gz\n";
      if(-s "$build_dir/variants/snvs.annotated.vcf.gz"){
          my $db_cmd = Genome::Model::Tools::Annotate::AddRsid->create(
              anno_file => $snv_file,
              output_file => "$snv_file.rsid",
              vcf_file => "$build_dir/variants/snvs.annotated.vcf.gz",
              );
          unless ($db_cmd->execute) {
              die "Failed to add dbsnp anno to file $snv_file.\n";
          }
          $snv_file = "$snv_file.rsid";
          #pad indel file with tabs to match - if we ever start annotating with indels from dbsnp, replace this section
          open(OUTFILE,">$indel_file.rsid");
          my $inFh = IO::File->new( "$indel_file" ) || die "can't open file2d\n";
          while( my $line = $inFh->getline ){
              chomp($line);
              print OUTFILE $line . "\t\t\n"
          }
          close($inFh);
          close(OUTFILE);

          $snv_file = "$snv_file.rsid";
          $indel_file = "$indel_file.rsid";
      } else {
          print STDERR "Warning: couldn't find annotated SNV file in build, skipping dbsnp anno\n";
      }
  }


  #-----------------------------------------------------
  #add tier info as a column
  my %tiers;
  foreach my $tier ("tier1","tier2","tier3"){
      my $inFh = IO::File->new( "$tsnvfile.$tier" ) || die "can't open file2d\n";
      while( my $line = $inFh->getline ){
          chomp($line);
          my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
          $tiers{bedToAnno(join("\t",($chr, $start, $stop, $ref, $var)))} = $tier
      }
      close($inFh);
      $inFh = IO::File->new( "$tindelfile.$tier" ) || die "can't open file2e\n";
      while( my $line = $inFh->getline ){
          chomp($line);
          my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
          $tiers{bedToAnno(join("\t",($chr, $start, $stop, $ref, $var)))} = $tier
      }
      close($inFh);
  }

  open(OUTFILE,">$snv_file.tiered");
  $inFh = IO::File->new( "$snv_file" ) || die "can't open file2f\n";
  while( my $line = $inFh->getline ){
      chomp($line);
      if($line =~ /^chromo/){
          print OUTFILE $line . "\ttier\n";
          next;
      }
      my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
      my $key = join("\t",($chr, $start, $stop, $ref, $var));
      if(defined($tiers{$key})){
          print OUTFILE $line . "\t" . $tiers{$key} . "\n";
      } else {
          print STDERR "WARNING: site $key not tiered\n";
          print OUTFILE $line . "\t\n";
      }
  }
  $snv_file = "$snv_file.tiered";
  close(OUTFILE);
  open(OUTFILE,">$indel_file.tiered");
  $inFh = IO::File->new( "$indel_file" ) || die "can't open file2g\n";
  while( my $line = $inFh->getline ){
      chomp($line);
      if($line =~ /^chromo/){
          print OUTFILE $line . "\ttier\n";
          next;
      }
      my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
      my $key = join("\t",($chr, $start, $stop, $ref, $var));
      if(defined($tiers{$key})){
          print OUTFILE $line . "\t" . $tiers{$key} . "\n";
      } else {
          print STDERR "WARNING: site $key not tiered\n";
          print OUTFILE $line . "\t\n";
      }
  }
  $indel_file = "$indel_file.tiered";
  close(OUTFILE);


  #-------------------------------------------------------
  #get readcounts
  if($self->get_readcounts){
      print STDERR "Getting readcounts...\n";
      if( -s "$snv_file" ){
          #get readcounts from the normal and tumor bams
          my $rc_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
              bam_files => "$normal_bam,$tumor_bam",
              output_file => "$snv_file.rcnt",
              variant_file => "$snv_file",
              genome_build => $ref_seq_fasta,
              header_prefixes => "Normal,Tumor",
              indel_size_limit => 5,
              );
          unless ($rc_cmd->execute) {
              die "Failed to obtain readcounts for file $snv_file.\n";
          }
      }
      if( -s "$indel_file" ){
          #get readcounts from the normal and tumor bams
          my $rc_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(
              bam_files => "$normal_bam,$tumor_bam",
              output_file => "$indel_file.rcnt",
              variant_file => "$indel_file",
              genome_build => $ref_seq_fasta,
              header_prefixes => "Normal,Tumor",
              indel_size_limit => 5,
              );
          unless ($rc_cmd->execute) {
              die "Failed to obtain readcounts for file $indel_file.\n";
          }
      }
      $snv_file = "$snv_file.rcnt";
      $indel_file = "$indel_file.rcnt";
  }



  #------------------------------------------------------
  # combine the files into one master table
  `head -n 1 $snv_file >$output_dir/$sample_name/snvs.indels.annotated`;
  `tail -n +2 $indel_file >>$output_dir/$sample_name/snvs.indels.annotated.tmp`;
  `tail -n +2 $snv_file >>$output_dir/$sample_name/snvs.indels.annotated.tmp`;
  `joinx sort -i $output_dir/$sample_name/snvs.indels.annotated.tmp >>$output_dir/$sample_name/snvs.indels.annotated`;
  `rm -f $output_dir/$sample_name/snvs.indels.annotated.tmp`;

  # convert master table to excel
  my $workbook  = Spreadsheet::WriteExcel->new("$output_dir/$sample_name/snvs.indels.annotated.xls");
  my $worksheet = $workbook->add_worksheet();

  my $row=0;
  $inFh = IO::File->new( "$output_dir/$sample_name/snvs.indels.annotated" ) || die "can't open file\n";
  while( my $line = $inFh->getline )
  {
      chomp($line);
      my @F = split("\t",$line);
      for(my $i=0;$i<@F;$i++){
          $worksheet->write($row, $i, $F[$i]);
      }
      $row++;
  }
  close($inFh);
  $workbook->close();


  #------------------------------------------------------
  #now get the files together for review
  if($self->create_review_files){
      print "Generating Review files...\n";
      my $revfile = "$snv_file";

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
      $inFh = IO::File->new( $indel_file ) || die "can't open file12\n";
      while( my $line = $inFh->getline )
      {
          my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
          print OUTFILE2 annoToBed(join("\t",($chr, $start, $stop, $ref . "/" . $var ))) . "\n";
      }
      close(OUTFILE2);
      `joinx sort $output_dir/review/$sample_name.tmp > $output_dir/review/$sample_name.bed`;
      `rm -f $output_dir/review/$sample_name.tmp`;


      my $bam_files;
      my $labels;
      $bam_files = join(",",($normal_bam,$tumor_bam));
      $labels = join(",",("normal $sample_name","tumor $sample_name"));

      my $igv_reference_name = "b37";
      if(defined($self->igv_reference_name)){
          $igv_reference_name = $self->igv_reference_name;
      } else {
          print STDERR "WARNING: No IGV reference name supplied - defaulting to build 37\n";
      }

      #create the xml file for review
      my $dumpXML = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
          bams => "$bam_files",
          labels => "$labels",
          output_file => "$output_dir/review/$sample_name.xml",
          genome_name => $sample_name,
          review_bed_file => "$output_dir/review/$sample_name.bed",
          reference_name => $igv_reference_name,
          );
      unless ($dumpXML->execute) {
          die "Failed to dump IGV xml for poorly covered sites.\n";
      }

      print STDERR "\n--------------------------------------------------------------------------------\n";
      print STDERR "Sites to review are here:\n";
      print STDERR "$output_dir/review/$sample_name.bed\n";
      print STDERR "IGV XML file is here:";
      print STDERR "$output_dir/review/$sample_name.xml\n\n";
  }

  #------------------------------------------------
  # tar up the files to be sent to collaborators
  #
  if($self->create_archive){
      mkdir("$output_dir/$sample_name/$sample_name");

      chdir("$output_dir/$sample_name/");
      #VCF files
      if($self->include_vcfs_in_archive){
          if(-e "$build_dir/variants/indels.detailed.vcf.gz"){
              `ln -s $build_dir/variants/indels.detailed.vcf.gz $sample_name/indels.vcf.gz`;
          } elsif(-e "$build_dir/variants/indels.vcf.gz") {
              `ln -s $build_dir/variants/indels.vcf.gz $sample_name/indels.vcf.gz`;
          } else {
              print STDERR "WARNING: no indel VCF file available. If this is an older model, a rebuild may fix this\n";
          }
          if(-e "$build_dir/variants/snvs.annotated.vcf.gz"){
              `ln -s $build_dir/variants/snvs.annotated.vcf.gz  $sample_name/snvs.vcf.gz`;
          } elsif (-e "$build_dir/variants/snvs.vcf.gz"){
              `ln -s $build_dir/variants/snvs.vcf.gz  $sample_name/snvs.vcf.gz`;
          } else {
              print STDERR "WARNING: no snv VCF file available. If this is an older model, a rebuild may fix this\n";
          }
      }
      #annotated snvs and indels
      `ln -s ../snvs.indels.annotated $sample_name/snvsAndIndels.annotated`;
      #same in excel format
      `ln -s ../snvs.indels.annotated.xls $sample_name/snvsAndIndels.annotated.xls`;
      #sv calls
      if($self->process_svs){
          `ln -s $sv_file $output_dir/$sample_name/$sample_name/svs`;
          #`ln -s $sv_file $output_dir/$sample_name/$sample_name/svs.annotated`;
      }
      #cnv calls - todo

      #tar it up
      `tar -czvfh $sample_name.tar.gz $sample_name`;
  }

  return 1;
}

1;
