package Genome::Model::Tools::Validation::ProcessSomaticValidation;

use warnings;
use strict;
use IO::File;
use Genome;
use Sort::Naturally qw(nsort);
use Genome::Info::IUB;
use Spreadsheet::WriteExcel;
use File::Basename;
use File::Path qw(make_path remove_tree);
use Genome::Utility::Text qw(sanitize_string_for_filesystem);

class Genome::Model::Tools::Validation::ProcessSomaticValidation {
  is => 'Command',
  has_input => [
      output_dir => {
          is => 'Text',
          doc => "Directory where output will be stored (under a subdirectory with the sample name)",
      },
      ],

  has_optional_input => [
      somatic_validation_model_id => {
          is => 'Text',
          doc => "ID of SomaticValidation model - either this or the --somatic-validation-build-id are required",
          is_optional => 1,
      },

      somatic_validation_build_id => {
          is => 'Text',
          doc => "ID of SomaticValidation build - either this or the --somatic-validation-model-id are required",
          is_optional => 1,
      },

      somatic_variation_model_id => {
          is => 'Text',
          doc => "Optional. If provided, will create xml sessions for review that display both the original bams from this model along with the capture validation bams",
          is_optional => 1,
      },

      tumor_only => {
          is => 'Boolean',
          is_optional => 1,
          default => 0,
          doc => "Model only has tumor data",
      },

      use_assembled_indels => {
          is => 'Boolean',
          is_optional => 1,
          default => 1,
          doc => 'Use the assembled indel calls from the build. Turning this off uses indels straight out of discovery',
      },

      restrict_to_target_regions =>{
          is => 'Boolean',
          is_optional => 1,
          default => 1,
          doc => "only keep new calls within the target regions. If --target-regions isn't defined, these are pulled from the build if possible",
      },

      target_regions =>{
          is => 'String',
          is_optional => 1,
          doc => "path to a target file region. Used in conjunction with --restrict-to-target-regions to limit sites to those appearing in these regions",
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

      add_tiers =>{
          is => 'Boolean',
          is_optional => 1,
          default => 1,
          doc => "add tiering information to the output",
      },

      add_dbsnp_and_gmaf => {
          is => 'Boolean',
          is_optional => 1,
          default => 0,
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

      get_readcounts =>{
          is => 'Boolean',
          is_optional => 1,
          default => 1,
          doc => "add readcounts to the final variant list",
      },

      variant_list_is_one_based => {
          is => 'Boolean',
          is_optional => 1,
          doc => "The variant list you provide is in annotation (one-based) format, instead of bed. This flag fixes that.",
          default => 0,
      },

      create_review_files => {
          is => 'Boolean',
          is_optional => 1,
          doc => "create xml and bed files for manual review",
          default => 0,
      },

      igv_reference_name =>{
          is => 'Text',
          is_optional => 1,
          doc => "name of the igv reference to use",
          example_values => ["hg38","mm10","b37","mm9","reference_build36"],
      },

      tiers_to_review => {
          is => 'String',
          is_optional => 1,
          doc => "comma-separated list of tiers to include in review. (e.g. 1,2,3 will create bed files with tier1, tier2, and tier3).  Anything other than 1,2,3,4 (all tiers) requires that the --add-tiers option be specified",
          default => "1,2,3,4",
      },

      max_indel_size => {
          is => 'Integer',
          is_optional => 1,
          doc => "set max indel size to extract readcounts",
          default => "4",
      },

      statuses_to_review => {
          is => 'String',
          is_optional => 1,
          doc => "comma-separated list of statuses to include in review. Accepted values: [validated,newcall,notvalidated]",
          default => "validated,newcall",
      },

      create_archive => {
          is => 'Boolean',
          is_optional => 1,
          doc => "create an archive suitable for passing to collaborators",
          default => 0,
      },

      sample_name =>{
          is => 'Text',
          is_optional => 1,
          doc => "override the sample name on the build and use this name instead",
      },

      reference_transcripts => {
          is => 'Text',
          is_optional => 1,
          doc => "use this reference transcript build instead of the one specified in the model (e.g. NCBI-mouse.ensembl/67_37)",
      },

      bam_readcount_version => {
          is => 'Text',
          doc => "the version of bam-readcount to use if generating counts",
          is_optional => 1,
      },

      # include_vcfs_in_archive => {
      #     is => 'Boolean',
      #     is_optional => 1,
      #     doc => "include full vcf files in archive (very large files)",
      #     default => 0,
      # },

      ],
};


sub help_detail {
  return <<HELP;
Given a SomaticValidation model, this tool will gather the resulting variants, remove
off-target sites, tier the variants, optionally filter them, etc. Calls are prepped for
manual review in the review/ directory.
HELP
}

sub _doc_authors {
    return <<AUTHS;
    Chris Miller
AUTHS
}

sub bedFileToAnnoFile{
    my ($file,$outfile) = @_;


    #remove bed from name
    my $newfile = $file;
    $newfile =~ s/\.bed//g;

    if($outfile){
        $newfile = $outfile;
    }

    open(OUTFILE,">$newfile");
    my $inFh = IO::File->new( $file ) || die "can't open file2\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        if($line =~ /^chrom/){
            print OUTFILE $line . "\n";
            next;
        }
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        if($ref =~ /\//){
            my @alleles = split("/",$ref);
            print OUTFILE bedToAnno(join("\t",($chr, $start, $stop, $alleles[0], $alleles[1]))) . "\t" . join("\t",($var,@rest)) . "\n";
        } else {
            print OUTFILE bedToAnno(join("\t",($chr, $start, $stop, $ref, $var))) . "\t" . join("\t",@rest) . "\n";
        }
    }
    close(OUTFILE);
    close($inFh);
    return($newfile);
}


sub annoFileToBedFile{
    my ($file) = @_;
    #add bed to name
    my $newfile = $file . ".bed";

    open(OUTFILE,">$newfile");
    my $inFh = IO::File->new( $file ) || die "can't open file2\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        if($line =~ /^chrom/){
            print OUTFILE $line . "\n";
            next;
        }
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        print OUTFILE annoToBed(join("\t",($chr, $start, $stop, $ref, $var))) . "\t" . join("\t",@rest) . "\n";
    }
    close(OUTFILE);
    close($inFh);
    return($newfile);
}

sub annoFileToSlashedBedFile{
    my ($file,$outfile) = @_;

    my $newfile = $file . ".bed";
    if($outfile){
        $newfile = $outfile;
    }

    open(OUTFILE,">$newfile");
    my $inFh = IO::File->new( $file ) || die "can't open file2\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        if($line =~ /^chrom/){
            next;
        }
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        my $bed = annoToBed(join("\t",($chr, $start, $stop, $ref, $var)));
        my @bedline = split(/\t/,$bed);
        print OUTFILE join("\t",(@bedline[0..2],"$bedline[3]/$bedline[4]")) . "\n";
    }
    close(OUTFILE);
    close($inFh);
    return($newfile);
}

sub bedToAnno{
    my ($chr,$start,$stop,$ref,$var) = split("\t",$_[0]);
    if ($ref =~ /^[-0*]/){ #indel INS
        $stop = $stop+1;
    } else { #indel DEL or SNV
        $start = $start+1;
    }
    return(join("\t",($chr,$start,$stop,$ref,$var)));
}


sub annoToBed{
    my ($chr,$start,$stop,$ref,$var) = split("\t",$_[0]);
    if ($ref =~ /^[-*0]/){ #indel INS
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

sub addName{
    my ($file,$name) = @_;
    if($file =~ /\.bed$/){
        $file =~ s/\.bed//g;
        return($file . "." . $name . ".bed");
    } else {
        return($file . "." . $name);
    }
}

#read in the file, output a cleaned-up version
sub cleanFile{
    my ($file) = @_;

    my $newfile = addName($file,"clean");

    my %dups;
    open(OUTFILE, ">$newfile");

    my $inFh = IO::File->new( $file ) || die "can't open file $file\n";
    while( my $line = $inFh->getline ){
        chomp($line);
        my ( $chr, $start, $stop, $ref, $var, @rest ) = split( /\t/, $line );
        if($ref =~ /\//){
            ( $ref, $var ) = split(/\//, $ref);
        }

        $ref =~ s/0/-/g;
        $var =~ s/0/-/g;
        $ref =~ s/\*/-/g;
        $var =~ s/\*/-/g;

        my @vars = ($var);
        unless($ref =~ /-/ || $var =~ /-/){ #fixiub doesn't handle indels
            @vars = fixIUB($ref, $var);
        }

        foreach my $v (@vars){
            unless (exists($dups{join("\t",($chr, $start, $stop, $ref, $v ))})){
                print OUTFILE join("\t",($chr, $start, $stop, $ref, $v )) . "\n";
            }
            $dups{join("\t",($chr, $start, $stop, $ref, $v ))} = 1;
        }
    }
    close(OUTFILE);
    close($inFh);
    `joinx sort -i $newfile >$newfile.tmp`;
    `mv -f $newfile.tmp $newfile`;
    return($newfile)
}


sub getFilterSites{
    my $filter_string = @_;
    my @filters = split(",",$filter_string);
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
                $ref =~ s/\*/-/g;
                $var =~ s/\*/-/g;

                my @vars = fixIUB($ref, $var);
                foreach my $v (@vars){
                    $filterSites{join("\t",($chr, $start, $stop, $ref, $v ))} = 0;
                }
            }
            close($inFh);
        } else {
            die("filter sites file does not exist: " . $filterfile);
        }
    }
    return(\%filterSites);
}


sub removeFilterSites{
    my ($file,$filterSites) = @_;

    my $newfile = addName($file,"filtered");
    #handle zero size files
    if( -z $file ){
        `touch $newfile`;
        return($newfile);
    }

    open(FILFILE,">$newfile") || die ("couldn't open filter output file");
    my $inFh = IO::File->new( $file ) || die "can't open filter input file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $line );
        if($ref =~ /\//){
            ( $ref, $var ) = split(/\//, $ref);
        }
        unless (defined($filterSites->join("\t",($chr, $start, $stop, $ref, $var )))){
            print FILFILE $line . "\n";
        }
    }
    close(FILFILE);
    return($newfile)
}


sub doAnnotation{
    my ($file,$annotation_build_name) = @_;

    if($file =~ /.bed$/){
        $file = bedFileToAnnoFile($file);
    }

    #handle zero size files
    if( -z $file ){
        `touch $file.anno`;
        return($file . ".anno");
    }
    my $anno_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
        variant_file => $file,
        output_file => $file . ".anno",
        reference_transcripts => $annotation_build_name,
        annotation_filter => "top",
        );
    unless ($anno_cmd->execute) {
        die "Failed to annotate variants successfully.\n";
    }
    return($file . ".anno");
}


sub addTiering{
    my ($file, $tier_file_location) = @_;

    unless($file =~ /\.bed/){
        $file = annoFileToBedFile($file);
    }

    my $newfile = addName($file, "tiered");

    #handle zero size files
    if( -z $file ){
        `touch $newfile`;
        return($newfile);
    }

    my $tier_cmd = Genome::Model::Tools::FastTier::AddTiers->create(
        input_file => $file,
        output_file => $newfile,
        tier_file_location => $tier_file_location,
        );
    unless ($tier_cmd->execute) {
        die "Failed to tier variants successfully.\n";
    }
    return($newfile);
}


sub getReadcounts{
    my ($file, $ref_seq_fasta, $bams, $labels, $bam_readcount_version,$indel_size) = @_;
    #todo - should check if input is bed and do coversion if necessary

    my $output_file = "$file.rcnt";
    if( -s $file ){
        my $bamfiles = join(",",@$bams);
	my $header = join(",",@$labels);

        #get readcounts from the tumor bam only
        my %params = (
            bam_files => $bamfiles,
            output_file => $output_file,
            variant_file => $file,
            genome_build => $ref_seq_fasta,
            header_prefixes => $header,
            indel_size_limit => $indel_size,
            );
        if($bam_readcount_version){
            $params{bam_readcount_version} = $bam_readcount_version;
        }
        my $rc_cmd = Genome::Model::Tools::Analysis::Coverage::AddReadcounts->create(%params);
        unless ($rc_cmd->execute) {
            die "Failed to obtain readcounts for file $file.\n";
        }
    } else {
        `touch $output_file`;
    }

    return($output_file);
}
###############################################################################################################
sub execute {

  $DB::single = 1;
  my $self = shift;
  my $model;
  my $build;
  my $bam_readcount_version = $self->bam_readcount_version;
  

  #validate that the inputs are valid

  #check statuses
  my @statuses = split(",",$self->statuses_to_review);
  foreach my $i (@statuses){
      unless($i =~ /^newcall$|^validated$|^notvalidated$/){
	  die $self->error_message("status $i is not a valid status to review. Status must be one of [validated,notvalidated,newcall]");
      }
  }

  #check that a build or model was specified
  unless(defined($self->somatic_validation_build_id) || defined($self->somatic_validation_model_id)){
      die $self->error_message("must specify either somatic-validation-build-id or somatic-validation-model-id");
  }
  
  # Check on the input model and params
  if(defined($self->somatic_validation_model_id)){
      $model = Genome::Model->get( $self->somatic_validation_model_id );
      unless( defined $model ){
	  $self->error_message("Could not find a model with ID: %s", $self->somatic_validation_model_id);
          return undef;
      }
      
      $build = $model->last_succeeded_build;
      unless( defined($build) ){
	  $self->warning_message("Model %s has no succeeded builds", $model->id);
          return undef;
      }

  } elsif(defined($self->somatic_validation_build_id)){
      $build = Genome::Model::Build->get( $self->somatic_validation_build_id );
      unless( defined $build ){
	  $self->error_message("Could not find a build with ID: %s", $self->somatic_validation_build_id );
          return undef;
      }
      $model = $build->model;      
  }

  #get the target-region file
  if($self->restrict_to_target_regions){
      unless($model->target_region_set || $self->target_regions){        
          die $self->error_message("No target regions provided and no target_region_set defined on model. Can't use the --restrict-to-target-regions option on this model.");
      }
  }


  #make the output dir if it doesn't exist
  my $output_dir = $self->output_dir;
  $output_dir =~ s/(\/)+$//; # Remove trailing forward-slashes if any  

  if(!-d $output_dir){
      make_path($output_dir);
  }



  #make sure specified tiers are valid and that if anything other than
  # all (tiers 1,2,3,4) is specified, that add-tiers is on
  my @tiers = split(",",$self->tiers_to_review);
  my $tstring = "";  
  foreach my $t (sort(@tiers)){
      unless ($t =~/^[1234]$/){
	  die $self->error_message("$t is not a valid tier to review. tiers-to-review should be a comma-separated list containing only [1,2,3,4]");
                    
      }
      $tstring .= $t
  }
  if(($tstring ne "1234") && !($self->add_tiers)){
      die $self->error_message("if a combination of tiers-to-review other than 1,2,3,4 is specified, --add-tiers must be specified. (otherwise I don't know how to grab the tiers you want for review!");
  }


  my $ref_seq_build_id = $model->reference_sequence_build->build_id;
  my $ref_seq_build = Genome::Model::Build->get($ref_seq_build_id);
  my $ref_seq_fasta = $ref_seq_build->full_consensus_path('fa');
  my $annotation_build_name = $model->annotation_build->name;
  if(defined $self->reference_transcripts){
      $self->status_message("Model's annotation build overridden. Using %s", $self->reference_transcripts);
      $annotation_build_name = $self->reference_transcripts;
  }

  my $tiering_files = $model->annotation_build->data_directory . "/annotation_data/tiering_bed_files_v3/";
  my $sample_name;
  if(!defined($self->sample_name)){
      $sample_name = $model->subject->name;
  } else {
      $sample_name = $self->sample_name;
  }
  $sample_name = Genome::Utility::Text::sanitize_string_for_filesystem($sample_name);

  
  $self->status_message("processing model with sample_name: $sample_name\n");
  #retrieve validation BAMs
  my $tumor_bam = $build->tumor_bam;
  my $normal_bam;
  unless($self->tumor_only){
      $normal_bam = $build->normal_bam;
  }
  #retrieve wgs/wxs BAMs from somatic_variation model
  my $som_var_bam_files=[];
  my $som_var_labels=[];
  if(defined $self->somatic_variation_model_id){
      ($som_var_bam_files,$som_var_labels) = $self->process_somatic_variation_models($self->somatic_variation_model_id,$sample_name);
  
  }


  my $build_dir = $build->data_directory;

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
  make_path("$output_dir/$sample_name");
  make_path("$output_dir/$sample_name/snvs") unless( -d "$output_dir/$sample_name/snvs" );
  make_path("$output_dir/$sample_name/indels") unless( -d "$output_dir/$sample_name/indels" );
  make_path("$output_dir/review") unless( -d "$output_dir/review" || !($self->create_review_files));
  `ln -s $build_dir $output_dir/$sample_name/build_directory`;



  #get the files to operate on
  my @orig_snv_files = ("$build_dir/validation/snvs/snvs.validated","$build_dir/validation/snvs/snvs.notvalidated","$build_dir/validation/snvs/snvs.newcalls");
  my @snv_files;
  for my $file (@orig_snv_files){
      unless( -e $file ){
          die $self->error_message("SNV results for $sample_name not found at $file");
      }
      `cp $file $output_dir/$sample_name/snvs/`;
      #only operate on non-zero files
      if( -s $file ){
          push(@snv_files,"$output_dir/$sample_name/snvs/". basename($file));
      }
  }
  my ($small_indel_file, $large_indel_file) = $self->indel_files( $build_dir );
  #cat all the indels together into one file
  if($large_indel_file =~ /variants/){ #no validated indels, use the bed files
      bedFileToAnnoFile($large_indel_file,"$output_dir/$sample_name/indels/large.indels");

      if($small_indel_file eq $large_indel_file){
          `joinx sort -i $output_dir/$sample_name/indels/large.indels >$output_dir/$sample_name/indels/indels.hq`;
          `rm -f $output_dir/$sample_name/indels/large.indels`;
      } else {
          `grep -w Somatic $small_indel_file | cut -f 1-5 >$output_dir/$sample_name/indels/small.indels`;
          bedFileToAnnoFile("$output_dir/$sample_name/indels/small.indels.bed");
          `cat $output_dir/$sample_name/indels/large.indels $output_dir/$sample_name/indels/small.indels | joinx sort >$output_dir/$sample_name/indels/indels.hq`;
      }
  } else { #use the validated indels
      `cut -f 1-5 $large_indel_file >$output_dir/$sample_name/indels/large.indels`;
      `grep -w Somatic $small_indel_file | cut -f 1-5 >$output_dir/$sample_name/indels/small.indels`;
      `cat $output_dir/$sample_name/indels/small.indels $output_dir/$sample_name/indels/large.indels | joinx sort >$output_dir/$sample_name/indels/indels.hq`;
  }
  my $indel_file = "$output_dir/$sample_name/indels/indels.hq";
  

  #SVs
  my $process_svs = $self->process_svs;
  my $sv_file = "$build_dir/validation/sv/assembly_output.csv.merged.readcounts.somatic.wgs_readcounts.somatic";
  if($process_svs){
      unless( -e $sv_file ){
          $self->warning_message("SV results for $sample_name not found, skipping SVs");
          $process_svs = 0;
      }
  }
  if($process_svs){
      `mkdir $output_dir/$sample_name/svs`;
      `cp $sv_file $output_dir/$sample_name/svs/svs.hq` unless( -e "$output_dir/$sample_name/svs/$sv_file");
  }


  
  #clean the files
  my $i=0;
  for($i=0;$i<@snv_files;$i++){
      $snv_files[$i] = cleanFile($snv_files[$i])
  }
  $indel_file = cleanFile($indel_file);

  #-------------------------------------------------
  #filter out the off-target regions, if target regions are available
  if($self->restrict_to_target_regions){
      $self->status_message("Filtering out off-target regions...");
      
      my $featurelist;
      if($self->target_regions) {
          $featurelist = $self->target_regions;
      } else { 
          $featurelist = $build->target_region_set->file_path;
      }


      if(defined($featurelist) && (-s $featurelist)){
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
          my $cmd;
          $cmd = "joinx sort '" . $output_dir . "/" . $sample_name . "/featurelist.tmp'" . " >'" . $output_dir . "/" . $sample_name . "/featurelist'";
          `$cmd`;
          $cmd = "rm -f '" . $output_dir . "/" . $sample_name . "/featurelist.tmp'";
          `$cmd`;

          #do the intersection
          my @new_snv_files;
          foreach my $file (@snv_files){
              my $new_file = addName($file,"ontarget");
              $cmd = "joinx intersect -a '" . $file . "' -b '" . $output_dir . "/" . $sample_name . "/featurelist' >'" . $new_file . "'";
              `$cmd`;
              $file = "$new_file";
              push(@new_snv_files,$file)
          }
          @snv_files = @new_snv_files;

          my $new_file = addName($indel_file,"ontarget");
          $cmd = "joinx intersect -a '" . $indel_file . "' -b '" . $output_dir . "/" . $sample_name . "/featurelist' >'" . $new_file . "'";
          `$cmd`;
          $indel_file = "$new_file";

      } else {
          $self->warning_message("feature list not found or target regions not specified; No target region filtering being done even though --restrict-to-target-regions set.");
      }
  }

  #-------------------------------------------------
  #remove filter sites specified by the user
  if(defined($self->filter_sites)){
      $self->status_message("Applying user-supplied filter...");
      my $filterSites = getFilterSites($self->filter_sites);

      for($i=0;$i<@snv_files;$i++){
          $snv_files[$i] = removeFilterSites($snv_files[$i],$filterSites);
      }
      $indel_file = removeFilterSites($indel_file,$filterSites);
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

      $self->status_message("Removing user-specified filter regions...");
      for($i=0;$i<@snv_files;$i++){
          my $snv_file = $snv_files[$i];
          my $cmd = "joinx intersect --miss-a $snv_file.filteredReg -a $snv_file -b $temp_file >/dev/null";
          my $result = Genome::Sys->shellcmd(
              cmd => "$cmd",
              );
          unless($result) {
              die $self->error_message("Failed to execute joinx: Returned $result");
          }
          $snv_files[$i] = addName($snv_files[$i],"filteredReg");
      }

      my $cmd = "joinx intersect --miss-a $indel_file.filteredReg -a $indel_file -b $temp_file >/dev/null";
      my $result = Genome::Sys->shellcmd(
          cmd => "$cmd",
          );
      unless($result) {
          die $self->error_message("Failed to execute joinx: Returned $result");
      }
      $indel_file = addName($indel_file,"filteredReg");
  }

  ##------------------------------------------------------
  # do annotation

  my $x=1;
  for($i=0;$i<@snv_files;$i++){
      $snv_files[$i] = doAnnotation($snv_files[$i], $annotation_build_name);
      my $x=1;
  }
  $indel_file = doAnnotation($indel_file, $annotation_build_name);


  #-------------------------------------------------------
  #add tiers
  if($self->add_tiers){
      $self->status_message("Adding tiers...");
      #do annotation
      for($i=0;$i<@snv_files;$i++){
          $snv_files[$i] = addTiering($snv_files[$i], $tiering_files);
          #convert back to annotation format (1-based)
          $snv_files[$i] = bedFileToAnnoFile($snv_files[$i]);
      }
      $indel_file = addTiering($indel_file, $tiering_files);
      #convert back to annotation format (1-based)
      $indel_file = bedFileToAnnoFile($indel_file);
  }


  #----------------------------------------------------
  # add dbsnp/gmaf

  if ($self->add_dbsnp_and_gmaf){
      $self->status_message("Adding dbsnp info...");

      if(defined($model->previously_discovered_variations_build)){
          my $dbsnpPath = $model->previously_discovered_variations_build->snv_result->path;
          my $dbsnpVcf = $dbsnpPath . "snvs.hq.vcf";
          if(-s $dbsnpVcf){
              for($i=0;$i<@snv_files;$i++){
                  my $snv_file = $snv_files[$i];
                  my $newfile = addName($snv_files[$i],"rsid");

                  #handle zero size files
                  if( -z $snv_file ){
                      `touch $newfile`;
                  } else {
                      my $db_cmd = Genome::Model::Tools::Annotate::AddRsid->create(
                          anno_file => $snv_file,
                          output_file => $newfile,
                          vcf_file => "$dbsnpVcf",
                          );
                      unless ($db_cmd->execute) {
                          die $self->error_message("Failed to add dbsnp anno to file $snv_file.");
                      }
                  }
                  $snv_files[$i] = $newfile;
              }
              #pad indel file with tabs to match - if we ever start annotating with indels from dbsnp, replace this section
              my $newfile = addName($indel_file,"rsid");
              open(OUTFILE,">$newfile");
              my $inFh = IO::File->new( "$indel_file" ) || die "can't open file2d\n";
              while( my $line = $inFh->getline ){
                  chomp($line);
                  if($line =~ /^chrom/){
                      print OUTFILE $line . "\trsid\tGMAF\n";
                      next;
                  }
                  print OUTFILE $line . "\t-\t-\n";
              }
              close($inFh);
              close(OUTFILE);
              $indel_file = $newfile;
          } else {
              $self->warning_message("couldn't find VCF with dbsnp annotations - skipping dbsnp anno");
          }
      } else {
	  $self->warning_message("Warning: couldn't find VCF with dbsnp annotations - skipping dbsnp anno");
      }
  }



  #-------------------------------------------------------
  #get readcounts
  my $indel_size=$self->max_indel_size;
  if($self->get_readcounts){
      $self->status_message("Getting readcounts...");
      if(!(defined($tumor_bam)) || !(-s $tumor_bam)){
          $self->warning_message("Tumor bam not found - skipping readcounting");
      } 
      else {
	  #getting readcounts for SNV files
          for($i=0;$i<@snv_files;$i++){
	      my @bamfiles=();
	      my @labels=();
	      if(scalar(@$som_var_bam_files)) { #if user specified som_var_model; grab readcounts from that first
		  push(@bamfiles, @$som_var_bam_files);
		  push(@labels,('somvar_Normal','somvar_Tumor')); #manually label the header #default is more suited for IGV xml crap
	      }

              if(!defined($normal_bam) || $self->tumor_only){
		  push(@bamfiles,$tumor_bam);
		  push(@labels,'val_Tumor');
              } elsif(defined($normal_bam) && !$self->tumor_only) {
		  push(@bamfiles,($normal_bam,$tumor_bam));
		  push(@labels,('val_Normal','val_Tumor'));
              }

	      $snv_files[$i] = getReadcounts($snv_files[$i], $ref_seq_fasta, \@bamfiles, \@labels, $bam_readcount_version,$indel_size)
	      
          }

	  #grab readcount for indels from the som_varmodel
	  if(scalar(@$som_var_bam_files)) { #if user specified som_var_model
	      my @bamfiles=();
	      my @labels=();
	      push(@bamfiles, @$som_var_bam_files);
	      push(@labels,('somvar_Normal','somvar_Tumor')); #manually label the header #default is more suited for IGV xml crap
	      $indel_file = getReadcounts($indel_file, $ref_seq_fasta, \@bamfiles,\@labels, $bam_readcount_version,$indel_size);

	  }

          #we can grab more accurate indel readcounts from the intermediate files, but only if the indel realignment was done
          #realigned
	  my $get_tumor_only;
	  if(!defined($normal_bam) || $self->tumor_only){
	      $get_tumor_only=1;
	  }elsif(defined($normal_bam) && !$self->tumor_only) {
	      $get_tumor_only=0;
	  }

          
          if($small_indel_file =~ /final_output/){
              my %counts;
              #first small indels
              my $inFh = IO::File->new( $small_indel_file ) || die "can't open small indel file\n";
              while( my $line = $inFh->getline ){
                  chomp($line);
                  my @F = split("\t",$line);

                  next if($line =~ /^chrom/ || !($line =~ /Somatic/));
                  my $key = join("\t",@F[0..4]);
                  $F[9] =~ s/%//g;
                  $F[13] =~ s/%//g;

		  if($get_tumor_only){
		      $counts{$key} = join("\t",@F[11..13]);  #only grab tumor readcount
		  }else {
		      $counts{$key} = join("\t",(@F[7..9],@F[11..13]));  #grab both normal and tumor read count
		  }
                  #$counts{$key} = join("\t",(@F[7..9],@F[11..13]));
              }
              #then large indels
              $large_indel_file =~ s/\.adapted//g;
              $inFh = IO::File->new( $large_indel_file ) || die "can't open large indel file\n";
              while( my $line = $inFh->getline ){
                  chomp($line);
                  my @F = split("\t",$line);
                  next if($line =~ /^chrom/ || !($line =~ /Somatic/));

                  my @comp = split("_",$F[0]);
                  my $key = $comp[0];
                  
                  if($comp[1] =~ /\d+\((\d+)\)/){
                      $key = $key . "\t" . $1;
                  } else {
                      $self->warning_message("unable to parse \"$F[0]\"");
                      next;
                  }
                  if($comp[2] =~ /\d+\((\d+)\)/){
                      $key = $key . "\t" . $1;
                  } else {
                      $self->warning_message("unable to parse \"$F[0]\"");
                      next;
                  }
                  $key = $key . "\t" . $comp[5] . "\t" . $comp[6];
		  $self->status_message($key);
                  if($get_tumor_only){
		      $counts{$key} = join("\t",($F[27]-$F[21], $F[21], $F[28]*100));
		  }else {
		      $counts{$key} = join("\t",($F[22]-$F[11], $F[11], $F[23]*100, $F[27]-$F[21], $F[21], $F[28]*100));
		  }

              }

	      my $header = join("\t",qw(val_Normal_ref_count val_Normal_var_count val_Normal_VAF val_Tumor_ref_count val_Tumor_var_count val_Tumor_VAF));
              open(OUTFILE,">$indel_file.rcnt");
	      $header = join("\t",qw(val_Tumor_ref_count val_Tumor_var_count val_Tumor_VAF)) if($get_tumor_only);

              $inFh = IO::File->new( $indel_file ) || die "can't open small indel file\n";
              while( my $line = $inFh->getline ){
                  chomp($line);
                  my @F = split("\t",$line);

                  if($line =~ /^chrom/){
                      print OUTFILE $line . "\t" . $header . "\n";
                      next;
                  }
                  my $key = join("\t",@F[0..4]);

                  if(defined($counts{$key})){
                      print OUTFILE join("\t",(@F)) . "\t" . $counts{$key} . "\n";
                  } else {
		      my $NA_string;
		      if($get_tumor_only) { 
			  $NA_string = "\tNA\tNA\tNA\n";
		      }else {
			  $NA_string = "\tNA\tNA\tNA\tNA\tNA\tNA\n";
		      }
		      print OUTFILE join("\t",(@F)) . $NA_string;
                  }
              }
	      if($indel_file =~ /\.rcnt$/) { #if a file name already ends in rcnt change the new file name to avoid *.rcnt.rcnt
		  Genome::Sys->rename("$indel_file.rcnt", $indel_file);
	      }else {
		  $indel_file = addName($indel_file,"rcnt"); #if a file name does not end in rcnt, append rcnt suffix
	      }

          } else { #not realigned, use bam-readcount
	      my @bamfiles=();
	      my @labels=();

              if(!defined($normal_bam) || $self->tumor_only){
		  push(@bamfiles,$tumor_bam);
		  push(@labels,'val_Tumor');
              } elsif(defined($normal_bam) && !$self->tumor_only) {
		  push(@bamfiles,($normal_bam,$tumor_bam));
		  push(@labels,('val_Normal','val_Tumor'));
              }
	      $indel_file = getReadcounts($indel_file, $ref_seq_fasta, \@bamfiles, \@labels, $bam_readcount_version, $indel_size);
          }
      }
  }

  #------------------------------------------------------
  # combine the files into one master table and append type (newcalls, validated, unvalidated)

  my $indel_type = "validated";
  if($small_indel_file =~ /variants/){
      $indel_type = "newcall";
  }
  open(OUTFILE,">$output_dir/$sample_name/snvs.indels.annotated");
  open(OUTFILE2,">$output_dir/$sample_name/snvs.indels.annotated.tmp");

  my $headerAdded=0;
  my $inFh = IO::File->new( $indel_file ) || die "can't open indel file\n";
  while( my $line = $inFh->getline ){
      chomp($line);
      my @F = split("\t",$line);

      if($line =~ /^chrom/){
          print OUTFILE $line . "\tstatus\n";
          $headerAdded=1;
          next;
      }
      print OUTFILE2 join("\t",(@F,$indel_type)) . "\n";
  }

  for($i=0;$i<@snv_files;$i++){
      my $snv_file = $snv_files[$i];
      my $snv_type = "newcall";
      if($snv_file =~ /\.validated\./){
          $snv_type = "validated";
      } elsif ($snv_file =~ /notvalidated/){
          $snv_type = "notvalidated";
      }

      $inFh = IO::File->new( $snv_file ) || die "can't open snv file\n";
      while( my $line = $inFh->getline ){
          chomp($line);
          if($line =~ /^chrom/){
              unless($headerAdded){
                  print OUTFILE $line . "\tstatus\n";
                  $headerAdded=1;
              }
              next;
          }

          my @F = split("\t",$line);
          print OUTFILE2 join("\t",(@F,$snv_type)) . "\n";
      }
  }
  close(OUTFILE);
  close(OUTFILE2);

  my $result = Genome::Sys->shellcmd(
      cmd => "joinx sort -i $output_dir/$sample_name/snvs.indels.annotated.tmp >>$output_dir/$sample_name/snvs.indels.annotated",
      );
  unless($result) {
      die $self->error_message("Failed to execute joinx: Returned $result");
  }
  $result = Genome::Sys->shellcmd(
      cmd => "rm -f $output_dir/$sample_name/snvs.indels.annotated.tmp"
      );
  unless($result) {
      die $self->error_message("Failed to execute rm: Returned $result");
  }

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
      $self->status_message("Generating Review files...");
      my @tiers = split(",",$self->tiers_to_review);
      my $tierstring = join("",@tiers);
      for my $i (@tiers){
          `grep -w tier$i $output_dir/$sample_name/snvs.indels.annotated >>$output_dir/$sample_name/snvs.indels.annotated.tier$tierstring.tmp`;
      }
      
      my @statuses = split(",",$self->statuses_to_review);
      foreach my $i (@statuses){
          `grep -w $i $output_dir/$sample_name/snvs.indels.annotated.tier$tierstring.tmp >>$output_dir/$sample_name/snvs.indels.annotated.tier$tierstring.tmp2`;
      }

      `head -n1 $output_dir/$sample_name/snvs.indels.annotated | grep "chromosome_name" >$output_dir/$sample_name/snvs.indels.annotated.tier$tierstring`;
      `joinx sort -i $output_dir/$sample_name/snvs.indels.annotated.tier$tierstring.tmp >>$output_dir/$sample_name/snvs.indels.annotated.tier$tierstring`;
      annoFileToSlashedBedFile("$output_dir/$sample_name/snvs.indels.annotated.tier$tierstring.tmp2","$output_dir/review/$sample_name.bed");
      `rm -f $output_dir/$sample_name/snvs.indels.annotated.tier$tierstring.tmp`;
      `rm -f $output_dir/$sample_name/snvs.indels.annotated.tier$tierstring.tmp2`;

      my $bam_files;
      my $labels;

      #add bam files from som-var model, if specified
      if(defined $self->somatic_variation_model_id){
	  push(@$bam_files,@$som_var_bam_files);
	  push(@$labels,@$som_var_labels);
      }

      
     # add bam files from this model
      if($self->tumor_only){
          push(@$bam_files, $tumor_bam);
          push(@$labels, "tumor $sample_name");
      } else {
          push(@$bam_files, $normal_bam);
          push(@$bam_files, $tumor_bam);

          push(@$labels, "normal $sample_name");
          push(@$labels, "tumor $sample_name");
      }

      my $igv_reference_name = "b37"; #default
      if(defined($self->igv_reference_name)){
          $igv_reference_name = $self->igv_reference_name;
      } else {
          $self->warning_message("No IGV reference name supplied - defaulting to build 37");
      }


      #create the xml file for review
      my $dumpXML = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
          bams => join(",",@$bam_files),
          labels => join(",",@$labels),
          output_file => "$output_dir/review/$sample_name.xml",
          genome_name => $sample_name,
          review_bed_files => "$output_dir/review/$sample_name.bed",
          reference_name => $igv_reference_name,
          );
      unless ($dumpXML->execute) {
          die $self->error_message("Failed to dump IGV xml.");
      }

      my $str = "\n--------------------------------------------------------------------------------\n"
              .  "Sites to review are here:\n"
              ."$output_dir/review/$sample_name.bed\n"
              ."IGV XML file is here:"
              ."$output_dir/review/$sample_name.xml\n\n";
      $self->status_message($str);
  }

  #------------------------------------------------
  # tar up the files to be sent to collaborators
  #
  if($self->create_archive){
      mkdir("$output_dir/$sample_name/$sample_name");

      chdir("$output_dir/$sample_name/");
      # #VCF files
      # if($self->include_vcfs_in_archive){
      #     if(-e "$build_dir/variants/indels.detailed.vcf.gz"){
      #         `ln -s $build_dir/variants/indels.detailed.vcf.gz $sample_name/indels.vcf.gz`;
      #     } elsif(-e "$build_dir/variants/indels.vcf.gz") {
      #         `ln -s $build_dir/variants/indels.vcf.gz $sample_name/indels.vcf.gz`;
      #     } else {
      #         print STDERR "WARNING: no indel VCF file available. If this is an older model, a rebuild may fix this\n";
      #     }
      #     if(-e "$build_dir/variants/snvs.annotated.vcf.gz"){
      #         `ln -s $build_dir/variants/snvs.annotated.vcf.gz  $sample_name/snvs.vcf.gz`;
      #     } elsif (-e "$build_dir/variants/snvs.vcf.gz"){
      #         `ln -s $build_dir/variants/snvs.vcf.gz  $sample_name/snvs.vcf.gz`;
      #     } else {
      #         print STDERR "WARNING: no snv VCF file available. If this is an older model, a rebuild may fix this\n";
      #     }
      # }
      #annotated snvs and indels
      `ln -s ../snvs.indels.annotated $sample_name/snvsAndIndels.annotated`;
      #same in excel format
      `ln -s ../snvs.indels.annotated.xls $sample_name/snvsAndIndels.annotated.xls`;
      #sv calls
      if($process_svs){
          `ln -s $sv_file $output_dir/$sample_name/$sample_name/svs`;
          #`ln -s $sv_file $output_dir/$sample_name/$sample_name/svs.annotated`;
      }
      #cnv calls - todo

      #tar it up
      `tar -chzvf $sample_name.tar.gz $sample_name`;
  }

  return 1;
}

sub indel_files {
    my ($self, $build_dir) = @_;
    my $small_indel_file = "$build_dir/validation/small_indel/final_output";
    unless( -e $small_indel_file && $self->use_assembled_indels ){
        $self->warning_message("realigned small indels not found (tumor only build?) or not requested. Using raw calls from validation data");
        $small_indel_file = "$build_dir/variants/indels.hq.bed";
    }
    my $large_indel_file = "$build_dir/validation/large_indel/combined_counts.csv.somatic.adapted";
    unless( -e $large_indel_file && $self->use_assembled_indels ){
        $self->warning_message("realigned large indels not found (tumor only build?) or not requested. Using raw calls from validation data");
        $large_indel_file = "$build_dir/variants/indels.hq.bed";
        $self->warning_message($large_indel_file);
    }
    $self->status_message("indel1: " . $large_indel_file);
    return ($small_indel_file,$large_indel_file);
}



sub process_somatic_variation_models {
#takes a somatic var model ID and returns the path to the 
#parent BAM files and labels

    my $self = shift;
    my $som_var_model_id = shift;
    my $sample_name = shift;   #som_var probably has its own sample name, but need to be consistent with validation model

    my $var_model = Genome::Model->get($som_var_model_id);


    my @bam_files;
    my @labels;
    my $som_var_obj={};
    if(!defined($var_model)){
	$self->error_message("Could not find a model with ID: $som_var_model_id");
    }
    else {
	my $tvar_build = $var_model->tumor_model->last_succeeded_build;
	my $nvar_build = $var_model->normal_model->last_succeeded_build;
	if(!defined($nvar_build) || !defined($tvar_build) ){
	    $self->error_message("Could not find a succeeded refalign builds from model ID: $som_var_model_id");
	}
	else {                  
	    my $tbam = $tvar_build->whole_rmdup_bam_file;
	    my $nbam = $nvar_build->whole_rmdup_bam_file;
	    if (!-s $tbam && !-s $nbam){
		$self->error_message("couldn't resolve bam files for somatic variation model");
	    } else {
		push(@bam_files, ($nbam,$tbam));
		push(@labels, "original normal $sample_name","original tumor $sample_name" );
	    }
	}
    }


    return (\@bam_files,\@labels);

}


1;
