package Genome::Model::SomaticVariation::Command::SimplifyVcf;
#Written by Malachi Griffith

use strict;
use warnings;
use Genome;
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use File::Path;

class Genome::Model::SomaticVariation::Command::SimplifyVcf {
    is => 'Command::V2',
    has_input => [
        builds => { 
              is => 'Genome::Model::Build::SomaticVariation',
              is_many => 1,
              shell_args_position => 1,
              require_user_verify => 0,
              doc => 'Somatic variation build to summarize',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written', 
        },
        include_and_merge_indels => {
            is => 'Boolean',
            default => 0,
            doc => 'If set to true, simplify the indel file as well and merge the indels with snvs',
        },
    ],
    doc => 'Create a simplified VCF that contains only column for tumor/normal and only passing somatic SNVs',
};

sub help_synopsis {
    return <<EOS

genome model somatic-variation simplify-vcf  --outdir=/tmp/  129399189

genome model somatic-variation simplify-vcf  --outdir=/tmp/  id=129399189

genome model somatic-variation simplify-vcf  --outdir=/tmp/  model.id=2888675853

genome model somatic-variation simplify-vcf  --outdir=/tmp/  'id in [129399189,129396794]'

EOS
}

sub help_detail {
    return <<EOS

Create a VCF file that is a simplified version of the default put out by the pipeline.  

Remove all variants except those that pass: based on the FT value of the *tumor* sample column

EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless (-e $self->outdir && -d $self->outdir) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['outdir'],
	                                          desc => RED . "Outdir: " . $self->outdir . " not found or not a directory" . RESET,
                                          );
  }
  return @errors;
}

sub execute {
    my $self = shift;
    my @somvar_builds = $self->builds;
    my $outdir = $self->outdir;

    unless ($outdir =~ /\/$/){
        $outdir .= "/";
    }

    #Go through each somatic build and create XML session files
    my $somvar_build_count = scalar(@somvar_builds);

    foreach my $somvar_build (@somvar_builds){
        my $somvar_build_id = $somvar_build->id;
        my $subject = $somvar_build->subject;
        my $subject_name = $subject->name;
        my $subject_common_name = $subject->common_name;
        my $subject_tissue_desc = $subject->tissue_desc || "NULL";
        my $patient_common_name = $subject->individual_common_name || "NULL";

        $self->debug_message("\nProcessing somatic-variation build: $somvar_build_id ($patient_common_name, $subject_name, $subject_tissue_desc, $subject_common_name)");

        #Get the tumor and normal sample names
        my $tumor_sample = $somvar_build->tumor_model->subject->name;
        unless($tumor_sample){
            $self->error_message("Unable to resolve tumor sample name from build");
            return;
        }
        my $normal_sample = $somvar_build->normal_model->subject->name;
        unless($normal_sample){
            $self->error_message("Unable to resolve normal sample name from build");
            return;
        }

        for my $type ("snv", "indel") {
            if (!$self->include_and_merge_indels and $type eq "indel") {
                next;
            }
            #Find the VCF SNV file in the build dir
            my $find_method = "find_" . $type . "_vcf_file";
            my $original_vcf_path = $self->$find_method('-build'=>$somvar_build);

            #Name the output file in a human readable way ($patient_common_name . $subject_name . $somatic_variation_build_id . passing.somatic.snvs.vcf
            my $resolve_method = "resolve_" . $type . "_vcf_filename";
            $DB::single=1;
            my $new_vcf_path = $self->$resolve_method($somvar_build);

            #Parse the VCF file and dump simplified version
            $self->simplify_vcf('-infile'=>$original_vcf_path, '-outfile'=>$new_vcf_path, '-tumor_sample'=>$tumor_sample, '-normal_sample'=>$normal_sample);
        }

        if ($self->include_and_merge_indels) {
            $self->merge_simplified_vcfs($somvar_build);
        }
    }

    return 1;
}

sub merge_simplified_vcfs {
    my ($self, $build) = @_;

    my $snv_vcf = $self->resolve_snv_vcf_filename($build);
    my $indel_vcf = $self->resolve_indel_vcf_filename($build);
    my $output_vcf = $self->resolve_combined_vcf_filename($build);

    my $joinx_command = Genome::Model::Tools::Joinx::VcfMerge->create(
        input_files => [$snv_vcf, $indel_vcf],
        output_file => $output_vcf,
        merge_samples => 1,
        use_bgzip => 0,
    );

    unless ($joinx_command->execute) {
        die $self->error_message("Failed to run joinx vcf merge");
    }

    $self->debug_message("$snv_vcf and $indel_vcf merged into $output_vcf");

    return 1;
}

sub resolve_indel_vcf_filename {
    my ($self, $build) = @_;
    return $self->resolve_vcf_filename($build, 'indels');
}

sub resolve_snv_vcf_filename {
    my ($self, $build) = @_;
    return $self->resolve_vcf_filename($build, 'snvs');
}

sub resolve_combined_vcf_filename {
    my ($self, $build) = @_;
    return $self->resolve_vcf_filename($build, 'combined');
}

sub resolve_vcf_filename{
  my ($self, $build, $type) = @_;

  my $build_id = $build->id;
  my $subject = $build->subject;
  my $subject_name = $subject->name;
  my $subject_common_name = $subject->common_name;
  $subject_name=~s/\(/_/g; $subject_name=~s/\)/_/g;
  my $subject_tissue_desc = $subject->tissue_desc;
  my $patient_common_name = $subject->individual_common_name;

  #Name the output file in a human readable way ($patient_common_name . $subject_name . $somatic_variation_build_id . passing.somatic.snvs.vcf
  my @names;
  push(@names, $patient_common_name) if ($patient_common_name);
  push(@names, $subject_name) if ($subject_name);
  push(@names, $build_id);
  push(@names, "passing.somatic." . $type . ".vcf");
  my $filename = join("_", @names);

  $self->debug_message("Writing simplified VCF to file: $filename");
  return $self->outdir . "/$filename";
}

sub find_snv_vcf_file {
    my ($self, %args) = @_;
    my $build = $args{'-build'};
    return $self->find_vcf_file('-build'=>$build,'-type'=>'snv');
}

sub find_indel_vcf_file {
    my ($self, %args) = @_;
    my $build = $args{'-build'};
    return $self->find_vcf_file('-build'=>$build,'-type'=>'indel');
}

sub find_vcf_file {
    my ($self, %args) = @_;
    my $build = $args{'-build'};
    my $type = $args{'-type'};

    my $data_dir = $build->data_directory;
    my $vcf_file;

    if (-e $data_dir . "/variants/" . $type . "s.annotated.vcf.gz"){
        $vcf_file = $data_dir . "/variants/" . $type . "s.annotated.vcf.gz";
    }elsif(-e $data_dir . "/variants/" . $type . "s.vcf.gz"){
        $vcf_file = $data_dir . "/variants/" . $type ."s.vcf.gz";
    }else{
        $self->error_message("Could not find VCF file in $data_dir");
        return;
    }
    $self->debug_message("Found VCF file: $vcf_file");
    return $vcf_file;  
}

sub simplify_vcf{
  my $self = shift;
  my %args = @_;
  my $infile = $args{'-infile'};
  my $outfile = $args{'-outfile'};
  my $tumor_sample_name = $args{'-tumor_sample'};
  my $normal_sample_name = $args{'-normal_sample'};

  open (VCF1, "zcat $infile |") || die "\n\nCould not open input VCF file: $infile\n\n";
  open (VCF2, ">$outfile") || die "\n\nCould not open output VCF file: $outfile\n\n";
  my %header;
  my $header_found = 0;
  my $variant_count = 0;
  my $passing_variant_count = 0;
  while(<VCF1>){
    chomp($_);
    my $record = $_;
    my @line = split("\t", $record);

    #Deal with all header lines
    if ($record =~ /^\#/){
      #Get the column header line info
      if ($line[0] =~ /CHROM/){
        $header_found = 1;
        my $p = 0;
        foreach my $head (@line){
          $header{$head}{pos} = $p;
          $p++;
        }
        unless ($header{$tumor_sample_name}){
          $self->error_message("Could not find tumor column for sample: $tumor_sample_name");
          return;
        }
        unless ($header{$normal_sample_name}){
          $self->error_message("Could not find normal column for sample: $normal_sample_name");
          return;
        }
        unless ($header{'FORMAT'}){
          $self->error_message("Could not find FORMAT column in header line");
          return;
        }
        #Print a header with only one column for normal and tumor
        print VCF2 "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[$header{$normal_sample_name}{pos}]\t$line[$header{$tumor_sample_name}{pos}]\n";
        next;
      }else{
        print VCF2 "$record\n";
        next;
      }
    }

    #Filter out entries that do not pass all filters in the tumor
    my $format_string = $line[$header{'FORMAT'}{pos}];
    my @format = split(":", $line[$header{'FORMAT'}{pos}]);
    my $tumor_string = $line[$header{$tumor_sample_name}{pos}];
    my @tumor = split(":", $line[$header{$tumor_sample_name}{pos}]);
    my $normal_string = $line[$header{$normal_sample_name}{pos}];
    my @normal = split(":", $line[$header{$normal_sample_name}{pos}]);
    my $p;
    my $ft_pos;
    foreach my $format_tag (@format){
      $ft_pos = $p if ($format_tag eq "FT");
      $p++;
    }
    unless (defined($ft_pos)){
      #TODO: If 'FT' is not defined lq/hq status of a variant can not be determined.  
      $self->error_message("Could not determine position of FT tag in FORMAT column:\n@line");
      die($self->error_message);
    }
    my $tumor_filter_value = $tumor[$ft_pos];
    $variant_count++;
    unless ($tumor_filter_value){
      print Dumper @line;
    }
    next unless ($tumor_filter_value eq "PASS");
    $passing_variant_count++;
    print VCF2 "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[$header{$normal_sample_name}{pos}]\t$line[$header{$tumor_sample_name}{pos}]\n";

  }
  close(VCF1);
  close(VCF2);

  $self->debug_message("Found " . $variant_count . " total variants, of which " . $passing_variant_count . " pass all filters in the tumor");

  unless ($header_found){
    $self->error_message("Could not find VCF header line in file $infile");
    return;
  }

  return;
}


1;

