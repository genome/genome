package Genome::Model::Tools::Annotate::ReannotateTranscriptVariants;

use strict;
use warnings;
use Genome;
use Command;
use IO::File;

class Genome::Model::Tools::Annotate::ReannotateTranscriptVariants {
    is => 'Command',
    has => [
        annotation_file => {
            is => 'FilePath',
            is_input => 1,
            is_optional => 0,
            doc => "File of variants. First five columns: chromosome_name start stop reference variant.",
        },
        output_file => {
            is => 'Text',
            is_input => 1,
            is_output=> 1,
            doc => "Store annotation in the specified file. Defaults to STDOUT if no file is supplied.",
            default => "STDOUT",
        },
        ],
    has_optional => [
        use_version => {
            is_input => 1,
            is => 'Text',
            default_value => __PACKAGE__->default_annotator_version,
            doc => 'Annotator version to use',
            valid_values => [0..4],
        },
        no_headers => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'Exclude headers in report output',
        },
        # Transcript params
        annotation_filter => {
            is => 'String',
            is_optional => 1,
            is_input => 1,
            default => 'gene',
            doc => 'The type of filtering to use on the annotation results. There are currently 3 valid options:
                    "none" -- This returns all possible transcript annotations for a variant. All transcript status are allowed including "unknown" status.
                    "gene" -- This returns the top transcript annotation per gene. This is the default behavior.
                    "top" -- This returns the top priority annotation for all genes. One variant in, one annotation out.',
        },
        reference_transcripts => {
            is => 'String',
            is_input => 1,
            is_optional => 1,
            doc => 'provide name/version number of the reference transcripts set you would like to use ("NCBI-human.combined-annotation/0").  Leaving off the version number will grab the latest version for the transcript set, and leaving off this option and build_id will default to using the latest combined annotation transcript set. Use this or --build-id to specify a non-default annoatation db (not both). See full help output for a list of available reference transcripts.'
        },
        ],
};

sub help_synopsis {
    return <<EOS
gmt annotate reannotate-transcript-variants --annotated-file variants.anno --output-file variants.newversion.anno --annotation-filter top --reference-transcripts NCBI-human.combined-annotation/70_38_v5 --user-version 4
EOS
}

sub default_annotator_version {
    return 4;
}

sub help_detail {
    my $self = shift;

    #Generate the currently available annotation models on the fly
    my @currently_available_models = Genome::Model->get(subclass_name => 'Genome::Model::ImportedAnnotation');
    my $currently_available_builds; 
    foreach my $model (@currently_available_models) {
        next unless $model;
        foreach my $build ($model->builds) {
            if($build and $build->status eq 'Succeeded' and $build->name =~ /ensembl/) {  #probably implicit in the loops, but in case we get undefs in our list
                 $currently_available_builds .= "\t" . $model->name . "/" . $build->version . "\n" if $build->version !~ /old/ and $model->name and $build->version;
            }
        }
    }

    return <<EOS
takes a annotated file, strips out the variants and reannotates using the specified reference transcripts and annotator version. Preserves any trailing fields

COMMONLY USED REFERENCE TRANSCRIPTS WITH VERSIONS
$currently_available_builds
EOS
}

############################################################

sub execute {
  my $self = shift;

  my %trailing_fields;

  #first, store all the extra fields and create a temprorary 5-col file for annotation
  my ($tifh,$temp_infile) = Genome::Sys->create_temp_file;
  my ($tofh,$temp_outfile) = Genome::Sys->create_temp_file;
  close($tofh);

  my $ifh = Genome::Sys->open_file_for_reading($self->annotation_file);
  while (my $line = $ifh->getline) {
      chomp $line;
      
      my ($chr, $start, $stop, $ref, $var, $type, $gene_name, $transcript_name, $transcript_species, $transcript_source, $transcript_version, $strand, $transcript_status, $trv_type, $c_position, $amino_acid_change, $ucsc_cons, $domain, $all_domains, $deletion_substructures, $transcript_error, @rest) = split(/\t/, $line);
      my $key = join("\t", ($chr, $start, $stop, $ref, $var));
      if(@rest){
          $trailing_fields{$key} = join("\t",@rest);
      }
  
      unless($chr =~ /chromosome_name/){
          print $tifh $key . "\n";
      }
  }
  close($tifh);


  my $anno_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
      variant_file => $temp_infile,
      output_file => $temp_outfile,
      reference_transcripts => $self->reference_transcripts,
      annotation_filter => $self->annotation_filter,
      use_version => $self->use_version,
      no_headers => $self->no_headers,
      );
  unless ($anno_cmd->execute) {
      die "Failed to annotate variants successfully.\n";
  }

  #match up the newly annotated sites with their trailing fields, if defined.
  my $ofh = Genome::Sys->open_file_for_writing($self->output_file);

  $ifh = Genome::Sys->open_file_for_reading($temp_outfile);
  while (my $line = $ifh->getline) {
      chomp $line;
      my ($chr, $start, $stop, $ref, $var, @rest) = split /\t/, $line;
      
      my $key = join("\t", ($chr, $start, $stop, $ref, $var));
      if(defined($trailing_fields{$key})){
          print $ofh $line . "\t" . $trailing_fields{$key} . "\n";
      } else {
          print $ofh $line . "\n";
      }
  }
  close($ifh);
  close($ofh);
  
  return 1;
}


