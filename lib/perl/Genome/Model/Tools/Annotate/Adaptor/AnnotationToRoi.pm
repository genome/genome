package Genome::Model::Tools::Annotate::Adaptor::AnnotationToRoi;

use strict;
use warnings;
use Genome;
use IO::File;
use Carp 'confess';

class Genome::Model::Tools::Annotate::Adaptor::AnnotationToRoi {
    is => 'Command',
    has => [
        annotation_file => {
            is => 'Path',
            doc => 'Path to annotation output file to be parsed, can have extra columns that will be ignored',
        },
        roi_file => {
            is => 'Path',
            doc => 'Path to the ROI output file, format will be chromosome, start, stop (tab delimited)',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'Reference sequence build whose coordinates were used in the annotation file',
        }
    ],
};

sub help_brief {
    return "Figures out transcript start and stop positions from an annotation output file";
}

sub help_detail {
    return <<EOS
Given an annotation output file, this tool will pull out the transcript name and version and get 
more information (for now, chromosome, start, and stop) that will be written to an ROI file.
EOS
}

sub annotation_columns {
    my @cols = qw/
        chromosome_name
        start
        stop
        reference
        variant
        type
        gene_name
        transcript_name
        species
        transcript_source
        transcript_version
        strand
        transcript_status
        trv_type
        c_position
        amino_acid_change
        conservation
        domain
        all_domains
        deletions_structures
        transcript_error
    /;
    return \@cols;
}

sub execute {
    my $self = shift;
    my $reference_build_id = $self->reference_build->id;
    
    unless (-e $self->annotation_file) {
        confess "No annotation file found at " . $self->annotation_file;
    }
    unless (-s $self->annotation_file) {
        confess "Annotation file has no size at " . $self->annotation_file;
    }

    if (-e $self->roi_file) {
        $self->warning_message("File already exists at ROI output file location " . $self->roi_file . ", removing it!");
        unlink $self->roi_file;
    }

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        headers => $self->annotation_columns,
        separator => "\t",
        is_regex => 1,
        ignore_extra_columns => 1,
        input => $self->annotation_file,
    );
    confess 'Could not create separated value reader object for annotation file!' unless $reader;

    my $roi_fh = IO::File->new($self->roi_file, "w");
    confess 'Could not get file handle for ' . $self->roi_file unless defined $roi_fh;

    my ($version, $species);
    my @annotation_data_dirs;

    LINE: while (my $line = $reader->next) {
        # Need to figure out if headers have been included in the file and, if so, skip the first line
        if ($reader->line_number == 1) {
            my $first_column = $self->annotation_columns->[0];
            if ($line->{$first_column} eq $first_column) {
                $self->debug_message("Skipping first line due to headers.");
                next LINE;
            }
        }

        # If version and species are not defined, need to figure that out and get data directories
        unless (defined $version and defined $species) {
            $version = $line->{transcript_version};
            $species = $line->{species};
            @annotation_data_dirs = $self->annotation_directories_for_version_and_species($version, $species);
            confess 'Could not determine annotation data directories!' unless @annotation_data_dirs;
            $self->debug_message("Found annotation data for version $version and species $species!");
        }

        # Get the transcript using information in the output file
        my @transcripts = Genome::Transcript->get(
            data_directory => \@annotation_data_dirs,
            transcript_name => $line->{transcript_name},
            reference_build_id => $reference_build_id,
        );
        confess 'No transcripts found with name ' . $line->{transcript_name} unless @transcripts;
        confess 'Multiple transcripts found with name ' . $line->{transcript_name} if @transcripts > 1;

        my $transcript = shift @transcripts;
        $roi_fh->print(join("\t", $transcript->chrom_name, $transcript->transcript_start, $transcript->transcript_stop) . "\n");    
    }

    $roi_fh->close;
    $self->debug_message("ROI file containing chromosome, start, and stop of transcripts successfully created!");
    return 1;
}
     
# Finds the annotation data directories given a version and species
sub annotation_directories_for_version_and_species {
    my ($self, $version, $species) = @_;
    confess 'Not given version!' unless defined $version;
    confess 'Not given species!' unless defined $species;

    my $model_name = "NCBI-$species.combined-annotation";
    my $model = Genome::Model->get(name => $model_name);
    confess "Could not get model with name $model_name" unless $model;

    my $build = $model->build_by_version($version);
    confess "Could not get build with version $version from model $model_name" unless $build;

    return $build->determine_data_directory;
}

1;

