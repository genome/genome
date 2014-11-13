package Genome::FeatureList::Command::Create;

use strict;
use warnings;

use Genome;


class Genome::FeatureList::Command::Create {
    is => 'Command::V2',
    has_input => [
        name => { is => 'Text', doc => 'The name of the feature-list' },
        format => { is => 'Text', doc => 'Indicates whether the file follows the BED spec.', valid_values => Genome::FeatureList->__meta__->property('format')->valid_values },
        file_path => { is => 'Text', doc => 'Path to the BED file on the file system (will be copied into an allocation)' },
        content_type => { is => 'Text', doc => 'the kind of information in the BED file', valid_values => Genome::FeatureList->__meta__->property('content_type')->valid_values },
    ],
    has_optional_input => [
        source => { is => 'Text', len => 64, doc => 'Provenance of this feature list. (e.g. Agilent)', },
        reference => { is => 'Genome::Model::Build::ReferenceSequence', doc => 'reference sequence build for which the features apply' },
        subject => { is => 'Genome::Model::Build', doc => 'subject to which the features are relevant' },
        description => { is => 'Text', doc => 'General description of the BED file' },
    ],
};

sub help_brief {
    "Create a new feature-list.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt feature-list create --name example-region --format true-BED --source WUGC --file-path path/to/file.bed
EOS
}

sub help_detail {
    return <<EOS 
Create a new feature-list.
EOS
}

sub execute {
    my $self = shift;

    my %create_params = (
        name => $self->name,
        format => $self->format,
    );

    for my $property (qw(source reference subject content_type description)) {
        my $value = $self->$property;
        $create_params{$property} = $value if $value; 
    }

    my $content_hash = Genome::Sys->md5sum($self->file_path);
    $create_params{file_path} = $self->file_path;
    $create_params{file_content_hash} = $content_hash;

    my $feature_list = Genome::FeatureList->create( %create_params );

    unless($feature_list) {
        $self->error_message('Failed to create feature-list.');
        return;
    }

    $self->status_message('Created feature-list "' . $feature_list->name . '" with ID: ' . $feature_list->id);
    return $feature_list;
}

1;
