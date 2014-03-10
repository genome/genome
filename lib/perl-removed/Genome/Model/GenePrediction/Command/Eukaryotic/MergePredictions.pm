package Genome::Model::GenePrediction::Command::Eukaryotic::MergePredictions;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::GenePrediction::Command::Eukaryotic::MergePredictions {
    is => 'Command',
    has => [
        temp_prediction_directories => {
            is => 'ARRAY',
            is_input => 1,
            doc => 'An array of temporary directories containing prediction files that need merged',
        },
        prediction_directory => {
            is => 'Path',
            is_input => 1,
            is_output => 1,
            doc => 'The directory that all predictions should be merged into',
        },
    ],
};

sub help_brief {
    return "Merges many gene predictions in temp directories into a single directory";
}

sub help_synopsis {
    return "Merges many gene predictions in temp directories into a single directory";
}

sub help_detail {
    return "Merges many gene predictions in temp directories into a single directory";
}

sub prediction_types {
    my $self = shift;
    return qw/ 
        Genome::Prediction::RNAGene 
        Genome::Prediction::Exon 
        Genome::Prediction::Transcript 
        Genome::Prediction::CodingGene 
        Genome::Prediction::Protein 
    /;
}

sub execute {
    my $self = shift;

    my $low = 20_000;
    my $high = 100_000;
    UR::Context->object_cache_size_lowwater($low);
    UR::Context->object_cache_size_highwater($high);

    $self->debug_message("Merging predictions into files in " . $self->prediction_directory);

    # Get meta object for each prediction type, grab attributes of the object (except for directory)
    TYPE: for my $type ($self->prediction_types) {
        $self->debug_message("*** Working on $type ***");
        my $meta = $type->__meta__;
        my @attributes = map { $_->property_name } $meta->_legacy_properties;
        @attributes = grep { $_ ne 'directory' } @attributes;

        TEMP_DIR: for my $temp_dir (@{$self->temp_prediction_directories}) {
            unless (-d $temp_dir) {
                $self->warning_message("No directory found at $temp_dir, this is most distressing!");
                next TEMP_DIR;
            }

            $self->debug_message("Looking at predictions in $temp_dir");
            my @temp_objects = $type->get(
                directory => $temp_dir,
            );
            next TEMP_DIR unless @temp_objects;

            # Create new objects in the merge directory that are exactly the same as the originals... except for directory!
            TEMP_OBJECT: for my $temp_object (@temp_objects) {
                my %properties;
                map { $properties{$_} = $temp_object->{$_} } @attributes;
                my $object = $type->create(
                    %properties,
                    directory => $self->prediction_directory,
                );
                unless (defined $object) {
                    confess "Could not create object of type $type with attributes:\n" . 
                        join("\n", map { join(",", $_, $properties{$_}) } @attributes );
                }
            }

            # Commit now so cache can be pruned later
            my $rv = UR::Context->commit;
            confess "Could not commit changes for $type!" unless defined $rv and $rv;
        }

        $self->debug_message("Done with $type!");
    }

    $self->debug_message("Successfully merged all files!");
    return 1;
}

1;

