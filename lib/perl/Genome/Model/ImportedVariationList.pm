package Genome::Model::ImportedVariationList;

use strict;
use warnings;

use Genome;

class Genome::Model::ImportedVariationList {
    is => 'Genome::ModelDeprecated',
    has => [
        server_dispatch => {
            is_constant => 1,
            is_class_wide => 1,
            value => 'inline',
            doc => 'lsf queue to submit the launcher or \'inline\''
        },
        job_dispatch => {
            is_constant => 1,
            is_class_wide => 1,
            value => 'inline',
            doc => 'lsf queue to submit jobs or \'inline\' to run them in the launcher'
        },
        reference_id => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'reference', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence' ],
            is_many => 0,
            is_mutable => 1,
            doc => 'reference sequence to align against'
        },
        reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            id_by => 'reference_id',
        },
        source_name => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'source_name', value_class_name => 'UR::Value::Text' ],
            is_many => 0,
            is_mutable => 1,
            is_optional => 1,
            doc => 'The name of the source of the imported variants (e.g., dbsnp, 1kg-wgs)',
        }
    ],
};

sub dbsnp_model_for_reference {
    my ($class, $reference) = @_;
    my @models = $class->get(name => "dbSNP-" . $reference->name);
    if (!@models and defined $reference->coordinates_from) {
        @models = $class->get(name => "dbSNP-" . $reference->coordinates_from->name);
    }
    return if @models != 1 or !$models[0]->reference->is_compatible_with($reference);
    return $models[0];
}

sub dbsnp_build_for_reference {
    my $self = shift;
    my $reference = shift;
    my $id;
    if ($reference->id eq "102671028" or $reference->id eq "106942997") {
        $id = 127786607;
    }
    elsif ($reference->id eq "107494762" or $reference->id eq "104420993") {
        $id = 127031265;
    }
    elsif ($reference->id eq "101947881") {
        $id = 106227442;
    }
    return Genome::Model::Build->get($id);
}

sub _execute_build {
    my ($self, $build) = @_;
    unless($build->model) {
        $self->error_message("Couldn't find model for build id " . $build->build_id . ".");
        return;
    }
    $self->status_message("Done.");
    return 1;
}

sub _resolve_disk_group_name_for_build {
    return $ENV{GENOME_DISK_GROUP_REFERENCES};
}

1;

