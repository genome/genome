package Genome::Test::Factory::Model::ImportedAnnotation;
use Test::MockObject;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::ProcessingProfile::ImportedAnnotation;
use Genome::Test::Factory::Taxon;

our @required_params = qw(subject);

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::ImportedAnnotation->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::ImportedAnnotation->setup_object();
    return $p->id;
}

sub create_subject {
    return Genome::Test::Factory::Taxon->setup_object();
}

sub create_mock_build {
    my ($self, $model) = @_;

    $model = $self->setup_object if not $model;

    my $build = Test::MockObject->new;
    $build->set_always('model', $model);

    my $existing_file = Genome::Sys->write_file(Genome::Sys->create_temp_file_path('existing'), 1);
    my $existing_file_no_ref_build = Genome::Sys->write_file(Genome::Sys->create_temp_file_path('existing-no-ref-build'), 1);

    for my $file_method ( Genome::InstrumentData::AlignmentResult::Command::PicardRnaSeqMetrics->required_files_from_annotation_build) {
        my $existing_file_method = $file_method.'_exists';
        $build->mock(
            $existing_file_method,
            sub{ $build->mock($file_method, sub{ $existing_file } ); },
        );
        $build->$existing_file_method; # set to existing
        $build->mock(
            $file_method.'_does_not_exist',
            sub{ $_[0]->mock($file_method, sub{ '/dev/null' }); },
        );
        $build->mock(
            $file_method.'_exists_without_reference_build',
            sub{ $_[0]->mock($file_method, sub{ return if defined $_[2]; $existing_file_no_ref_build; } );
            },
        );
    }

    return $build;
}

1;
