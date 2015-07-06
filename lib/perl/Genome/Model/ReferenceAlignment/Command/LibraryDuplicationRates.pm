package Genome::Model::ReferenceAlignment::Command::LibraryDuplicationRates;

use strict;
use warnings;

class Genome::Model::ReferenceAlignment::Command::LibraryDuplicationRates {
    is => 'Command::V2',
    has => [
        builds => {
            is_many => 1,
            is => 'Genome::Model::Build::ReferenceAlignment',
            shell_args_position => 1,
            doc => 'List of builds to report on.',
        },
    ],
};

sub execute {
    my $self = shift;

    my $writer;
    for my $build ($self->builds) {
        my $library_metrics = $build->mark_duplicates_library_metrics_hash_ref;
        my @instrument_data = $build->instrument_data;
        for my $library_name (keys %{$library_metrics}) {
            my $data = $library_metrics->{$library_name};
            my @library_instrument_data = grep {$_->library_name eq $library_name} @instrument_data;
            $data->{BUILD_ID} = $build->id;
            $data->{MODEL_ID} = $build->model->id;
            $data->{SUBJECT_NAME} = $build->model->subject->name;
            $data->{COUNT_INSTRUMENT_DATA} = scalar(@library_instrument_data);
            unless ($writer) {
                my @headers = keys %{$data};
                $writer = Genome::Utility::IO::SeparatedValueWriter->create(
                    separator => "\t",
                    headers => \@headers,
                );
            }
            $writer->write_one($data);
        }
    }
    $writer->output->close;

    return 1;
}


1;
