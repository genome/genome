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

sub help_brief {
    'A command to print Picard MarkDuplicates duplication metrics for a reference alignment build.'
}

sub help_synopsis {
    'Dump reference alignment build duplication metrics to standard out in tab-delimited format.'
}

sub help_detail{
    return <<"EOS"
All reference alignment builds currently use Picard MarkDuplicates to flag PCR duplicates after merging the individual instrument data alignments.
One output from Picard MarkDuplicates is a metrics file that provides a per library breakdown of PCR and optical duplicate reads.  This tool simply parses the metrics file and dumps those metrics along with basic build/model information in tab-delimited form.
EOS
}

sub execute {
    my $self = shift;

    my @data;
    for my $build ($self->builds) {
        my @instrument_data = $build->instrument_data;
        my @library_names = List::MoreUtils::uniq(sort map { $_->library->name } @instrument_data);

        for my $library_name (@library_names) {
            my @library_instrument_data = grep {$_->library_name eq $library_name} @instrument_data;
            my $data = {
                'BUILD_ID' => $build->id,
                'MODEL_ID' => $build->model->id,
                'SUBJECT_NAME' => $build->model->subject->name,
                'COUNT_INSTRUMENT_DATA' => scalar(@library_instrument_data),
            };

             my @library_metrics = Genome::Model::Metric->get(
                build_id => $build->id,
                name => { operator => 'like', value => $library_name .'%' },
            );
            unless (@library_metrics) {
                $self->error_message('Missing MarkDuplicates metrics for library('. $library_name .') and build('. $build->id .')');
                die($self->error_message);
            }
            for my $library_metric (@library_metrics) {
                my $regex = $library_name .'_([A-Z_]+)';
                $library_metric->name =~ /^$regex/;
                my $key = $1;
                my $value = $library_metric->value;
                $data->{$key} = $value;
            }
            push @data, $data;
        }
    }
    unless (@data) {
        $self->warning_message('MarkDuplicates metric files were not found for builds!');
        return 1;
    }
    my @headers = List::MoreUtils::uniq(sort map { keys %$_ } @data);
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@headers,
    );
    for (@data) {
        $writer->write_one($_);
    }
    $writer->output->close;

    return 1;
}


1;
