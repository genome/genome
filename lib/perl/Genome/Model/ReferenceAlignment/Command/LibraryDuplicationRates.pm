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
All reference alignment builds with duplicate marking currently use Picard MarkDuplicates to flag PCR duplicates after merging the individual instrument data alignments.
One output from Picard MarkDuplicates is a metrics file that provides a per library breakdown of PCR and optical duplicate reads.  This tool returns the database stored build metrics or parses the metric files in the alignment directory. The tab-delimited output also includes basic model and build information.
EOS
}

sub execute {
    my $self = shift;

    my @data;
    for my $build ($self->builds) {
        my @instrument_data = $build->instrument_data;
        my %library_counts;
        for my $instrument_data (@instrument_data) {
            $library_counts{$instrument_data->library_name}++;
        }

        for my $library_name (sort keys %library_counts) {
            my $data = {
                'BUILD_ID' => $build->id,
                'MODEL_ID' => $build->model->id,
                'SUBJECT_NAME' => $build->model->subject->name,
                'COUNT_INSTRUMENT_DATA' => $library_counts{$library_name},
            };
            my $library_metrics = $self->_resolve_library_metrics_from_build_metrics($build,$library_name);
            unless ($library_metrics) {
                $library_metrics = $self->_resolve_library_metrics_from_build_files($build,$library_name);
                unless ($library_metrics) {
                    $self->error_message('Failed to resolve the MarkDuplicate metrics for library \''. $library_name. '\' and build \''. $build->id .'\'');
                    die($self->error_message);
                }
            }
            for my $key (keys %{$library_metrics}) {
                $data->{$key} = $library_metrics->{$key};
            }
            push @data, $data;
        }
    }
    unless (@data) {
        $self->warning_message('MarkDuplicate metric files were not found for builds!');
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

sub _resolve_library_metrics_from_build_metrics {
    my $self = shift;
    my ($build,$library_name) = @_;

    my @library_metrics = Genome::Model::Metric->get(
        build_id => $build->id,
        name => { operator => 'like', value => $library_name .'%' },
    );
    unless (@library_metrics) {
        $self->warning_message('Missing MarkDuplicate database metrics for library('. $library_name .') and build('. $build->id .')');
        return;
    }
    my $library_metrics;
    for my $library_metric (@library_metrics) {
        my $regex = $library_name .'_([A-Z_]+)';
        $library_metric->name =~ /^$regex/;
        my $key = $1;
        my $value = $library_metric->value;
        $library_metrics->{$key} = $value;
    }
    return $library_metrics;
}

sub _resolve_library_metrics_from_build_files {
    my $self = shift;
    my ($build,$library_name) = @_;

    my $all_library_metrics = $build->mark_duplicates_library_metrics_hash_ref;
    unless ($all_library_metrics) {
        $self->warning_message('Missing MarkDuplicate metric files for library('. $library_name .') and build('. $build->id .')');
        return;
    }
    return $all_library_metrics->{$library_name};
}


1;
