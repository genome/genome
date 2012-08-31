package Genome::Model::Event::Build::AmpliconAssembly::Reports;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Event::Build::AmpliconAssembly::Reports {
    is => 'Genome::Model::Event',
};

sub execute {
    my $self = shift;

    for my $report_type ('stats', 'composition') {
        $self->_generate_and_save_report($report_type);
    }

    return 1;
}

sub _generate_and_save_report {
    my ($self, $type) = @_;

    my $class = 'Genome::Model::AmpliconAssembly::Report::'.Genome::Utility::Text::string_to_camel_case($type);
    my $generator = $class->create(
        build_id => $self->build->id,
    );
    unless ( $generator ) {
        $self->error_message(
            sprintf(
                'Could not create %s report generator (MODEL <Name:%s Id:%s> BUILD <Id:%s>)', 
                $type,
                $self->model->name,
                $self->model->id,
                $self->build->id,
            )
        );
        return;
    }
    my $report = $generator->generate_report;
    unless ( $report ) {
        $self->error_message(
            sprintf(
                'Could not generate %s report (MODEL <Name:%s Id:%s> BUILD <Id:%s>)', 
                $type,
                $self->model->name,
                $self->model->id,
                $self->build->id,
            )
        );
        return;
    }

    unless ( $self->build->add_report($report) ) {
        $self->error_message(
            sprintf(
                'Could not save %s report (MODEL <Name:%s Id:%s> BUILD <Id:%s>)', 
                $type,
                $self->model->name,
                $self->model->id,
                $self->build->id,
            )
        );
    }

    my @datasets = $report->get_datasets;
    unless ( @datasets ) { # not ok
        $self->error_message(
            sprintf(
                'No datasets in %s report (MODEL <Name:%s Id:%s> BUILD <Id:%s>)', 
                $type,
                $self->model->name,
                $self->model->id,
                $self->build->id,
            )
        );
        return;
    }
    my $file_base = sprintf(
        '%s/%s',
        $self->build->reports_directory,
        $report->name_to_subdirectory($report->name),
    );

    for my $dataset ( @datasets ) {
        my $dataset_name = $dataset->name;
        my $file = sprintf(
            '%s/%s.%s.tsv',
            $file_base,
            $self->model->subject_name,
            $dataset_name,
        );
        unlink $file if -e $file;
        my $fh = Genome::Sys->open_file_for_writing($file)
            or return; # not ok
        my ($svs) = $dataset->to_separated_value_string(separator => "\t");
        unless ( $svs ) { # not ok
            $self->error_message(
                sprintf(
                    'Could not get dataset (%s) for %s report (MODEL <Name:%s Id:%s> BUILD <Id:%s>)', 
                    $dataset_name,
                    $type,
                    $self->model->name,
                    $self->model->id,
                    $self->build->id,
                )
            );
            return;
        }
        $fh->print($svs);
        $fh->close;
    }

    return $report;
}

1;

#$HeadURL$
#$Id$
