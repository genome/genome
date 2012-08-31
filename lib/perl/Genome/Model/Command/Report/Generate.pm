package Genome::Model::Command::Report::Generate;
#:adukes check, further comments on where Report/Generator types for Model classes go

use strict;
use warnings;

use Genome;

require Genome::Model::Report;
use Data::Dumper 'Dumper';
use IO::Handle;

class Genome::Model::Command::Report::Generate {
    is => 'Genome::Model::Command::Report',
    has_optional => [
    directory => {
        is => 'Text',
        doc => 'Directory to save the report to.  If not given, the report will added to the build.',
    },
    force => {
        is => 'Boolean',
        default_value => 0,
        doc => 'If the report exists, overwrite it with newly generated report',
    },
    ],
};

#< Helps >#
sub help_brief {
    'Generate a report';
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
EOS
}

#< Command >#
sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    if ( $self->directory ) {
        unless ( Genome::Sys->validate_directory_for_write_access($self->directory) ) {
            $self->error_message('Requested to save to directroy XSLT file is required to create');
            return;
        }
    }

    return $self;
}

sub execute {
    my $self = shift;
    
    my $generator_class = Genome::Model::Report::get_report_class_for_type_name_and_report_name(
        $self->model->type_name,
        $self->report_name,
    ) or return;

    my $generator = $generator_class->create(
        build_id => $self->build_id,
    );

    my $report = $generator->generate_report;
    unless ( $report ) {
        $self->error_message("Can't generate report ");
        return;
    }

    return $self->_save_report($report);
}

sub _save_report {
    my ($self, $report) = @_;

    if ( $self->directory ) {
        unless ( $report->save($self->directory, $self->force) ) {
            $self->error_message( sprintf('Could not save report to directory (%s)', $self->directory) );
            return;
        }
    }
    else {
        unless ( $self->build->add_report($report) ) {
            $self->error_message( sprintf('Could not add report to build (%s)', $self->build_id) );
            return;
        }
    }

    return 1;
}

1;

