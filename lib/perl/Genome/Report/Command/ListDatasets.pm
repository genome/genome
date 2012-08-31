package Genome::Report::Command::ListDatasets;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use IO::Handle;

class Genome::Report::Command::ListDatasets {
    is => 'Genome::Report::Command',
};

#< Helps >#
sub help_brief {
    'List a report\'s datasets';
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
sub execute {
    my $self = shift;

    # report
    my $report = $self->report;
    unless ( $report ) {
        $self->error_message(
            sprintf('Report name (%s) not found for build (%s)', $self->report_name, $self->build_id)
        );
        return;
    }

    # names
    my @names = $report->get_dataset_names;
    unless ( @names ) {
        $self->error_message('No datasets found in report ('.$report->name.')');
        return;
    }

    return print join("\n", @names)."\n";
}

1;

#$HeadURL$
#$Id$
