package Genome::ModelGroup::Command::GenerateReport;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::GenerateReport {
    is => 'Genome::Command::Base',
    has => [
        group => {
            shell_args_position => 1,
            is => 'Genome::ModelGroup',
        },
        output_file => {
            is => 'String',
            doc => 'Location to write the report',
        },
        base_url => {
            is => 'Text',
            doc => 'Server with CSS and JS resources',
            default_value => $ENV{GENOME_SYS_SERVICES_SEARCH_URL},
        },
        perspective => {
            is => 'Text',
            default_value => 'coverage',
            valid_values => ['coverage', 'status'],
        },
    ],
};

sub help_brief {
    "Generate a copy of the coverage report or status view for a model-group",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome model convergence generate-report --group 12345 --output /tmp/report.html --base-url $ENV{GENOME_SYS_SERVICES_FILES_URL}
EOS
}

sub help_detail {
    return <<EOS
Generates a copy of the desired report as is produced by the HTML view of ModelGroups.
EOS
}

sub execute {
    my $self = shift;

    my $output_file = $self->output_file;
    Genome::Sys->validate_file_for_writing($output_file); #This will croak on errors

    my $group = $self->group;
    my $view = $group->create_view(perspective => $self->perspective, toolkit => 'html', xsl_root => Genome->base_dir. '/xsl' );
    my $content = $view->content;

    my $base = $self->base_url;
    $content =~ s!<head>!<head><base href="$base"/>!;

    Genome::Sys->write_file($self->output_file, $content);
    return 1;
}


1;
