package Genome::ModelGroup::Command::GenerateCoverageReport;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::GenerateCoverageReport {
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
        },
    ],
};

sub help_brief {
    "Generate a copy of the coverage report",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome model convergence generate-coverage-report --group 12345 --output /tmp/report.html
EOS
}

sub help_detail {
    return <<EOS
Generates a copy of the coverage report as is produced by the coverage HTML view of ModelGroups.
EOS
}

sub execute {
    my $self = shift;

    my $output_file = $self->output_file;
    Genome::Sys->validate_file_for_writing($output_file); #This will croak on errors

    my $group = $self->group;
    my $view = $group->create_view(perspective => 'coverage', toolkit => 'html', xsl_root => Genome->base_dir. '/xsl' );
    my $content = $view->content;

    my $base = $self->base_url;
    $content =~ s!<head>!<head><base href="$base"/>!;

    Genome::Sys->write_file($self->output_file, $content);
    return 1;
}


1;
