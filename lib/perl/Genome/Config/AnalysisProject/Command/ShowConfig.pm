package Genome::Config::AnalysisProject::Command::ShowConfig;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::ShowConfig {
    is => 'Command::V2',
    has_input => [
       analysis_project  => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis projects to put on hold',
            shell_args_position => 1,
        }
    ],
};

sub help_brief {
    return 'print out a summary of the configuration profile of an analysis project';
}

sub help_synopsis {
    return "genome config analysis-project show-config <analysis-project>";
}

sub help_detail {
    return <<"EOS"
This will print out a summary of the current configuration assigned to an analysis project
EOS
}

sub execute {
    my $self = shift;

    print $self->_header_line;

    for my $profile_item ($self->analysis_project->config_items) {
        print $self->_format_profile_item($profile_item);
    }

    return 1;
}

sub _header_line {
    my $self = shift;

    return sprintf($self->_line_template(),
        'File Path', 'Concrete',
        'Menu Item', 'Menu Item Name');
}

sub _format_profile_item {
    my $self = shift;
    my $item = shift;

    return sprintf($self->_line_template(100),
        $item->file_path,
        $self->_format_boolean($item->_is_concrete),
        $self->_format_boolean($item->analysis_menu_item),
        $item->analysis_menu_item ? $item->analysis_menu_item->name : '');

}

sub _format_boolean {
    my $arg = shift;
    return $arg ? '1' : '0';
}

sub _line_template {
    my $self = shift;
    my $padding = shift || 0;

    return "%-${padding}s\t%s\t%s\t%s\n";
}

1;
