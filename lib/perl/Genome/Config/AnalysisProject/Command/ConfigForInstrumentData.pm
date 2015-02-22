package Genome::Config::AnalysisProject::Command::ConfigForInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::ConfigForInstrumentData {
    is => ['Genome::Command::WithColor'],
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            shell_args_position => 1,
            is_many => 1,
        },
        color => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
            doc => 'Display report in color.'
        },
    ],
    doc => 'compare instdata to AnP configs to show matches, models and mismatches',
};

sub execute {
    my $self = shift;
    my @instrument_data = $self->instrument_data;

    for my $id (@instrument_data) {
        my @analysis_projects = $id->analysis_projects;
        unless(@analysis_projects) {
            $self->warning_message('No analysis-projects associated with instrument data %s.', $id->__display_name__);
            return 1;
        }

        for my $ap (@analysis_projects) {
            $self->_display_matches($id, $ap);
        }
    }

    return 1;
}


sub _bridge_status_colors {
    return (
        # InstrumentDataBridge statuses
        failed    => "red",
        new       => "white",
        processed => "green",
        skipped   => "magenta",
    );
}

sub _display_matches {
    my $self = shift;
    my $instrument_data = shift;
    my $analysis_project = shift;

    my $bridge = Genome::Config::AnalysisProject::InstrumentDataBridge->get(
        instrument_data => $instrument_data,
        analysis_project => $analysis_project,
    );
    my $config = $analysis_project->get_configuration_profile();
    my @maps = $config->config_rule_maps();

    print sprintf("Comparing Instrument Data %s to Analysis Project %s:\n",
        $self->_color($instrument_data->id,'cyan'),
        $self->_color($analysis_project->name,'magenta')
    );

    my $cqid_status = $self->_colorize_text_by_map($bridge->status, $bridge->status, _bridge_status_colors());
    printf("  CQID Status: %s\n", $cqid_status);

    for my $map (@maps) {
        my @rules = $map->rules;

        my @mismatches;
        for my $rule (@rules) {
            my ($match, $actual) = $map->evaluate_rule($rule,$instrument_data);
            unless($match) {
                push @mismatches, [$rule, $actual];
            }
        }

        if(@mismatches) {
            print sprintf("  Mismatches for config %s:\n",
                $self->_color($map->config->file_path,'red')
            );
            for my $m (@mismatches) {
                print sprintf("    got %s for rule %s\n",
                    $self->_color($m->[1],'cyan'),
                    $self->_color($m->[0]->__display_name__,'magenta')
                );
            }

        } else {
            if ($map->config->has_model_for($instrument_data)) {
                print sprintf("  Matches config %s.\n",
                    $self->_color($map->config->file_path,'green')
                );
            } else {
                printf("  Missing expected model for matched config %s.\n",
                    $self->_color($map->config->file_path, 'blue'));
            }
        }
    }
}

1;
