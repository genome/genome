package Genome::Config::AnalysisProject::Command::ConfigForInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::ConfigForInstrumentData {
    is => ['Command::V2'],
    roles => ['Genome::Role::CommandWithColor'],
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
    has_optional => [
        config_file => {
            is => 'Path',
            doc => 'a specific config file to evaluate instead of existing AnP config',
        },
    ],
    doc => 'compare instdata to AnP configs to show matches, models and mismatches',
};

sub execute {
    my $self = shift;
    my @instrument_data = $self->instrument_data;

    for my $id (@instrument_data) {
        if (my $config_file = $self->config_file) {
            my $parsed_file = Genome::Config::Parser->parse($config_file);
            my @rules = Genome::Config::Rule->create_from_hash($parsed_file->{rules});

            print sprintf("Comparing Instrument Data %s to supplied config file:\n",
                $self->_color($id->id,'cyan'),
            );

            $self->_display_match($id, $config_file, 0, @rules);
        } else {
            my @analysis_projects = $id->analysis_projects;
            unless(@analysis_projects) {
                $self->warning_message('No analysis-projects associated with instrument data %s.', $id->__display_name__);
                return 1;
            }

            for my $ap (@analysis_projects) {
                $self->_display_matches($id, $ap);
            }
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
        $self->_display_match($instrument_data, $map->config->file_path, !$map->config->has_model_for($instrument_data), $map->rules);
    }
}

sub _display_match {
    my $self = shift;
    my $instrument_data = shift;
    my $file_path = shift;
    my $model_missing = shift;
    my @rules = @_;

    my @mismatches;
    for my $rule (@rules) {
        my ($match, $actual) = Genome::Config::RuleModelMap->evaluate_rule($rule,$instrument_data);
        unless($match) {
            push @mismatches, [$rule, $actual];
        }
    }

    if(@mismatches) {
        print sprintf("  Mismatches for config %s:\n",
            $self->_color($file_path,'red')
        );
        for my $m (@mismatches) {
            print sprintf("    got %s for rule %s\n",
                $self->_color($m->[1],'cyan'),
                $self->_color($m->[0]->__display_name__,'magenta')
            );
        }

    } else {
        unless ($model_missing) {
            print sprintf("  Matches config %s.\n",
                $self->_color($file_path,'green')
            );
        } else {
            printf("  Missing expected model for matched config %s.\n",
                $self->_color($file_path, 'blue'));
        }
    }
}

1;
