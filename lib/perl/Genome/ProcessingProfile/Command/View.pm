package Genome::ProcessingProfile::Command::View;

use strict;
use warnings;
use feature 'switch';

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Genome;
use Term::ReadKey 'GetTerminalSize';
use List::MoreUtils "uniq";

class Genome::ProcessingProfile::Command::View {
    is => 'Command::V2',
    has => [
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            shell_args_position => 1,
            doc => 'Genome::ProcessingProfile',
        },
        color => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
            doc => 'Display report in color.'
        },
        max_num_models_shown => {
            is => 'Int',
            is_optional => 1,
            default_value => 10,
            doc => 'Display no more than this many models (-1 for all).'
        },
        show_strategies_as => {
            is => 'Text',
            is_optional => 1,
            default => 'condensed-tree',
            valid_values => ['raw', 'raw-tree', 'condensed-tree'],
            doc => 'How the strategy parameters (if there are any) will be displayed',
        }
    ],
};

sub help_synopsis {
    my $self = shift;
    my $result .= <<EOP;
    Displays basic information about a processing-profile.
EOP
    return $result;
}

sub help_detail {
    my $self = shift;
    my $result .= <<EOP;
Displays information about the parameters for a processing-profile as well as
what models are using it.
EOP
    return $result;
}

sub execute {
    my ($self) = @_;

    my ($screen_width) = GetTerminalSize();
    print $self->_get_basic_info($screen_width);
    print $self->_get_parameters($screen_width);
    print $self->_get_strategies($screen_width);
    print $self->_get_models($screen_width);

    1;
}

sub get_report {
    my ($self, $width) = @_;

    my $report = $self->_get_basic_info($width);
    $report .= $self->_get_parameters($width);
    $report .= $self->_get_strategies($width);
    $report .= $self->_get_models($width);

    return $report;
}

sub _get_basic_info {
    my ($self, $width) = @_;
    my $pp = $self->processing_profile;

    my $result = '';
    $result .= sprintf("\nName: %s\n", $self->_color($pp->name, 'bold'));
    $result .= sprintf("ID: %s   ", $self->_color($pp->id, 'bold'));
    $result .= sprintf("Type Name: %s\n", $self->_color($pp->type_name, 'bold'));

    return $result;
}

sub _get_parameters {
    my ($self, $width) = @_;

    my $result = '';

    my @params = $self->processing_profile->params;
    @params = grep(not ($_->name =~ m/_strategy$/), @params);
    my @columns;
    my @justifications;
    my @fills;

    my $name_column = join("\n", map {$_->name} @params);
    push(@columns, $name_column);
    push(@justifications, "left");
    push(@fills, ".");

    my @value_class_names = map {$_->value_class_name} @params;
    if(scalar(uniq(@value_class_names) > 1)) {
        push(@columns, join("\n", @value_class_names));
        push(@justifications, "left");
        push(@fills, " ");
    }

    my $value_id_column = join("\n",
            map {$self->_color_first_line($_->value_id, 'bold')} @params);
    push(@columns, $value_id_column);
    push(@justifications, "left");
    push(@fills, " ");

    $result .= Genome::Utility::Text::side_by_side(\@columns,
            separator => ' ',
            justification => \@justifications,
            fill => \@fills,
            max_width => $width,
    );
    $result .= "\n";
    return $result;
}

sub _get_strategies {
    my ($self, $width) = @_;

    my $result = "\n";
    my @fss;
    my @params = $self->processing_profile->params;
    @params = grep($_->name =~ m/_strategy$/, @params);
    for my $param (@params) {
        my $name = $param->name;
        my $strategy_str = $param->value_id;

        my $formatted_strategy_str;
        given($self->show_strategies_as) {
            when('raw') {
                $formatted_strategy_str = sprintf("%s: %s\n", $name,
                        $self->color($strategy_str, 'bold'));
            }
            when('raw-tree') {
                my $tree = _get_strategy_tree($strategy_str);
                $formatted_strategy_str =
                        Genome::Utility::Text::tree_to_string($tree);
                $formatted_strategy_str = "$name\n" . $formatted_strategy_str;
            }
            when('condensed-tree') {
                my $tree = _get_strategy_tree($strategy_str);
                $formatted_strategy_str =
                        Genome::Utility::Text::tree_to_condensed_string(
                        _format_tree($tree));
                $formatted_strategy_str = "$name\n" . $formatted_strategy_str;
            }
        }
        $formatted_strategy_str = Genome::Utility::Text::strip_color(
            $formatted_strategy_str) unless $self->color;
        push(@fss, $formatted_strategy_str);
    }

    $result .= Genome::Utility::Text::side_by_side(\@fss,
            separator => ' | ',
            justification => 'left',
            max_width => $width,
            stack => 1,
    );
    $result .= "\n";
    return $result;
}

sub _get_models {
    my ($self, $width) = @_;

    my $result = '';

    my @models = $self->processing_profile->models;
    my $num_models = scalar(@models);

    my $models_shown = 0;
    my @model_names;
    my @model_ids;
    my @subject_ids;
    my @subject_names;
    for my $model (@models) {
        if($models_shown >= $self->max_num_models_shown and
           $self->max_num_models_shown != -1) {
            last;
        }

        my $color = 'bold';
        if($self->color) {
            my $latest_build = $model->latest_build;
            if(defined($latest_build)) {
                my $status = $latest_build->status;
                if(defined($status)) {
                    if($status eq 'Succeeded') {
                        $color = 'green bold';
                    } else {
                        $color = 'red bold';
                    }
                }
            }
        }
        push(@model_names,   $self->_color($model->name, 'bold'));
        push(@model_ids,     $model->id);
        if(defined($model->subject)){
            push(@subject_ids,   $model->subject->id);
            push(@subject_names, $self->_color($model->subject->name, 'bold'));
        } else {
            push(@subject_ids, $self->_color("undef", 'bold red'));
            push(@subject_names, $self->_color('undef', 'bold red'));
        }

        $models_shown += 1;
    }

    my $model_name_column = join("\n", "Model Name",
            map {$self->_color($_, 'bold green')} @model_names);
    my $model_id_column = join("\n", "Model ID", @model_ids);
    my $subject_id_column = join("\n", "Sample ID", @subject_ids);
    my $subject_name_column = join("\n", "Sample Name",
            map {$self->_color($_ . '', 'bold')} @subject_names);

    my $models = Genome::Utility::Text::side_by_side(
            [$model_id_column, $model_name_column],
            separator => ' ',
            justification => 'left',
            max_width => $width,
    );
    my $subjects = Genome::Utility::Text::side_by_side(
            [$subject_id_column, $subject_name_column],
            separator => ' ',
            justification => 'left',
            max_width => $width,
    );

    $result .= sprintf("\n=== %s [%s of %s shown] ===\n",
        $self->_color('Models', 'bold'), $models_shown, $num_models);

    $result .= Genome::Utility::Text::side_by_side(
            [$models, $subjects],
            separator => ' ',
            justification => 'left',
            max_width => $width,
    );
    $result .= "\n";

    if($models_shown >= $self->max_num_models_shown and
       $self->max_num_models_shown != -1) {
        $result .= sprintf("... %s of %s models shown" .
                " (see max-num-models-shown option)\n",
                $models_shown, $num_models);
    }

    return $result;
}

sub _color {
    my $self = shift;
    my $string = shift;

    if($self->color) {
        return Term::ANSIColor::colored($string, @_);
    } else {
        return $string;
    }
}

sub _color_first_line {
    my $self = shift;
    my $string = shift;

    my @lines = split("\n", $string);
    return join("\n", $self->_color($lines[0], 'bold'), @lines[1..$#lines]);
}

sub _get_strategy_tree {
    my ($strategy_str) = @_;

    my $strategy = Genome::Model::Tools::DetectVariants2::Strategy->create($strategy_str);
    my $tree = $strategy->parse($strategy_str);
    $strategy->delete();
    return $tree;
}

sub _format_tree {
    my ($tree) = @_;

    my @result;
    my @keys = keys %{$tree};
    for my $key (sort(@keys)) {
        if($key eq 'detector') {
            my $detector_info = $tree->{$key};
            delete $tree->{detector};
            my ($new_key, $new_value) = _format_detector($detector_info);
            $tree->{$new_key} = $new_value;
        } else {
            my @steps = @{$tree->{$key}};
            my @new_value;
            for my $step (@steps) {
                push(@new_value, _format_tree($step));
            }
            delete $tree->{$key};
            $key = Term::ANSIColor::colored($key, 'bold magenta');
            $tree->{$key} = \@new_value;
        }
    }
    return $tree;
}

sub _format_detector {
    my ($detector_info) = @_;

    # a detector is almost exactly like a filter but with filters...
    my ($new_key, $params) = _format_filter($detector_info, 'green');
    my @new_value = @{$params};
    my @filters = @{$detector_info->{filters}};
    for my $filter_info (@filters) {
        my ($filter_name, $filter_value) = _format_filter($filter_info, 'red');
        push(@new_value, {"filtered by $filter_name" => $filter_value});
    }
    return $new_key, \@new_value;
}

sub _format_filter {
    my ($filter_info, $color) = @_;

    my $name = $filter_info->{name};
    my $version = $filter_info->{version};
    my $param_str = $filter_info->{params};
    my $new_key = Term::ANSIColor::colored($name, "bold $color") .
                  Term::ANSIColor::colored(" $version", 'bold');
    my $new_value = [];
    if($param_str) {
        $new_value = _format_params($param_str);
    }
    return $new_key, $new_value;
}

sub _format_params {
    my ($param_str) = @_;

    my %params = Genome::Utility::Text::param_string_to_hash($param_str);
    my @result;
    for my $key (sort(keys %params)) {
        my $value = Term::ANSIColor::colored($params{$key}, 'bold');
        push(@result, sprintf("%s: %s", $key, $value));
    }
    return \@result;
}

1;
