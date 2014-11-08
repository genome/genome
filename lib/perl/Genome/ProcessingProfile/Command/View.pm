package Genome::ProcessingProfile::Command::View;

use strict;
use warnings;
use feature 'switch';

use Genome;
use List::MoreUtils "uniq";
use Genome::Utility::Text qw(strip_color
                             param_string_to_hash
                             tree_to_string
                             tree_to_condensed_string
                             side_by_side);

class Genome::ProcessingProfile::Command::View {
    doc => "Display basic information about a processing-profile.",
    is => ['Genome::Command::Viewer', 'Genome::Command::WithColor'],
    has => [
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            shell_args_position => 1,
            doc => 'Genome::ProcessingProfile',
        },
        models => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'Display which models are using this processing-profile.',
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
        },
        latest_build_status => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'Display latest build status for each model.'
        },
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

sub write_report {
    my ($self, $width, $handle) = @_;

    $self->_write_basic_info($width, $handle);
    $self->_write_parameters($width, $handle);
    $self->_write_strategies($width, $handle);
    $self->_write_models($width, $handle);
}

sub _write_basic_info {
    my ($self, $width, $handle) = @_;
    my $pp = $self->processing_profile;

    printf $handle "\nName: %s\n", $self->_color($pp->name, 'bold');
    printf $handle "ID: %s   ", $self->_color($pp->id, 'bold');
    printf $handle "Type Name: %s\n", $self->_color($pp->type_name, 'bold');
}

sub _write_parameters {
    my ($self, $width, $handle) = @_;

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
            map {$self->_color($_->value_id, 'bold')} @params);
    push(@columns, $value_id_column);
    push(@justifications, "left");
    push(@fills, " ");

    printf $handle "%s\n", side_by_side(\@columns,
            separator => ' ',
            justification => \@justifications,
            fill => \@fills,
            max_width => $width,
    );
}

sub _write_strategies {
    my ($self, $width, $handle) = @_;

    my @params = $self->processing_profile->params;
    @params = grep($_->name =~ m/_strategy$/, @params);

    return unless scalar(@params);

    print $handle "\n";
    my @fss;
    for my $param (@params) {
        my $name = $param->name;
        my $strategy_str = $param->value_id;

        my $formatted_strategy_str;
        given($self->show_strategies_as) {
            when('raw') {
                $formatted_strategy_str = sprintf("%s: %s\n", $name,
                        $self->_color($strategy_str, 'bold'));
            }
            when('raw-tree') {
                $formatted_strategy_str = $self->_format_tree(
                        $name, $strategy_str, 1);

            }
            when('condensed-tree') {
                $formatted_strategy_str = $self->_format_tree(
                        $name, $strategy_str, 0);
            }
        }
        $formatted_strategy_str = strip_color(
            $formatted_strategy_str) unless $self->color;
        push(@fss, $formatted_strategy_str);
    }

    printf $handle "%s\n", side_by_side(\@fss,
            separator => ' | ',
            justification => 'left',
            max_width => $width,
            stack => 1,
    );
}

our %STRATEGY_PARSERS = (
    alignment_strategy => \&_get_alignment_strategy_tree,
    snv_detection_strategy => \&_get_dv2_strategy_tree,
    sv_detection_strategy => \&_get_dv2_strategy_tree,
    indel_detection_strategy => \&_get_dv2_strategy_tree,
    cnv_detection_strategy => \&_get_dv2_strategy_tree,
);

our %STRATEGY_FORMATTERS = (
    alignment_strategy => \&_format_alignment_strategy_tree,
    snv_detection_strategy => \&_format_dv2_strategy_tree,
    sv_detection_strategy => \&_format_dv2_strategy_tree,
    indel_detection_strategy => \&_format_dv2_strategy_tree,
    cnv_detection_strategy => \&_format_dv2_strategy_tree,
);

sub _format_tree {
    my ($self, $name, $strategy_str, $raw) = @_;

    my $parser = $STRATEGY_PARSERS{$name};
    unless(defined($parser)) {
        my $warning = "I don't know how to parse this strategy type...";
        return sprintf("%s\n  %s\n  %s", $name,
                    $self->_color($warning, 'red'), $strategy_str);
    }

    my $tree = $parser->($strategy_str);
    my $result;
    if($tree) {
        if($raw) {
            $result = tree_to_string($tree);
        } else {
            my $formatter = $STRATEGY_FORMATTERS{$name};
            my $formatted_tree = $formatter->($tree, $strategy_str);
            $result = tree_to_condensed_string($formatted_tree);
        }
    } elsif($strategy_str) {
        $result = $self->_color("Error parsing strategy!", 'red');
    }
    $result = "$name\n" . $result;

    return $result;
}

sub _get_alignment_strategy_tree {
    my ($strategy_str) = @_;

    my $strategy = Genome::InstrumentData::Composite::Strategy->create(strategy=>$strategy_str);
    my $tree = $strategy->execute();
    $strategy->delete();
    return $tree;
}

sub _format_alignment_strategy_tree {
    my ($tree, $strategy_str) = @_;

    my $actions = [];
    _mine_tree_for_type($tree, 'type', $actions);

    my $names = [];
    _mine_tree_for_type($tree, 'name', $names);

    my %color_dict;
    for my $key (@{$actions}) {
        $color_dict{$key . 'd'} = 'bold magenta';
        $color_dict{$key . 'ed'} = 'bold magenta';
    }
    for my $key (@{$names}) {
        $color_dict{$key} = 'bold red';
    }

    my @words = split(/\s/, $strategy_str);
    my @formatted_words;
    for my $word (@words) {
        if($word eq 'then') {
            push(@formatted_words, sprintf("\n    %s", $word));
        } elsif($color_dict{$word}) {
            my $color = $color_dict{$word};
            push(@formatted_words, Term::ANSIColor::colored($word, $color));
        } else {
            push(@formatted_words, $word);
        }
    }
    return join(' ', @formatted_words);
}

sub _mine_tree_for_type {
    my ($tree, $name, $bucket) = @_;

    if(ref($tree) eq 'HASH') {
        for my $key (keys %{$tree}) {
            if($key eq $name) {
                push(@{$bucket}, $tree->{$key});
            } else {
                _mine_tree_for_type($tree->{$key}, $name, $bucket);
            }
        }
    } elsif(ref($tree) eq 'ARRAY') {
        for my $item (@{$tree}) {
            _mine_tree_for_type($item, $name, $bucket);
        }
    }
    return;
}

sub _get_dv2_strategy_tree {
    my ($strategy_str) = @_;

    my $strategy = Genome::Model::Tools::DetectVariants2::Strategy->create($strategy_str);
    my $tree = $strategy->parse($strategy_str);
    $strategy->delete();
    return $tree;
}

sub _format_dv2_strategy_tree {
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
                push(@new_value, _format_dv2_strategy_tree($step));
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

    my $name = $detector_info->{name};
    my $version = $detector_info->{version};
    my $param_str = $detector_info->{params};
    my $new_key = Term::ANSIColor::colored($name, "bold cyan") .
                  Term::ANSIColor::colored(" $version", 'bold');
    my @new_value;
    if($new_key =~ m/strelka/) { #TODO make this less specific
        # show params for filters differently than for detectors... apparently
        # in DV2 the detectors don't have a particular format for their params,
        # it is just detector specific, whereas filters are all command style
        # params.
        use Genome::Model::Tools::DetectVariants2::Strelka;
        my %params = Genome::Model::Tools::DetectVariants2::Strelka::parse_params($param_str);
        push(@new_value, _format_params(\%params));
    } else {
        if($param_str) {
            my $display_str = sprintf("params: %s",
                    Term::ANSIColor::colored($param_str, 'bold'));
            push(@new_value, $display_str);
        }
    }
    my @filters = @{$detector_info->{filters}};
    for my $filter_info (@filters) {
        my ($filter_name, $filter_value) = _format_filter($filter_info);
        push(@new_value, {"filtered by $filter_name" => $filter_value});
    }
    return $new_key, \@new_value;
}

sub _format_filter {
    my ($filter_info) = @_;

    my $name = $filter_info->{name};
    my $version = $filter_info->{version};
    my $param_str = $filter_info->{params};
    my $new_key = Term::ANSIColor::colored($name, "bold red") .
                  Term::ANSIColor::colored(" $version", 'bold');
    my $new_value = [];
    if($param_str) {
        $new_value = _format_params($param_str);
    }
    return $new_key, $new_value;
}

sub _format_params {
    my ($params) = @_;

    my %params;
    if(ref($params) eq 'HASH') {
        %params = %{$params};
    } else {
        %params = param_string_to_hash($params);
    }

    my @result;
    for my $key (sort(keys %params)) {
        my $value = Term::ANSIColor::colored($params{$key}, 'bold');
        push(@result, sprintf("%s: %s", $key, $value));
    }
    return \@result;
}

sub _write_models {
    my ($self, $width, $handle) = @_;

    return unless $self->models and $self->max_num_models_shown != 0;
    my @models = $self->processing_profile->models;
    my $num_models = scalar(@models);

    my $models_shown = 0;
    my @model_names;
    my @model_ids;
    my @latest_build_ids;
    my @latest_build_statuses;
    my @subject_ids;
    my @subject_names;
    for my $model (@models) {
        if($models_shown >= $self->max_num_models_shown and
           $self->max_num_models_shown != -1) {
            last;
        }

        if($self->latest_build_status) {
            my $color = 'bold';
            my $latest_build = $model->latest_build;
            if(defined($latest_build)) {
                my $status = $latest_build->status;
                if(defined($status)) {
                    if($status eq 'Succeeded') {
                        $color = 'cyan bold';
                    } else {
                        $color = 'red bold';
                    }
                }
                push(@latest_build_statuses, $self->_color($status, $color));
                push(@latest_build_ids, $latest_build->id);
            } else {
                push(@latest_build_statuses, $self->_color("N/A", 'cyan'));
                push(@latest_build_ids, $self->_color("N/A", 'cyan'));
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

    my $model_name_column = join("\n",
            $self->_color("Model Name", 'bold'), @model_names);
    my $model_id_column = join("\n", "Model ID", @model_ids);
    my $models = side_by_side(
            [$model_id_column, $model_name_column],
            separator => ' ',
            justification => 'left',
            max_width => $width,
    );
    my @all_columns = ($models);

    my $subject_id_column = join("\n", "Subject ID", @subject_ids);
    my $subject_name_column = join("\n",
            $self->_color("Subject Name", 'bold'), @subject_names);
    my $subjects = side_by_side(
            [$subject_id_column, $subject_name_column],
            separator => ' ',
            justification => 'left',
            max_width => $width,
    );
    push(@all_columns, $subjects);

    if($self->latest_build_status) {
        my $latest_build_id_column = join("\n", "Latest Build ID",
                @latest_build_ids);
        my $latest_build_status_column = join("\n",
                $self->_color("Status", 'bold'),
                @latest_build_statuses);
        my $builds = side_by_side(
                [$latest_build_id_column, $latest_build_status_column],
                separator => ' ',
                justification => 'left',
                max_width => $width,
        );
        push(@all_columns, $builds);
    }

    printf $handle "\n=== %s [%s of %s shown] ===\n",
            $self->_color('Models', 'bold'), $models_shown, $num_models;
    printf $handle "%s\n", side_by_side(
            \@all_columns,
            separator => ' ',
            justification => 'left',
            max_width => $width,
    );

    if($models_shown >= $self->max_num_models_shown and
       $self->max_num_models_shown != -1) {
        printf $handle "... %s of %s models shown" .
                " (see max-num-models-shown option)\n",
                $models_shown, $num_models;
    }
}


1;
