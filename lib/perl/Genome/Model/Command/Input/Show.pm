package Genome::Model::Command::Input::Show;

use strict;
use warnings;
use List::Util qw(max first);
require Term::ANSIColor;

use Genome;
use Genome::Utility::Text "justify";


class Genome::Model::Command::Input::Show {
    is => ['Genome::Command::Viewer', 'Genome::Command::WithColor'],
    has => [
        model => {
            is => 'Genome::Model',
            shell_args_position => 1,
        },
        show_display_names => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
            doc => "Show display_name instead ID for each input."
        },
    ],
    doc => 'Show the inputs of a model.',
};

sub write_report {
    my ($self, $width, $handle) = @_;
    my $model = $self->model;

    $self->write_inputs_for_model_or_build(
            'width' => $width,
            'handle' => $handle,
            'target' => $model,
            'target_type' => "model",
            'show_display_names' => $self->show_display_names,
            'color' => $self->color,
    );
}

sub _get_sorted_input_properties {
    my ($self, $target) = @_;

    if($target->can('real_input_properties')) {
        return $target->real_input_properties;
    } else {
        die $self->error_message('Could not load properties for target.'); #TODO support builds?
    }

}

sub write_inputs_for_model_or_build {
    my $self = shift;
    my %params = @_;
    my $width = $params{width};
    my $handle = $params{handle};
    my $target = $params{target};
    my $target_type = $params{target_type};
    my $show_display_names = $params{show_display_names};
    my $color = $params{color};

    my @properties = $self->_get_sorted_input_properties($target);

    unless(@properties) {
        printf $handle "No inputs found for %s %s",
                $target_type, $target->id;
        return;
    }

    my $name_header = 'Input Name';
    my $is_many_header = "Is Many";
    my $value_header = 'Value(s)';

    # get contents of table
    my %inputs;
    my %is_many;
    my @input_name_lengths = (length($name_header));
    my @input_value_lengths = (length($value_header));
    for my $property (@properties) {
        my $name = $property->{name};
        push(@input_name_lengths, length($name));

        my @values;
        for my $value ($target->$name) {
            if($value) {
                if($show_display_names and 
                        eval {$value->can('__display_name__')}) {
                    push(@values, $value->__display_name__);
                } elsif(eval {$value->can('class')} and 
                        eval {$value->can('id')}) {
                    push(@values, sprintf("%s(%s)", $value->class, $value->id));
                } else {
                    push(@values, $value);
                }
            } else {
                push(@values, 'undef');
            }
        }
        unless(@values) {
            @values = ('undef');
        }
        my @value_lengths = map { length($_) } @values;
        @input_value_lengths = (@input_value_lengths, @value_lengths);

        $inputs{$name} = \@values;
        $is_many{$name} = $property->{is_many};
    }

    # determine widths of columns of table
    my $max_name_length = max(@input_name_lengths);
    my $max_value_length = max(@input_value_lengths);

    # print out table header
    printf $handle "%s %s %s\n", 
            justify($name_header, 'right', $max_name_length),
            $is_many_header,
            justify($value_header, 'left', $max_value_length);

    printf $handle "%s %s %s\n", '-' x $max_name_length,
                                 '-' x length($is_many_header),
                                 '-' x $max_value_length;

    # print out table
    for my $name (sort keys %inputs) {
        my $name_part = justify($name, 'right', $max_name_length, " ", "");
        my $is_many_part = sprintf('  %s  ',
                $self->_format_is_many($is_many{$name}), $color);
        my $value_part = $self->_format_values($inputs{$name},
                $max_name_length + length($is_many_header) + 2,
                $max_value_length,
                $color);
        print $handle join(' ', $name_part, $is_many_part, $value_part);
    }
    print $handle "\n";
}


sub _format_is_many {
    my ($self, $is_many, $color) = @_;
    my ($pre, $post, $true, $false) = ('[', ']', 'X',' ');

    $pre = $self->_color($pre, 'white', $color);
    $post = $self->_color($post, 'white', $color);
    my $mid = $is_many ? $true : $false;
    return join('', $pre, $mid, $post);
}

sub _color {
    my ($self, $value, $color, $flag) = @_;
    if($flag) {
        return Term::ANSIColor::colored($value, $color);
    } else {
        return $value;
    }
}

sub _format_values {
    my ($self, $values, $left_padding, $right_size, $color) = @_;
    my @values = @{$values};

    my @formatted_values;
    my $first_value = 1;
    for my $value (@values) {
        my $padded_value;
        if($first_value) {
            $padded_value = justify($value, 'left', $right_size);
            $first_value = 0;
        } else {
            $padded_value = justify($value, 'left', $right_size);
            $padded_value = justify($padded_value, 'right',
                    $right_size + $left_padding);
        }
        if($value eq 'undef') {
            $padded_value = $self->_color($padded_value, 'cyan', $color);
        }
        push(@formatted_values, $padded_value);
    }

    my $result = join("\n", @formatted_values);
    $result .= "\n";
    return $result;
}

1;


