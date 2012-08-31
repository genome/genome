package Genome::Model::Command::Input::Show;

use strict;
use warnings;
use List::Util 'max';
require Term::ANSIColor;

use Genome;


class Genome::Model::Command::Input::Show {
    is => ['Command::V2'],
    has => [
        model => {
            is => 'Genome::Model',
            shell_args_position => 1,
        },
        color => {
            is => 'Boolean',
            default_value => 1,
            doc => 'Display in color.'
        },
    ],
    doc => 'Show the inputs of a model or type.',
};

sub execute {
    my ($self) = @_;
    my $model = $self->model;

    my @properties = eval{ $model->real_input_properties };
    return if not @properties;

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
        for my $value ($model->$name) {
            if($value) {
                if ( eval{ $value->can('__display_name__') } ) {
                    push(@values, $value->__display_name__);
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
    my $header1 = sprintf("%s %s %s\n", _pad_left($name_header, $max_name_length),
                                       $is_many_header,
                                       _pad_right($value_header, $max_value_length));
    $self->status_message($header1);
    my $header2 = sprintf("%s %s %s\n", '-' x $max_name_length,
                                       '-' x length($is_many_header),
                                       '-' x $max_value_length);
    $self->status_message($header2);

    # print out table
    for my $name (sort keys %inputs) {
        my $name_part = _pad_left($name, $max_name_length);
        my $is_many_part = '  ' . _format_is_many($is_many{$name}, $self->color) . '  ';
        my $value_part = _format_values($inputs{$name},
                                         $self->color,
                                         $max_name_length + length($is_many_header) + 2,
                                         $max_value_length);
        $self->status_message(join(' ', $name_part, $is_many_part, $value_part));
    }

    return 1;
}


sub _format_is_many {
    my ($is_many, $color) = @_;
    my ($pre, $post, $true, $false) = ('[', ']', '*',' ');

    if($color) {
        $pre = Term::ANSIColor::colored($pre, 'white');
        $post = Term::ANSIColor::colored($post, 'white');
    }
    my $mid = $is_many ? $true : $false;
    return join('', $pre, $mid, $post);
}

sub _pad_left {
    my ($arg, $length) = @_;
    return sprintf("% $length"."s", $arg);
}

sub _pad_right {
    my ($arg, $length) = @_;
    return sprintf("%-$length"."s", $arg);
}

sub _format_values {
    my ($values, $color, $left_padding, $right_size) = @_;
    my @values = @{$values};

    my @formatted_values;
    my $first_value = 1;
    for my $value (@values) {
        my $padded_value;
        if($first_value) {
            $padded_value = _pad_right($value, $right_size);
            $first_value = 0;
        } else {
            $padded_value = _pad_right($value, $right_size);
            $padded_value = _pad_left($padded_value, $right_size + $left_padding);
        }
        if($color and $value eq 'undef') {
            $padded_value = Term::ANSIColor::colored($padded_value, 'cyan');
        }
        push(@formatted_values, $padded_value);
    }

    my $result = join("\n", @formatted_values);
    $result .= "\n";
    return $result;
}

1;


