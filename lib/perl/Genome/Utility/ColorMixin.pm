package Genome::Utility::ColorMixin;

use strict;
use warnings;

use Genome;
use Genome::Utility::Text qw(justify);

class Genome::Utility::ColorMixin {
    has => [
        color => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
            doc => 'Use color in display.'
        },
    ],
};

sub _color {
    my $self = shift;
    my $string = shift;

    if(-t STDOUT and -t STDERR and $self->color and @_) {
        return Term::ANSIColor::colored($string, @_);
    } else {
        return $string;
    }
}

our %STATUS_COLORS = (
    new => "white",
    scheduled => "white",

    running => "cyan",
    'running*' => 'cyan',

    done      => "green",
    succeeded => "green",

    abandoned => "magenta",

    crashed     => "red",
    failed      => "red",
    unstartable => "red",
);

sub _status_colors {
    my $self = shift;
    my $status = shift;
    return $STATUS_COLORS{$status};
}

sub _status_color {
    my ($self, $text) = @_;
    return $self->_colorize_text_by_map($text, $text, %STATUS_COLORS);
}

sub _colorize_text_by_map {
    my ($self, $text, $color_key, %color_map) = @_;

    my $stripped_key = $self->_strip_key($color_key);
    if (exists $color_map{$stripped_key}) {
        return $self->_color($text, $color_map{$stripped_key});
    }

    return $text;
}

sub _strip_key {
    my ($self, $text) = @_;

    my $stripped_text = $text;
    $stripped_text =~ tr/A-Z/a-z/;
    $stripped_text =~ s/ //g;

    return $stripped_text;
}

sub _color_heading {
    my ($self, $text) = @_;
    return $self->_color_dim('=== ') . $self->_color($text, 'bold') .
        $self->_color_dim(' ===');
}

sub _color_pair {
    my ($self, $key, $value) = @_;
    return $self->_color_dim($key.':') . ' ' . ($value || '');
}

sub _color_dim {
    my ($self, $text) = @_;
    return $self->_color($text, 'white');
}

sub _column_width {
    return Genome::Command::Viewer->get_terminal_width / 2;
}

sub _write_pairs_line {
    my ($self, $handle, $l_label, $l_value, $r_label, $r_value) = @_;

    if ($r_label and $r_value) {
        print $handle justify($self->_color_pair($l_label, $l_value), 'left',
            $self->_column_width), " ", $self->_color_pair($r_label, $r_value), "\n";

    } else {
        print $handle $self->_color_pair($l_label, $l_value), "\n";
    }
}

1;
