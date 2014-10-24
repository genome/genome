package Genome::File::BedPe::Entry;

use strict;
use warnings;

use Genome;
use Carp qw/confess/;
use List::AllUtils qw/min first/;

use Genome::File::BedPe::Header qw(
    @REQUIRED_FIELDS @EXTRA_FIELDS @ALL_FIELDS
    %FIELD_INDICES
    );

sub new {
    my ($class, $header, $line) = @_;
    my $self = {
        header => $header,
        custom => [],
    };

    bless $self, $class;
    $self->_parse($line);
    return $self;
}

sub new_from_fields {
    my ($class, $header, %fields) = @_;


    my %defaults = map {$_ => '.'} keys %FIELD_INDICES;
    my $self = {
        header => $header,
        custom => [],
        %defaults,
    };


    my %standard_mapping;
    my %custom_mapping;
    my %required_seen;
    @required_seen{@REQUIRED_FIELDS} = undef;

    my $max_custom = -1;
    for my $f (keys %fields) {
        if (exists $FIELD_INDICES{$f}) {
            delete $required_seen{$f};
            $self->{$f} = $fields{$f};
        }
        else {
            my $cidx = $header->custom_field_index($f);
            if (!defined $cidx) {
                confess "Unknown BedPe attribute $f specified. Custom fields can be added " .
                        "to the header with \$header->set_custom_fields(...)";
            }
            $custom_mapping{$f} = $cidx;
            $max_custom = $cidx > $max_custom ? $cidx : $max_custom;
        }
    }

    my @missing = keys %required_seen;
    confess sprintf("Required fields missing: %s", join(", ", @missing)) if @missing;

    $self->{custom} = [('.') x $max_custom];

    for my $f (keys %custom_mapping) {
        $self->{custom}->[$custom_mapping{$f}] = $fields{$f};
    }

    bless $self, $class;

    return $self;
}

sub _parse {
    my ($self, $line) = @_;
    my @fields = split("\t", $line);
    if (scalar @fields < scalar @REQUIRED_FIELDS) {
        confess "Too few fields in bedpe record (" . scalar @fields . ")";
    }

    my $last_field = min($#fields, $#ALL_FIELDS);

    @{$self}{@ALL_FIELDS[0..$last_field]} = @fields[0..$last_field];
    if ($#fields > $last_field) {
        $self->{custom} = [ @fields[$last_field + 1 .. $#fields] ];
    }
}

sub custom_by_name {
    my ($self, $name) = @_;
    my $index = $self->{header}->custom_field_index($name);
    if (!defined $index) {
        print Data::Dumper::Dumper($index);
        confess "Unknown custom value $name";
    }
    return $self->{custom}->[$index];
}

sub validate {
    my $self = shift;
    # numeric fields, (-1 means not known)
    my @numeric_fields = qw/start1 start2 end1 end2/;

    for my $f (@numeric_fields) {
        next if $f eq '-1';
        if ($self->{$f} !~ /^\d+$/) {
            my $line = $self->to_string;
            confess "in entry $line: $f is not numeric";
        }
    }
}

sub to_string {
    my $self = shift;
    my $last_field = first { !defined $self->{$ALL_FIELDS[$_]} } 0..$#ALL_FIELDS;
    if (!defined $last_field) {
        $last_field = $#ALL_FIELDS;
    }
    else {
        --$last_field;
    }

    return join("\t", @{$self}{@ALL_FIELDS[0..$last_field]}, @{$self->{custom}});
}

1;
