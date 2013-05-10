package Genome::Model::Tools::Sx::Functions;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::Sx::Functions {};

#< Quality Calculation >#
sub calculate_quality {
    my ($self, $quality_string) = @_;
    Carp::confess('No quality string to calculate average.') if not defined $quality_string and not length $quality_string;
    my $total = 0;
    for my $q ( split('', $quality_string) ) {
        $total += ord($q) - 33;
    }
    return $total;
}

sub calculate_average_quality {
    my ($self, $quality_string) = @_;
    Carp::confess('No quality string to calculate average.') if not defined $quality_string and not length $quality_string;
    return sprintf('%.0f', __PACKAGE__->calculate_quality($quality_string) / length($quality_string));
}

sub calculate_qualities_over_minumum {
    my ($self, $quality_string, $min) = @_;
    Carp::confess('No quality string to calculate qualities over minimum.') if not defined $quality_string and not length $quality_string;
    Carp::confess('No minimum calculate qualities over minimum.') if not defined $min;
    my $total = 0;
    for my $q ( split('', $quality_string) ) {
        next if (ord($q) - 33) < $min;
        $total++;
    }
    return $total;
}

sub maximum_quality {
    my ($self, $quality_string) = @_;
    Carp::confess('No quality string to get maximum quality.') if not defined $quality_string and not length $quality_string;
    my ($highest) = sort { $b <=> $a } map { ord($_) - 33 } split('', $quality_string);
    return $highest;
}

sub minimum_quality {
    my ($self, $quality_string) = @_;
    Carp::confess('No quality string to get minimum quality.') if not defined $quality_string and not length $quality_string;
    my ($lowest) = sort { $a <=> $b } map { ord($_) - 33 } split('', $quality_string);
    return $lowest;
}
#<>#

#<>#
sub config_to_hash {
    my ($self, $config) = @_;

    if ( not $config ) {
        $self->error_message('No config to convert to hash!');
        return;
    }

    my %hash;
    for my $part ( split(':', $config) ) {
        my ($key, $value) = split('=', $part);
        if ( not defined $value ) {
            $value = $key;
            $key = 'file';
        }
        if ( exists $hash{$key} ) { # could support...
            $self->error_message("Duplicate key ($key) in config! $config");
            return;
        }
        $hash{$key} = $value;
    }

    return %hash;
}


sub hash_to_config {
    my ($self, %hash) = @_;

    if ( not %hash ) {
        $self->error_message('No hash to convert to config!');
        return;
    }

    my @configs;
    for my $key ( sort { $a cmp $b } keys %hash ) {
        if ( my $ref = ref $hash{$key} ) { # could support...
            $self->error_message("Unsupported data type ($ref) in hash! ".Data::Dumper::Dumper(\%hash));
            return;
        }
        push @configs, join('=', $key, $hash{$key});
    }

    return join(':', @configs);
}
#<>#

#< Complement >#
sub complement_sequence {
    my ($self, $seq) = @_;

    Carp::confess('No sequence to complement.') if not $seq;

    $seq->{seq} =~ tr/[ACGTacgt]/[TGCAtgca]/;
    $seq->{seq} = reverse $seq->{seq};
    $seq->{qual} = reverse $seq->{qual};

    return $seq;
}
#<>#

1;

