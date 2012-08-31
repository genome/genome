package Genome::Model::Tools::Sx::Filter::ByAvgQual;

use strict;
use warnings;

use Genome;

require List::Util;
use Regexp::Common;

class Genome::Model::Tools::Sx::Filter::ByAvgQual {
    is  => 'Genome::Model::Tools::Sx::Filter::Base',
    has => [
        quality => {
            is => 'Integer',
            doc => 'The minimum average quality of the sequence.',
        },    
     ],

};

sub help_brief {
    return 'Filter sequences by average quality';
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;
    if ( $self->quality !~ /^$RE{num}{int}$/ or $self->quality < 0 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ quality /],
            desc => 'Quality is not a integer greater than or equal to 0 => '.$self->quality,
        );
    }
    return @errors;
}

sub _create_filters {
    my $self = shift;

    my $quality = $self->quality;
    return sub{
        my $seqs = shift;
        for my $seq ( @$seqs ) {
            return if Genome::Model::Tools::Sx::Functions->calculate_average_quality($seq->{qual}) < $quality;
        }
        return 1;
    }
}

1;

