package Genome::Model::Tools::Sx::Trim::ByAvgQual;

use strict;
use warnings;

use Genome;            

use Regexp::Common;

class Genome::Model::Tools::Sx::Trim::ByAvgQual {
    is => 'Genome::Model::Tools::Sx::Base',
    has => [
        quality => {
            is => 'Integer',
            doc => 'The minimum quality of the entire seqeunce. Bases will be trimmed from the end until quality reaches this average.',
        },    
     ],
};

sub help_brief {
    return "Trim sequence from 3' end until avg qual is above the threshold";
}

sub help_detail {
    return "Trim one base at a time from the 3' end until the average quality of the sequence reaches or exceeds the quality threshold. If the average qualtity of the sequence does not reach or execeed the threshold, the sequence and quality will be set to empty strings. Emtpy sequences are not removed or filtered."; 
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;
    my $quality = $self->quality;
    if ( $quality !~ /^$RE{num}{int}$/ or $quality < 1 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ quality /],
            desc => 'Quality is not a integer greater than 0 => '.$quality,
        );
    }
    return @errors;
}

sub _create_evaluator {
    my $self = shift;

    my $quality = $self->quality;
    return sub{
        SEQ: for my $seq ( @{$_[0]} ) {
            while ( Genome::Model::Tools::Sx::Functions->calculate_average_quality($seq->{qual}) < $quality ) {
                chop $seq->{seq};
                chop $seq->{qual};
                next SEQ if length($seq->{seq}) == 0;
            }
        }
        return 1;
    }
}

1;

