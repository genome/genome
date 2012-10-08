package Genome::Model::Tools::Sx::Trim::ByAvgQualWindow;

use strict;
use warnings;

use Genome;            

use Regexp::Common;

class Genome::Model::Tools::Sx::Trim::ByAvgQualWindow {
    is => 'Genome::Model::Tools::Sx::Base',
    has => [
        window => {
            is => 'Integer',
            doc => 'The size of the window.',
        },
        quality => {
            is => 'Integer',
            doc => 'The minimum quality threshold of the window. Bases will be trimmed from the end until quality reaches this average.',
        },    
     ],
};

sub help_brief {
    return "Trim sequence from 3' end until avg qual of a window is above the threshold";
}

sub help_detail {
    return "Trim one base at a time from the 3' end until the average quality of the sequence reaches or exceeeds the quality threshold. If the length of the sequence is is shorter than the window, the sequence and quality will be set to empty strings. Emtpy sequences are not removed or filtered."; 
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;
    for my $property (qw/ quality window /) {
        my $value = $self->$property;
        if ( $value !~ /^$RE{num}{int}$/ or $value < 1 ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [ $property ],
                desc => ucfirst($property).' is not a integer greater than 0 => '.$value,
            );
        } 
    }
    return @errors;
}

sub _create_evaluator {
    my $self = shift;

    my $window = $self->window;
    my $quality = $self->quality;
    return sub{
        SEQ: for my $seq ( @{$_[0]} ) {
            while ( 1 ) {
                if ( length($seq->{seq}) < $window ) { 
                    $seq->{seq} = '';
                    $seq->{qual} = '';
                    next SEQ;
                }
                my $offset = length($seq->{qual}) - $window;
                my $qual = substr($seq->{qual}, $offset, $window);
                my $score = Genome::Model::Tools::Sx::Functions->calculate_average_quality($qual);
                if ( $score >= $quality ) { 
                    next SEQ;
                }
                chop $seq->{seq};
                chop $seq->{qual};
            }
        }
        return 1;
    }
}

1;

