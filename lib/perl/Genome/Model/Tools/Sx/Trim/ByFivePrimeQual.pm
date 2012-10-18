package Genome::Model::Tools::Sx::Trim::ByFivePrimeQual;

use strict;
use warnings;

use Genome;            

use Regexp::Common;

class Genome::Model::Tools::Sx::Trim::ByFivePrimeQual {
    is => 'Genome::Model::Tools::Sx::Base',
    has => [
        quality => {
            is => 'Integer',
            doc => 'The minimum quality of the sequecnce.',
        },    
     ],
};

sub help_brief {
    return "Trim sequence when avg qual from the 5' end falls below the threshold";
}

sub help_detail {
    return "Starting at the beginning of the sequence, calculate average quality one base at a time. When the quality falls below the threshold, the sequence is trimmed to the end from the base that caused the quality to fall below the threshold. If the average quality never exceeds the threshold, the sequence and quality will be set to empty strings. Emtpy sequences are not removed or filtered."; 
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;
    my $quality = $self->quality;
    if ( $quality !~ /^$RE{num}{int}$/ or $quality < 1 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ length /],
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
            my $orig_seq = $seq->{seq};
            my $orig_qual = $seq->{qual};
            my $i = 0;
            while ( 1 ) {
                $i++;
                $seq->{seq} = substr($orig_seq, 0, $i);
                $seq->{qual} = substr($orig_qual, 0, $i);
                if ( Genome::Model::Tools::Sx::Functions->calculate_average_quality($seq->{qual}) < $quality ) { 
                    chop $seq->{seq};
                    chop $seq->{qual};
                    next SEQ;
                }
                if ( length($seq->{seq}) == length($orig_seq) ) {
                    next SEQ;
                }
            }
        }
        return 1;
    }
}

1;

