package Genome::Model::Tools::Sx::Trim::ByQual;

use strict;
use warnings;

use Genome;            

use Regexp::Common;

class Genome::Model::Tools::Sx::Trim::ByQual {
    is => 'Genome::Model::Tools::Sx::Base',
    has => [
        quality => {
            is => 'Integer',
            doc => 'Trim bases from the right end if the quality is this value or less.',
        },    
     ],
};

sub help_brief {
    return "Trim right end if qual is less than/equal to the threshold";
}

sub help_detail {
    return "Trim bases from the right end of a sequence if the quality is less than or equal to the quality threshold. If all bases are trimmed, the sequence and quality will be set to empty strings. Use a filter to remove unwanted sequences (gmt sx filter -h)."; 
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;
    my $quality = $self->quality;
    if ( $quality !~ /^$RE{num}{int}$/ or not $quality >= 0 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ quality /],
            desc => 'Quality is not a integer greater than or equal to 0 => '.$quality,
        );
    }
    return @errors;
}

sub _eval_seqs {
    my ($self, $seqs) = @_;

    for my $seq (@$seqs) {
        my $qual = chop $seq->{qual};
        while ( $qual ne '' and Genome::Model::Tools::Sx::Functions->calculate_quality($qual) <= $self->quality ) {
            chop $seq->{seq};
            $qual = chop $seq->{qual};
        }
        $seq->{qual} .= $qual if defined $qual;
    }

    return 1;
}

1;

