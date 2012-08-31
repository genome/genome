package Genome::Model::Tools::Sx::Trim::BwaStyle;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use Regexp::Common;

class Genome::Model::Tools::Sx::Trim::BwaStyle {
    is  => 'Genome::Model::Tools::Sx::Base',
    has => [
        trim_qual_level => {
            is  => 'Integer',
            doc => 'Trim quality level.',
        },
    ],
};

sub help_brief {
    return 'Trim bwa style';
}

sub help_synopsis {
    return <<HELP
HELP
}

sub help_detail {
    return <<HELP 
HELP
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;

    my $trim_qual_level = $self->trim_qual_level;
    if ( not $trim_qual_level or $trim_qual_level !~ /^$RE{num}{int}$/ or $trim_qual_level < 0 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ trim_qual_level /],
            desc => 'Trim qual level must be a positive integer: '.$self->trim_qual_level,
        );
    }

    return @errors;
}

sub _eval_seqs {
    my ($self, $seqs) = @_;

    # sanger: qual_str => # qual_thresh 33
    # illumina: qual_str => B qual_thresh => 64
    for my $seq ( @$seqs ) {
        my $trimmed_length;
        my ($pos, $maxPos, $area, $maxArea) = (length $seq->{seq}, length $seq->{seq}, 0, 0);

        while ($pos > 0 and $area >= 0) {
            $area += $self->trim_qual_level - (ord(substr($seq->{qual}, $pos - 1, 1)) - 33);
            if ($area > $maxArea) {
                $maxArea = $area;
                $maxPos = $pos;
            }
            $pos--;
        }

        if ($pos == 0) { 
            # scanned whole read and didn't integrate to zero?  replace with "empty" read ...
            $seq->{seq}  = 'N';
            $seq->{qual} = '#';
            $maxPos = 1;  # reset to leave 1 base/qual as N/# there
        }
        else {  # integrated to zero?  trim before position where area reached a maximum (~where string of qualities were still below 20 ...)
            $seq->{seq}  = substr($seq->{seq},  0, $maxPos);
            $seq->{qual} = substr($seq->{qual}, 0, $maxPos);
        }
    }

    return 1;
}

1;

