package Genome::Model::Tools::Sx::Limit::ByBases;

use strict;
use warnings;

use Genome;            

use Regexp::Common;

class Genome::Model::Tools::Sx::Limit::ByBases {
    is => 'Genome::Model::Tools::Sx::Limit::Base',
    has => [
        bases => {
            is => 'Number',
            is_optional => 1,
            doc => 'The maximum number of total bases to write. When this amount is exceeded, writing will be concluded.',
        },
    ],
};

sub help_brief { return 'Limit sequences by bases'; }
sub help_detail { help_brief(); }
sub help_synopsis { help_brief(); }

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;

    my $bases = $self->bases;
    if ( not defined $bases or $bases !~ /^$RE{num}{int}$/ or $bases < 1 ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ bases /],
            desc => 'Bases ('.(defined $bases ? $bases : 'NULL').') must be a positive integer greater than 1',
        );
    }

    return @errors;
}

sub _create_limiters { 
    my $self = shift;

    my $bases = $self->bases;
    return unless defined $bases;

    return sub{
        for my $seq ( @{$_[0]} ) { 
            $bases -= length($seq->{seq});
        }
        return ( $bases > 0 ) ? 1 : 0 ;
    };
}

1;

