package Genome::Model::Tools::Sx::Rename;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Rename {
    is  => 'Genome::Model::Tools::Sx::Base',
    has => [
        matches => {
            is => 'Text',  
            is_many => 1,
            is_input => 1,
            doc => 'Match patterns and replace strings. Match strings can be regular text or regular expressions incased in qr//. Separate match and replace string w/ an equals (=). Separate multiple matches with commas.',
        },
        first_only => {
            is => 'Boolean',
            is_input => 1,
            is_optional => 1,
            default_value => 0,
            doc => 'If given multiple matches, stop after one match/replace is successful.',
        },
        _match_and_replace => {
            is => 'ARRAY',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    return 'Rename sequences';
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;
    
    my @match_and_replace;
    MATCH: for my $match ( $self->matches ) {
        my ($match, $replace) = split('=', $match);
        my $evald_match = eval($match);
        unless ( defined $evald_match ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ matches /],
                desc => "Failed to compile match ($match) string: $@",
            );
            next MATCH;
        }
        push @match_and_replace, [ $evald_match, $replace ];
    }
    $self->_match_and_replace(\@match_and_replace) if not @errors;

    return @errors;
}


sub _create_evaluator {
    my $self = shift;

    my @match_and_replace = @{$self->_match_and_replace};
    my $first_only = $self->first_only;
    return sub{
        for my $seq ( @{$_[0]} ) {
            MnR: for my $match_and_replace ( @match_and_replace ) {
                if ( $seq->{id} =~ s/$match_and_replace->[0]/$match_and_replace->[1]/g ) {
                    last MnR if $first_only;
                }
            }
        }
        return 1;
    }
}

1;

