package Genome::Model::Tools::Sx::Rename;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

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

sub create {
    my ($class , %params) = @_;

    my $self = $class->SUPER::create(%params)
        or return;

    my @matches = $self->matches;
    unless ( @matches ) {
        $self->error_message("No match and replace patterns given.");
        $self->delete;
        return;
    }

    my @match_and_replace;
    for my $match ( @matches ) {
        my ($match, $replace) = split('=', $match);
        my $evald_match = eval($match);
        unless ( defined $evald_match ) {
            $self->error_message("Can't compile match ($match) string: $@");
            return;
        }
        push @match_and_replace, [ $evald_match, $replace ];
    }

    $self->_match_and_replace(\@match_and_replace);

    return $self;
}

sub _eval_seqs {
    my ($self, $seqs) = @_;

    for my $seq ( @$seqs ) { 
        MnR: for my $match_and_replace ( @{$self->_match_and_replace} ) {
            if ( $seq->{id} =~ s/$match_and_replace->[0]/$match_and_replace->[1]/g ) {
                last MnR if $self->first_only;
            }
        }
    }

    return 1;
}

1;

