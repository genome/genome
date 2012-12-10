package Genome::Utility::MetagenomicClassifier;

use strict;
use warnings;

use Bio::Taxon;

sub error_message {
    my ($self, $msg) = @_;

    printf(
        "%s [ERROR] %s\n",
        __PACKAGE__,
        $msg,
    );

    return $msg;
}

my @RANKS = (qw/ domain kingdom phylum class order family genus species /); 
sub taxonomic_ranks {
    return @RANKS;
}

sub taxonomic_rank_at {
    return $RANKS[ $_[1] ];
}

sub validate_taxonomic_rank {
    my ($self, $rank) = @_;

    unless ( defined $rank ) {
        $self->error_message("No rank given");
        return;
    }

    unless ( grep { $rank eq $_ } @RANKS ) {
        $self->error_message("Rank ($rank) is not in the list of valid ranks: ".join(', ', @RANKS));
        return;
    }


    return 1;
}

my @DOMAINS = (qw/ archaea bacteria eukarya/);
sub domains {
    return @DOMAINS;
}

sub validate_domain {
    my ($self, $domain) = @_;

    unless ( defined $domain ) {
        $self->error_message("No domain given");
        return;
    }

    unless ( grep { $domain =~ /^$_$/i } $self->domains ) {
        $self->error_message("Domain ($domain) is not in the list of valid domains: ".join(', ', @DOMAINS));
        return;
    }
    
    return 1;
}

my @KINGDOMS = (qw/ animalia archaea eubacteria fungi plantae protista /);
sub kingdoms {
    return @KINGDOMS;
}

sub create_taxon {
    my ($self, %params) = @_;

    $params{id} =~ s/['"]//g;
    
    my $taxon = Bio::Taxon->new(
        '-id' => $params{id},
        '-rank' => $params{rank},
    );

    if ( $params{tags} ) {
        for my $key ( keys %{$params{tags}} ) {
            $taxon->add_tag_value($key, $params{tags}->{$key});
        }
    }

    $params{ancestor}->add_Descendent($taxon) if $params{ancestor};

    return $taxon;
}

1;

=pod

=head1 Name

ModuleTemplate

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

