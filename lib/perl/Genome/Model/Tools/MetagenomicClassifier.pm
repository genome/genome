package Genome::Model::Tools::MetagenomicClassifier;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::MetagenomicClassifier {
    is => 'Command',
};

sub help_brief {
    "Metagenomic classification tools",
}

sub help_synopsis {
    return <<"EOS"
genome-model tools metagenomic-classifier ...    
EOS
}

my @RANKS = (qw/ domain kingdom phylum class order family genus species /); 
sub taxonomic_ranks {
    return @RANKS;
}

sub taxonomic_rank_at {
    return $RANKS[ $_[1] ];
}

sub is_rank_valid {
    my ($class, $rank) = @_;

    unless ( defined $rank ) {
        $class->error_message("No rank given");
        return;
    }

    unless ( grep { $rank eq $_ } @RANKS ) {
        $class->error_message("Rank ($rank) is not in the list of valid ranks: ".join(', ', @RANKS));
        return;
    }


    return 1;
}

my %DOMAINS_AND_KINGDOMS = (
    archaea => [qw/ archaea /],
    bacteria => [qw/ bacteria /],
    eukarya => [qw/ animalia fungi plantae protista /],
);

sub domains {
    return keys %DOMAINS_AND_KINGDOMS;
}

sub is_domain_valid {
    my ($class, $domain) = @_;

    unless ( defined $domain ) {
        $class->error_message("No domain given");
        return;
    }

    unless ( exists $DOMAINS_AND_KINGDOMS{$domain} ) {
        $class->error_message("Domain ($domain) is not in the list of valid domains: ".join(', ', $class->domains));
        return;
    }
    
    return 1;
}

sub kingdoms_for_domain {
    my ($class, $domain) = @_;

    my $valid_domain = $class->is_domain_valid($domain);
    return if not $valid_domain;

    return @{$DOMAINS_AND_KINGDOMS{$domain}};
}

1;

