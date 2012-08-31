package Genome::Model::GenotypeMicroarray::Filter::ByWhitelist;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::Filter::ByWhitelist {
    has => [
        whitelist_snps_file => { is => 'Text' },
    ],
    has_transient_optional => [
        whitelist => { is => 'HASH', default => {} },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    my $snp_input = new FileHandle ($self->whitelist_snps_file);
    unless($snp_input) {
        $self->error_message("Unable to open ".$self->whitelist_snps_file);
        return;
    }

    while (my $line = <$snp_input>) {
        chomp($line);
        my ($id) = split(/\t/, $line);
        $self->whitelist->{$id} = 1;
    }

    return $self;
}

sub filter {
    my ($self, $variant) = @_;
    return 1 if $self->whitelist->{ $variant->{id} };
    return;
}

1;
