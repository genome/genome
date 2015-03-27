package Genome::Process::WithVcf;

use strict;
use warnings;
use Genome;
use Genome::Utility::Vcf;

class Genome::Process::WithVcf {
    is => 'Genome::Process',
};

sub special_compare_functions {
    my $self = shift;
    my @functions = $self->SUPER::special_compare_functions(@_);
    push @functions, qr(\.vcf$), sub {my ($a, $b) = @_; Genome::Utility::Vcf::compare_vcf($a, $b)};
    return @functions;
}

1;

