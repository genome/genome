package Genome::Process::WithVcf;

use strict;
use warnings;
use Genome;
use Genome::Utility::Vcf;

class Genome::Process::WithVcf {
    is => 'Genome::Process',
};

sub _get_handler_for_file {
    my $self = shift;
    my $file = shift;
    if (is_vcf($file)) {
        return "compare_vcf";
    }

    return $self->SUPER::_get_handler_for_file($file);
}

sub compare_vcf {
    my ($self, $first_file, $second_file) = @_;
    return Genome::Utility::Vcf::compare_vcf($first_file, $second_file);
}

sub is_vcf {
   my $file = shift;
   return $file =~ m/\.vcf$/;
}

1;

