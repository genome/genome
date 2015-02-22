package Genome::Model::Tools::DetectVariants2::Result::Base;

use strict;
use warnings;

use Genome;
use Carp qw(confess);

class Genome::Model::Tools::DetectVariants2::Result::Base {
    is => ['Genome::SoftwareResult::Stageable'],
    is_abstract => 1,
    doc => 'This class represents the result of a detect-variants operation. This base class just unites the various result types',
};

sub path {
    my $self = shift;
    my ($str) = @_;

    return join('/', $self->output_dir, $str);
}

sub get_vcf_result {
    my $self = shift;

    confess sprintf("Method 'get_vcf_result' is abstract and not defined for class (%s)",
        $self->class);
}

1;
