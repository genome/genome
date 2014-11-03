package Genome::VariantReporting::Framework::Test::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Basename qw(dirname);

class Genome::VariantReporting::Framework::Test::Adaptor {
    is => 'Genome::VariantReporting::Framework::Component::Adaptor',

    has_planned_output => [
        __planned__ => {},
    ],
    has_output => [
        __provided__ => {is_many => 1},
    ],
};

sub name {
    "__test__";
}

sub resolve_plan_attributes {
    my $self = shift;
    $self->SUPER::resolve_plan_attributes();
    my $input_vcf = dirname(__FILE__)."/Expert.t.d/input.vcf.gz";
    $self->__provided__([$input_vcf, $input_vcf]);
}

1;
