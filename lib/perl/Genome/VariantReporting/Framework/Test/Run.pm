package Genome::VariantReporting::Framework::Test::Run;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Basename qw(dirname);

class Genome::VariantReporting::Framework::Test::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_planned_transient => [
        __planned__ => {},
    ],
    has_transient => [
        __input__ => {
            is_many => 1,
        },
    ],
};

sub name {
    '__test__';
}

sub result_class {
    'Genome::VariantReporting::Framework::Test::RunResult';
}

sub resolve_plan_attributes {
    my $self = shift;
    $self->SUPER::resolve_plan_attributes();
    my $input_vcf = dirname(__FILE__)."/Expert.t.d/input.vcf.gz";
    my $tmp_vcf = Genome::Sys->create_temp_file_path;
    Genome::Sys->copy_file($input_vcf, $tmp_vcf);
    $self->__input__([$tmp_vcf, $tmp_vcf]);
}

1;
