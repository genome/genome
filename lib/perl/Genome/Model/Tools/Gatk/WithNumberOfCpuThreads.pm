package Genome::Model::Tools::Gatk::WithNumberOfCpuThreads;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gatk::WithNumberOfCpuThreads {
    is => 'UR::Object',
    is_abstract => 1,
    subclass_description_preprocessor => 'Genome::Model::Tools::Gatk::WithNumberOfCpuThreads::_preprocess_subclass_description',
};

sub _preprocess_subclass_description {
    my ($class, $desc) = @_;

    my $is = $desc->{is};
    my $base_class = 'Genome::Model::Tools::Gatk::Base';
    unless (grep { $_->isa($base_class) } @$is) {
        Carp::confess(__PACKAGE__ . ' can only be used with subclasses of ' . $base_class);
    }

    my $prop = 'number_of_cpu_threads';
    return $desc if exists $desc->{has}{$prop};

    $desc->{has}{$prop} = {
        is => 'Number',
        gatk_param_name => '-nct',
        doc => 'Controls the number of CPU threads allocated to each data thread',
        is_optional => 1,
        is_input => 1,
    };

    return $desc;
}

1;

