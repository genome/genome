package Genome::Model::Tools::Gatk::WithNumberOfThreads;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gatk::WithNumberOfThreads {
    is => 'UR::Object',
    is_abstract => 1,
    subclass_description_preprocessor => 'Genome::Model::Tools::Gatk::WithNumberOfThreads::_preprocess_subclass_description',
};

sub _preprocess_subclass_description {
    my ($class, $desc) = @_;

    my $is = $desc->{is};
    my $base_class = 'Genome::Model::Tools::Gatk::Base';
    unless (grep { $_->isa($base_class) } @$is) {
        Carp::confess(__PACKAGE__ . ' can only be used with subclasses of ' . $base_class);
    }

    my $prop = 'number_of_threads';
    return $desc if exists $desc->{has}{$prop};

    $desc->{has}{$prop} = {
        is => 'Number',
        gatk_param_name => '-nt',
        doc => 'Controls the number of data threads sent to the processor',
        is_optional => 1,
        is_input => 1,
    };

    return $desc;
}

1;

