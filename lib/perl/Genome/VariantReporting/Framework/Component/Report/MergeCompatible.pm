package Genome::VariantReporting::Framework::Component::Report::MergeCompatible;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Framework::Component::Report::MergeCompatible {
    is => 'Genome::SoftwareResult::StageableSimple',
    is_abstract => 1,
};

sub can_be_merged {
    return 0;
}

sub merge_parameters {
    return {};
}

sub report_path {
    my $self = shift;
    die sprintf("Method 'report_path' must be overridden in class %s",
        $self->class);
}

1;
