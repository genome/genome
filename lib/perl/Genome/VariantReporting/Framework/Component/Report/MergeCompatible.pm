package Genome::VariantReporting::Framework::Component::Report::MergeCompatible;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Framework::Component::Report::MergeCompatible {
    is => 'Genome::SoftwareResult::StageableSimple',
    is_abstract => 1,
    has_transient_optional => {
        has_size => {
            is => 'Boolean',
        },
    },
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

sub has_size {
    my $self = shift;

    return unless $self->report_path;

    unless (defined($self->__has_size)) {
        $self->__has_size(-s $self->report_path);
    }
    return $self->__has_size;
}

1;
