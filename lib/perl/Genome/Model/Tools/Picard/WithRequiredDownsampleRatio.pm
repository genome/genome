package Genome::Model::Tools::Picard::WithRequiredDownsampleRatio;

use strict;
use warnings;

use UR;
use Genome::Model::Tools::Picard::WithDownsampleRatio;

class Genome::Model::Tools::Picard::WithRequiredDownsampleRatio {
    is => 'UR::Object',
    is_abstract => 1,
    subclass_description_preprocessor => __PACKAGE__ . '::_preprocess_subclass_description',
};

sub _preprocess_subclass_description {
    my $desc = Genome::Model::Tools::Picard::WithDownsampleRatio::_preprocess_subclass_description(@_);
    if ( not $desc ) {
        die 'Failed to process subclass for with required downsample ratio!';
    }

    if ( not exists $desc->{has}{downsample_ratio} ) {
        die 'Downsample ratio does not exist in class description!';
    }

    $desc->{has}{downsample_ratio}->{is_optional} = 0;

    return $desc;
}

sub __errors__ {
    my $self = shift;

    my $downsample_ratio = $self->downsample_ratio;
    my @errors = $self->SUPER::__errors__;

    if ( defined $downsample_ratio ) {
        push @errors, Genome::Model::Tools::Picard::WithDownsampleRatio::__errors__($self); # call with downsample ratio __errors__
    }

    return @errors;
}

1;

