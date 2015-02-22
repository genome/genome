package Genome::Model::Tools::Picard::WithDownsampleRatio;

use strict;
use warnings;

use UR;

require List::MoreUtils;
use Regexp::Common;

class Genome::Model::Tools::Picard::WithDownsampleRatio {
    is => 'UR::Object',
    is_abstract => 1,
    subclass_description_preprocessor => __PACKAGE__ . '::_preprocess_subclass_description',
};

sub _preprocess_subclass_description {
    my ($class, $desc) = @_;

    my $class_name = $desc->{class_name};
    my $downsample_ratio_property = {
        class_name => $class_name,
        property_name => 'downsample_ratio',
        is => 'Float',
        doc => 'Ratio of reads to keep. A value of 0.01 means keep 1 in 100 reads. Specify a value to the hundreths place (0.00) that is greater than 0, and less than 1.',
    };
    $downsample_ratio_property->{is_optional} = 1 unless $class_name->can('_downsample_ratio_is_required')
        and $class_name->_downsample_ratio_is_required;
    $downsample_ratio_property->{is_input} = 1 if grep { $_ eq 'Command' or $_ eq 'Command::V2' } @{$desc->{is}};

    $desc->{has}{downsample_ratio} = $downsample_ratio_property;

    return $desc;
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    my $downsample_ratio_error = List::MoreUtils::any { $_ eq 'downsample_ratio' } map { $_->properties } @errors;

    my $downsample_ratio = $self->downsample_ratio;
    if ( not $downsample_ratio_error and ( $downsample_ratio <= 0 or $downsample_ratio >= 1 ) ) {
        push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ downsample_ratio /],
                desc => 'Must be greater than 0 and less than 1! '.$downsample_ratio,
        );
    }

    return @errors;
}

1;

