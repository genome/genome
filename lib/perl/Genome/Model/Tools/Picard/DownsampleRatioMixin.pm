package Genome::Model::Tools::Picard::DownsampleRatioMixin;

use strict;
use warnings;

use UR;

use Regexp::Common;

class Genome::Model::Tools::Picard::DownsampleRatioMixin {
};

sub downsample_ratio_property {
    return {
        is => 'Float',
        doc => 'Ratio of reads to keep. A value of 0.01 means keep 1 in 100 reads. Specify a value to the hundreths place (0.00) that is greater than 0, and less than 1.',
    };
}

sub __errors__ {
    my $obj_with_downsample_ratio = shift;

    my $downsample_ratio = $obj_with_downsample_ratio->downsample_ratio;
    return if not defined $downsample_ratio;

    if ( $downsample_ratio !~ /$RE{num}{real}/ ) {
        return (
            UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ downsample_ratio /],
                desc => 'Invalid number! '.$downsample_ratio,
            )
        );
    }

    if ( $downsample_ratio <= 0 or $downsample_ratio >= 1 ) {
        return (
            UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ downsample_ratio /],
                desc => 'Must be greater than 0 and less than 1! '.$downsample_ratio,
            )
        );
    }

    # set to the hundreths place
    $obj_with_downsample_ratio->downsample_ratio( sprintf('%.2f', $downsample_ratio) );

    return;
}

1;

