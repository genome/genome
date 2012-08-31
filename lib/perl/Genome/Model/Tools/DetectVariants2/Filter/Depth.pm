package Genome::Model::Tools::DetectVariants2::Filter::Depth;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::Filter::Depth {
    is => 'Genome::Model::Tools::DetectVariants2::Filter',
    doc => 'This module filters out variants that fall outside specified depth thresholds.',
    has_input => [
        minimum_depth => {
            is => 'Integer',
            doc => 'Keep variants that have at least a depth of this value',
            is_optional => 1,
        },
        maximum_depth => {
            is => 'Integer',
            doc => 'Keep variants that have at most a depth of this value',
            is_optional => 1,
        },
    ],
};

sub help_synopsis {
    'gmt detect-variants2 filter depth --minimum-depth 30 --maximum-depth 10000 ...'
}

sub help_detail {
    <<'EOS'
A simple filter on the depth of the read-support for a variant.  This filter is used
as part of a strategy for the dispatcher, e.g.:
    "filtered by depth v1 [--minimum-depth 30 --maximum-depth 10000]"
EOS
;
}

sub _variant_type { 'snvs' };

sub _filter_variants {
    my $self = shift;

    my $input_file = join('/',  $self->input_directory, 'snvs.hq.bed');
    my $hq_output_file = join('/', $self->_temp_staging_directory, 'snvs.hq.bed');
    my $lq_output_file = join('/', $self->_temp_staging_directory, 'snvs.lq.bed');

    my $input_fh = Genome::Sys->open_file_for_reading($input_file);
    my $hq_output_fh = Genome::Sys->open_file_for_writing($hq_output_file);
    my $lq_output_fh = Genome::Sys->open_file_for_writing($lq_output_file);

    while (my $line = <$input_fh>) {
        chomp $line;
        my @fields = split "\t", $line;
        my $depth = $fields[5];

        my $pass = 1;
        if(defined $self->minimum_depth and $depth < $self->minimum_depth) { $pass = 0; }
        if(defined $self->maximum_depth and $depth > $self->maximum_depth) { $pass = 0; }

        if($pass) {
            $hq_output_fh->print($line, "\n");
        } else {
            $lq_output_fh->print($line, "\n");
        }
    }

    close($input_fh);
    close($hq_output_fh);
    close($lq_output_fh);

    return 1;
}

1;
