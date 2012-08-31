package Genome::Model::Tools::DetectVariants2::Filter::VarscanPValue;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::Filter::VarscanPValue {
    is => 'Genome::Model::Tools::DetectVariants2::Filter',
    doc => 'This module filters out variants that do not meet a p-value cutoff.',
    has => [
        variant_p_value => {
            is => 'Number',
            doc => 'Keep variants that have at least this significant a variant p-value',
            is_optional => 1,
        },
        somatic_p_value => {
            is => 'Number',
            doc => 'Keep variants that have at least this significant a somatic p-value',
            is_optional => 1,
        },
    ],
};

sub help_synopsis {
    'gmt detect-variants2 filter varscan-p-value --variant-p-value 0.01 --somatic-p-value 0.02 ...'
}

sub help_detail {
    <<'EOS'
A simple filter on the p-value of the VarScan call.  This filter is used
as part of a strategy for the dispatcher, e.g.:
    "filtered by varscan-p-value v1 [--variant-p-value 0.01 --somatic-p-value 0.02]"
EOS
;
}

sub _variant_type { 'snvs' };

sub _filter_variants {
    my $self = shift;

    my $input_file = join('/',  $self->input_directory, 'snvs.hq');
    my $hq_output_file = join('/', $self->_temp_staging_directory, 'snvs.hq');
    my $lq_output_file = join('/', $self->_temp_staging_directory, 'snvs.lq');

    my $input_fh = Genome::Sys->open_file_for_reading($input_file);
    my $hq_output_fh = Genome::Sys->open_file_for_writing($hq_output_file);
    my $lq_output_fh = Genome::Sys->open_file_for_writing($lq_output_file);

    while (my $line = <$input_fh>) {
        chomp $line;
        my @fields = split "\t", $line;
        my $variant_p_value = $fields[13];
        my $somatic_p_value = $fields[14];

        my $pass = 1;
        if(defined $self->variant_p_value and $variant_p_value > $self->variant_p_value ) { $pass = 0; }
        if(defined $self->somatic_p_value and $somatic_p_value > $self->somatic_p_value) { $pass = 0; }

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

sub _create_bed_file {
    my $self = shift;
    my $detector_file = shift;
    my $bed_file = shift;

    my $convert = Genome::Model::Tools::Bed::Convert::Snv::VarscanSomaticToBed->create(
                    source => $detector_file,
                    output => $bed_file,
                    reference_build_id => $self->reference_build_id,
    );
    unless($convert->execute){
        die $self->error_message("Failed to convert detector to bed");
    }

    return 1;
}

1;
