package Genome::Model::Tools::Picard::Downsample;

use strict;
use warnings;

use Genome;
use IO::File;
use File::Basename;


class Genome::Model::Tools::Picard::Downsample {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input_file => {
            is  => 'String',
            doc => 'The SAM/BAM file to downsample.',
            is_optional => 0,
        },
        output_file => {
            is  => 'String',
            doc => 'The resulting downsampled SAM/BAM file.  File type is determined by suffix.',
            is_optional => 0,
        },
        downsample_ratio => {
            is => 'String',
            doc => 'ratio at which to keep reads in order to downsample. 0.01 means keep 1 in 100 reads. ',
            is_optional => 0,
        },
        use_version => {
            is => 'String',
            doc => 'Version must be 1.52 or greater, default is 1.52',
            default => '1.52',
        },
        max_records_in_ram => {
            is => 'String',
            doc => 'Set this to control the number of records in ram. Defaults to 500,000',
            default => '500000',
        },
        random_seed => {
            is => 'String',
            doc => 'Set this to attain reproducability',
            is_optional => 1,
        },
    ],
};

sub help_detail {
    return <<EOS
    Tool to downsample a BAM or SAM file using Picard.  For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#DownsampleSam
EOS
}

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path .'/DownsampleSam.jar';
    unless (-e $jar_path) {
        die('Failed to find '. $jar_path .'!  This command may not be available in version '. $self->use_version);
    }
    my $input_file = $self->input_file;
    my $output_file = $self->output_file;
    my $ratio = $self->downsample_ratio;
    my $sort_cmd = $jar_path .' net.sf.picard.sam.DownsampleSam O='. $output_file .' I='. $input_file .' PROBABILITY='.$ratio.' MAX_RECORDS_IN_RAM='. $self->max_records_in_ram;
    if(defined($self->random_seed)){
        my $seed = int($self->random_seed);
        $sort_cmd .= " RANDOM_SEED=".$seed;
    }
    $self->run_java_vm(
        cmd => $sort_cmd,
        input_files => [$input_file],
        output_files => [$output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
