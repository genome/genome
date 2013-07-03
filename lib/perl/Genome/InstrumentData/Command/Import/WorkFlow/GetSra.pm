package Genome::InstrumentData::Command::Import::WorkFlow::GetSra;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::GetSra { 
    is => 'Command::V2',
    has_input => [
        working_directory => {
            is => 'Text',
            doc => 'Detination directory for sra.',
        },
        source_sra_path => {
            is => 'Text',
            doc => 'Source bam to get.',
        },
    ],
    has_output => [
        sra_path => {
            calculate_from => [qw/ working_directory source_sra_base_name /],
            calculate => q( return $working_directory.'/'.$source_sra_base_name; ),
            doc => 'Final retrieved sra path.',
        }, 
    ],
    has_optional_calculated => [
        source_sra_base_name => {
            calculate_from => [qw/ source_sra_path /],
            calculate => q( return File::Basename::basename($source_sra_path); ),
        },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Get sra...');

    # TODO
    # add more methods to get
    # md5 verification

    my $get_sra = $self->_get_sra;
    return if not $get_sra;

    $self->status_message('Get sra...done');
    return 1;
}

sub _get_sra {
    my ($self, $fastq) = @_;

    my $source_sra_path = $self->source_sra_path;
    $self->status_message('Source SRA path: '.$source_sra_path);
    my $sra_path = $self->sra_path;
    $self->status_message('Destination SRA path: '.$sra_path);

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $copy_ok = $helpers->copy_file($source_sra_path, $sra_path);
    return if not $copy_ok;

    return 1;
}

1;

