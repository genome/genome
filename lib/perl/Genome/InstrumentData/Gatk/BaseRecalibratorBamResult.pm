package Genome::InstrumentData::Gatk::BaseRecalibratorBamResult;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Gatk::BaseRecalibratorResult;
require File::Path;

# recalibrator result
#  bam [from indel realigner]
#  ref [fasta]
#  known_sites [knownSites]
#  > grp [gatk report file]
#
# print reads
#  bam [from indel realigner]
#  ref [fasta]
#  grp [from recalibrator]
#  > bam
class Genome::InstrumentData::Gatk::BaseRecalibratorBamResult { 
    is => 'Genome::InstrumentData::Gatk::BaseWithKnownSites',
    has_transient_optional => [
        base_recalibrator_result => { is => 'Genome::InstrumentData::Gatk::BaseRecalibratorResult', },
    ],
};

sub resolve_allocation_kilobytes_requested {
    my $self = shift;
    my $kb_requested = -s $self->input_bam_path;
    return int($kb_requested / 1024 * 7);
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $base_recalibrator_result = $self->_get_or_create_base_recalibrator_result;
    if ( not $base_recalibrator_result ) {
        $self->delete;
        return;
    }
    $base_recalibrator_result->add_user(user => $self, label => 'recalibration table');

    my $print_reads = $self->_print_reads;
    if ( not $print_reads ) {
        $self->delete;
        return;
    }

    my $run_flagstat = $self->run_flagstat_on_output_bam_path;
    if ( not $run_flagstat ) {
        $self->delete;
        return;
    }

    my $run_md5sum = $self->run_md5sum_on_output_bam_path;
    if ( not $run_md5sum ) {
        $self->delete;
        return;
    }

    my $allocation = $self->disk_allocations;
    eval { $allocation->reallocate };

    return $self;
}

sub _get_or_create_base_recalibrator_result {
    my $self = shift;
    $self->status_message('Get or create base recalibrator result...');

    my %base_recalibrator_params = $self->base_recalibrator_params;
    my $base_recalibrator_result = Genome::InstrumentData::Gatk::BaseRecalibratorResult->get_or_create(%base_recalibrator_params);
    if ( not $base_recalibrator_result ) {
        $self->error_message('Failed to get or create base recalibrator result!');
        return;
    }
    $self->base_recalibrator_result($base_recalibrator_result);

    my $recalibration_table_file = $base_recalibrator_result->recalibration_table_file;
    if ( not -s $recalibration_table_file ) {
        $self->error_message('Got base recalibrator result, but failed to find the recalibration table file!');
        return;
    }
    $self->status_message('Recalibration table file: '.$recalibration_table_file);

    $self->status_message('Get or create base recalibrator result...done');
    return $base_recalibrator_result;
}

sub base_recalibrator_params {
    my $self = shift;

    my %base_recalibrator_params = (
        version => $self->version,
        bam_source => $self->bam_source,
        reference_build => $self->reference_build,
    );
    my @known_sites = $self->known_sites;
    $base_recalibrator_params{known_sites} = \@known_sites if @known_sites;

    return %base_recalibrator_params;
}

sub get_base_recalibrator_result {
    my $self = shift;

    return $self->base_recalibrator_result if $self->base_recalibrator_result;

    my %base_recalibrator_params = $self->base_recalibrator_params;
    my $base_recalibrator_result = Genome::InstrumentData::Gatk::BaseRecalibratorResult->get_with_lock(%base_recalibrator_params);
    return if not $base_recalibrator_result;

    return $self->base_recalibrator_result($base_recalibrator_result);
}

sub _print_reads {
    my $self = shift;
    $self->status_message('Print reads...');

    my $tmp_dir = $self->output_dir.'/tmp';
    my $mkdir = eval{ Genome::Sys->create_directory($tmp_dir); };
    if ( not $mkdir ) {
        $self->error_message('Failed to maketmp directory!');
        return;
    }

    my $bam_path = $self->bam_path;
    my $print_reads = Genome::Model::Tools::Gatk::PrintReads->create(
        version => $self->version,
        input_bams => [ $self->input_bam_path ],
        reference_fasta => $self->reference_fasta,
        output_bam => $bam_path,
        bqsr => $self->base_recalibrator_result->recalibration_table_file,
        tmp_dir => $tmp_dir,
    );
    if ( not $print_reads ) {
        $self->error_message('Failed to create print reads!');
        return;
    }
    if ( not $print_reads->execute ) {
        $self->error_message('Failed to execute print reads!');
        return;
    }

    if ( not -s $bam_path ) {
        $self->error_message('Ran print reads, but failed to create the output bam!');
        return;
    }
    $self->status_message('Bam file: '.$bam_path);

    File::Path::rmtree($tmp_dir);

    $self->status_message('Print reads...done');
    return 1;
}

1;

