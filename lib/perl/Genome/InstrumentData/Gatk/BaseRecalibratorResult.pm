package Genome::InstrumentData::Gatk::BaseRecalibratorResult;

use strict;
use warnings;

use Genome;

# recalibrator
#  bam [from indel realigner]
#  ref [fasta]
#  known_sites [knownSites]
#  > grp [gatk report file]
class Genome::InstrumentData::Gatk::BaseRecalibratorResult { 
    is => 'Genome::InstrumentData::Gatk::BaseWithKnownSites',
    has_output => [
        recalibration_table_file => {
            is_output => 1,
            calculate_from => [qw/ output_dir bam_source /],
            calculate => q| return $output_dir.'/'.$bam_source->id.'.bam.grp'; |,
        },
    ],
};

sub resolve_allocation_kilobytes_requested {
    my $self = shift;
    my $kb_requested = -s $self->input_bam_path;
    return int($kb_requested / 1024 * .2);
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    $self->status_message('Bam source: '.$self->bam_source->id);
    $self->status_message('Reference: '.$self->reference_build->id);
    for my $known_sites ( $self->known_sites ) {
        $self->status_message('Known sites: '.$known_sites->id);
    }

    my $run_recalibrator = $self->_run_base_recalibrator;
    if ( not $run_recalibrator ) {
        $self->delete;
        return;
    }

    my $allocation = $self->disk_allocations;
    eval { $allocation->reallocate };

    return $self;
}

sub _run_base_recalibrator {
    my $self = shift;
    $self->status_message('Run base recalibrator...');

    my $recalibration_table_file = $self->recalibration_table_file;
    my %base_recalibrator_params = (
        version => $self->version,
        input_bam => $self->input_bam_path,
        reference_fasta => $self->reference_fasta,
        output_recalibration_table => $recalibration_table_file,
    );
    $base_recalibrator_params{known_sites} = $self->known_sites_vcfs if @{$self->known_sites_vcfs};
    $self->status_message('Params: '.Data::Dumper::Dumper(\%base_recalibrator_params));

    my $base_recalibrator = Genome::Model::Tools::Gatk::BaseRecalibrator->create(%base_recalibrator_params);
    if ( not $base_recalibrator ) {
        $self->error_message('Failed to create base recalibrator creator!');
        return;
    }
    if ( not eval{ $base_recalibrator->execute; } ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to execute base recalibrator creator!');
        return;
    }

    if ( not -s $recalibration_table_file ) {
        $self->error_message('Ran base recalibrator creator, but failed to make a recalibration table file!');
        return;
    }
    $self->status_message('Recalibration table file: '.$recalibration_table_file);

    $self->status_message('Run base recalibrator...done');
    return 1;
}

1;

