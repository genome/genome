package Genome::InstrumentData::Gatk::IndelRealignerResult;

use strict;
use warnings;

use Genome;

require File::Temp;

# realigner target creator
#  bam [original]
#  ref [fasta]
#  known_sites
#  > intervals
#
# indel realigner
#  bam [original]
#  ref [fasta]
#  known_sites
#  intervals [from realigner target crreator]
#  > bam
class Genome::InstrumentData::Gatk::IndelRealignerResult { 
    is => 'Genome::InstrumentData::Gatk::BaseWithKnownSites',
    has_constant => [
        # outputs
        intervals_file  => { 
            is_output => 1,
            calculate_from => [qw/ output_dir bam_source /],
            calculate => q| return $output_dir.'/'.$bam_source->id.'.bam.intervals'; |,
        },
    ]
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $create_targets = $self->_create_targets;
    if ( not $create_targets ) {
        $self->delete;
        return;
    }

    my $realign_indels = $self->_realign_indels;
    if ( not $realign_indels ) {
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

sub _create_targets {
    my $self = shift;
    $self->status_message('Run realigner target creator...');

    my $intervals_file = $self->intervals_file;
    my %target_creator_params = (
        version => $self->version,
        input_bam => $self->input_bam_path,
        reference_fasta => $self->reference_fasta,
        output_intervals => $intervals_file,
    );
    $target_creator_params{known} = $self->known_sites_indel_vcfs if @{$self->known_sites_indel_vcfs};
    $self->status_message('Params: '.Data::Dumper::Dumper(\%target_creator_params));

    my $target_creator = Genome::Model::Tools::Gatk::RealignerTargetCreator->create(%target_creator_params);
    if ( not $target_creator ) {
        $self->error_message('Failed to create realigner target creator!');
        return;
    }
    if ( not eval{ $target_creator->execute; } ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to execute realigner target creator!');
        return;
    }

    if ( not -s $intervals_file ) {
        $self->error_message('Ran target creator, but failed to make an intervals file!');
        return;
    }
    $self->status_message('Intervals file: '.$intervals_file);

    $self->status_message('Run realigner target creator...done');
    return 1;
}

sub _realign_indels {
    my $self = shift;
    $self->status_message('Run indel realigner...');

    my $bam_path = $self->bam_path;
    my %realigner_params = (
        version => $self->version,
        input_bam => $self->input_bam_path,
        reference_fasta => $self->reference_fasta,
        target_intervals => $self->intervals_file,
        output_realigned_bam => $bam_path,
        index_bam => 1,
        target_intervals_are_sorted => 1,
    );
    $realigner_params{known} = $self->known_sites_indel_vcfs if @{$self->known_sites_indel_vcfs};
    $self->status_message('Params: '.Data::Dumper::Dumper(\%realigner_params));

    my $realigner = Genome::Model::Tools::Gatk::IndelRealigner->create(%realigner_params);
    if ( not $realigner ) {
        $self->error_message('Failed to create indel realigner!');
        return;
    }
    if ( eval{ not $realigner->execute; }) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to execute indel realigner!');
        return;
    }

    if ( not -s $bam_path ) {
        $self->error_message('Ran indel realigner, but failed to create the output bam!');
        return;
    }
    $self->status_message('Bam file: '.$bam_path);

    $self->status_message('Run indel realigner...done');
    return 1;
}

1;

