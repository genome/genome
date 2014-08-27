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
    $self->debug_message('Run realigner target creator...');

    my $target_creator = $self->_create_realigner_target_creator_with_threading;
    return if not $target_creator;
    my $rv = eval { $target_creator->execute; };
    if ( not $rv ) {
        my $error_message = "$@";
        $self->error_message($error_message);
        my $exit_code = $self->_resolve_exit_code_from_shellcmd_error_message($error_message);
        if ( $exit_code ne '134' ) {
            # 134 is typically for seg faults
            $self->error_message('Failed to execute realigner target creator with threading!');
            return;
        }
        # Try again w/o threading
        $self->debug_message("Recieved exit code $exit_code, which may be a seg fault due to ta threading issue. Attempting to rerun realigner target creator without threading.");
        my $target_creator = $self->_create_realigner_target_creator;
        return if not $target_creator;
        $rv = eval{ $target_creator->execute; };
        if ( not $rv ) {
            $self->error_message('Failed to execute realigner target creator!');
            return;
        }
    }

    my $intervals_file = $target_creator->output_intervals;
    $self->debug_message('Intervals file: '.$intervals_file);
    if ( not -s $intervals_file ) {
        $self->error_message('Ran target creator, but failed to make an intervals file!');
        return;
    }

    $self->debug_message('Run realigner target creator...done');
    return 1;
}

sub _create_realigner_target_creator_with_threading {
    my $self = shift;
    my $target_creator = $self->_create_realigner_target_creator;
    $target_creator->number_of_threads(8);
    return $target_creator;
}

sub _create_realigner_target_creator {
    my $self = shift;

    my $intervals_file = $self->intervals_file;
    my %target_creator_params = (
        version => $self->version,
        input_bam => $self->input_bam_path,
        reference_fasta => $self->reference_fasta,
        output_intervals => $intervals_file,
        max_memory => $self->max_memory_for_gmt_gatk,
    );
    $target_creator_params{known} = $self->known_sites_indel_vcfs if @{$self->known_sites_indel_vcfs};
    $self->debug_message('Params: '.Data::Dumper::Dumper(\%target_creator_params));

    my $target_creator = Genome::Model::Tools::Gatk::RealignerTargetCreator->create(%target_creator_params);
    if ( not $target_creator ) {
        $self->error_message('Failed to create realigner target creator!');
        return;
    }

    return $target_creator;
}

sub _resolve_exit_code_from_shellcmd_error_message {
    my ($self, $error_message) = @_;

    if ( $error_message =~ /^ERROR RUNNING COMMAND.  Exit code (\d+) from/ ) {
        return $1;
    }

    return 1;
}

sub _realign_indels {
    my $self = shift;
    $self->debug_message('Run indel realigner...');

    my $bam_path = $self->bam_path;
    my %realigner_params = (
        version => $self->version,
        input_bam => $self->input_bam_path,
        reference_fasta => $self->reference_fasta,
        target_intervals => $self->intervals_file,
        output_realigned_bam => $bam_path,
        index_bam => 1,
        target_intervals_are_sorted => 1,
        max_memory => $self->max_memory_for_gmt_gatk,
    );
    $realigner_params{known} = $self->known_sites_indel_vcfs if @{$self->known_sites_indel_vcfs};
    $self->debug_message('Params: '.Data::Dumper::Dumper(\%realigner_params));

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
    $self->debug_message('Bam file: '.$bam_path);

    $self->debug_message('Run indel realigner...done');
    return 1;
}

1;

