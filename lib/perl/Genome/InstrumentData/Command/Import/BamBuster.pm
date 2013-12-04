package Genome::InstrumentData::Command::Import::BamBuster;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::InstrumentData::Command::Import::BamBuster { 
    is => 'Genome::Command::Base',
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData::Imported',
            doc => 'The ID of the instrument data to bust, and make imported instrument data for each library.',
        },
        lsf_job_group_name => {
            is => 'String',
            is_many => 0,
            doc => 'The job group name, to throttle LSF jobs on the import'
        },
    ],
    has_optional => [
        _sample => {is_transient=>1},
        _bam => {is_transient=>1},
        _bam_size => {is_transient=>1},
        _disk_allocation => {is_transient=>1},
    ],
};

sub help_brief {
    return 'Bust and import the bam of an instrument data';
}

sub execute {
    my $self = shift;

    my $bam_ok = $self->_check_instrument_data;
    return if not $bam_ok;

    my $disk_allocation = $self->_create_disk_allocation_for_busting;
    return if not $disk_allocation;

    my $bust_bams = $self->_bust_bams;
    return if not $bust_bams;

    my %libraries_and_busted_bams = $self->_get_busted_bams_and_create_libraries;
    return if not %libraries_and_busted_bams;

    my $jobs_ok = $self->_launch_jobs_and_monitor(%libraries_and_busted_bams);
    return if not $jobs_ok;

    return 1;
}

sub _check_instrument_data {
    my $self = shift;

    my $instrument_data = $self->instrument_data;
    $self->status_message('Check instrument data');

    $self->status_message('Check sample');
    my $sample_id = $instrument_data->sample_id;
    if ( not defined $sample_id ) {
        $self->error_message('No sample id fo instrument data: '.$instrument_data->id);
        return;
    }
    my $sample = Genome::Sample->get($instrument_data->sample_id);
    if ( not $sample ) {
        $self->error_message('Cannot get instrument data sample for id: '.$instrument_data->sample_id);
        return;
    }
    $self->_sample($sample);
    $self->status_message('Check sample OK');

    $self->status_message('Check BAM');
    my $bam = $instrument_data->data_directory . "/all_sequences.bam";
    my $bam_size = -s $bam;
    if ( not -s $bam ) {
        $self->error_message("Did not find BAM file named 'all_sequences.bam' in main directory of instrument data: ".$instrument_data->data_directory);
        return;
    }
    $self->_bam($bam); 
    $self->_bam_size($bam_size); 
    $self->status_message('BAM OK');

    $self->status_message('Check instrument data OK');

    return 1; 
}

sub _create_disk_allocation_for_busting {
    my $self = shift;

    $self->status_message('Create disk allocation');

    # FIXME WHAT SHOULD OWN THIS ALLOCATION?
    my %disk_allocation_params = (
        disk_group_name     => $ENV{GENOME_DISK_GROUP_ALIGNMENTS},
        allocation_path     => 'alignment_data/busted_bam/' . $self->instrument_data->id,
        kilobytes_requested => sprintf('%.0f',($self->_bam_size * .002)), # reserve 2X the bam size
        owner_class_name    => 'Genome::InstrumentData::Bam',
        owner_id            => $self->instrument_data->id,
    );
    my $disk_allocation = Genome::Disk::Allocation->allocate(%disk_allocation_params);
    if ( not $disk_allocation ) {
        $self->error_message("Failed to create disk allocation with params:\n". Data::Dumper::Dumper(%disk_allocation_params));
        return;
    }
    $self->_disk_allocation($disk_allocation);

    $self->status_message('Create disk allocation OK');

    return 1 ;
}

sub _bust_bams {
    my $self = shift;

    $self->status_message('Bust BAM');

    my $bam = $self->_bam;
    $self->status_message('From: '.$bam);
    my $disk_allocation  = $self->_disk_allocation;
    $self->status_message('To: '.$disk_allocation->absolute_path);

    my $buster = Genome::Model::Tools::Sam::BamBuster->create(
        input => $bam,
        output_directory => $disk_allocation->absolute_path,
    );
    if ( not $buster ) {
        $self->error_message("Failed to create bam buster.");
        return;
    }
    my $execute = $buster->execute;
    if ( not $execute ) {
        $self->error_message("Failed to execute bam buster.");
        return;
    }

    $self->status_message('Bust BAM OK');

    return  1;
}

sub _get_busted_bams_and_create_libraries {
    my $self = shift;

    $self->status_message('Get busted bams');

    my $disk_allocation = $self->_disk_allocation;
    Carp::confess('No disk allocation!') if not $disk_allocation;

    my $absolute_path = $disk_allocation->absolute_path;
    Carp::confess('No abolute path for disk allocation: '.$disk_allocation->id) if not $absolute_path;

    my @subdirs = grep { -d } glob($absolute_path.'/*');
    if ( not @subdirs ) {
        $self->error_message('No sub directories in disk allocation absolute path after busting bams: '.$absolute_path);
        return;
    }

    my $sample = $self->_sample;
    my $library_cnt = $sample->libraries;
    my %libraries_and_busted_bams;
    for my $subdir ( @subdirs ) {
        $library_cnt++;
        my $library = Genome::Library->create(
            sample_id => $sample->id,
            library_name => $sample->name.'-extlib'.$library_cnt,
        );
        if ( not $library ) {
            $self->error_message("Cannot create library #$library_cnt for sample: ".$sample->name);
            return;
        }
        # check, commit
        my @bams = glob($absolute_path.'/'.$subdir.'/*.bam');
        if ( not @bams ) {
            $self->error_message('No bams in sub directory: '.$subdir);
            return;
        }
        $libraries_and_busted_bams{ $library->id } = [ map { $absolute_path.'/'.$subdir.'/'.$_ } @bams ];
    }

    $self->status_message('Commit libraries');
    if ( not UR::Context->commit ) {
        $self->error_message('Cannot commit libraries');
        return;
    }
    $self->status_message('Commit libraries OK');

    $self->status_message('Got '.scalar(values %libraries_and_busted_bams).' busted bams');

    return %libraries_and_busted_bams;
}

sub _launch_jobs_and_monitor {
    my ($self, %libraries_and_busted_bams)  = @_;

    my $instrument_data = $self->instrument_data;
    my $sample = $self->_sample;

    my %import_params = (
        '--sample' => $sample->id,
        '--import-source-name' => $instrument_data->import_source_name,
        '--import-format' => 'bam',
        '--description' => 'Busted bam from instrument data '.$instrument_data->id,
        '--species-name' => $instrument_data->species_name,
    );
    my $import_params_string = join(' ', %import_params);

    my %jobs;
    for my $library_id ( sort keys %libraries_and_busted_bams ) {
        for my $bam ( @{$libraries_and_busted_bams{$library_id}} ) {
            my $import_cmd = "genome instrument-data import bam --original-data-path $bam --library $library_id $import_params_string";
            my $bsub;
    
            $bsub = sprintf("bsub -g %s -u %s $import_cmd", $self->lsf_job_group_name, $ENV{'USER'} . '@genome.wustl.edu');
            print $bsub, "\n";
            system($bsub);
        }
    }

    # TODO MONITOR

    return 1;
}

1;

=pod

=head1 Disclaimer

Copyright (C) 2005 - 2010 Genome Center at Washington University in St. Louis

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$

