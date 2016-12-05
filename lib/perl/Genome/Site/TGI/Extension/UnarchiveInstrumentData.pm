package Genome::Config::AnalysisProject::Command::UnarchiveInstrumentData;

use strict;
use warnings;

use Genome;
use Genome::Config::AnalysisProject::Command::UnarchiveInstrumentData; #load real module first

use File::Basename qw();
use File::Spec qw();

use Sub::Install;

Sub::Install::reinstall_sub ({
    as => 'unarchive_additional_data',
    code => \&unarchive_lims_data,
});

__PACKAGE__->__meta__->property('volume')->is_optional(0);

sub unarchive_lims_data {
    my $self = shift;
    my $instrument_data = shift;

    my $lims_unarchive_script_path = $self->resolve_lims_unarchive_script_path();

    if ($instrument_data->class ne 'Genome::InstrumentData::Solexa') { return 1; }

    my $bam_path = $instrument_data->bam_path;
    if (defined($bam_path)) {
        if (-e $bam_path) {
            return 1;
        }
    } else {
        if (defined($instrument_data->archive_path)) {
            $self->warning_message('Skipping old data with archive_path: '. $instrument_data->__display_name__);
            return 1;
        }
    }

    $self->status_message('Working on missing BAM %s for instrument data %s', $bam_path, $instrument_data->id);

    my $dv = $self->volume;
    my ($bam_filename, $bam_dirname) = File::Basename::fileparse($bam_path);

    my ($unarchived_bam_path) = glob(File::Spec->join($dv->mount_path,'*','csf_*',$bam_filename));
    if (!defined($unarchived_bam_path)) {
        # Try one more tine allowing for an additional directory in the path, ie. "condensed"
        ($unarchived_bam_path) = glob(File::Spec->join($dv->mount_path,'*','*','csf_*',$bam_filename));
    }
    unless ($unarchived_bam_path) {
        Genome::Sys->shellcmd(
            cmd => [$lims_unarchive_script_path,$instrument_data->id,$dv->mount_path],
        );
    } else {
        my $allocation = $self->get_or_create_allocation($instrument_data);

        my ($unarchived_bam_filename,$unarchived_bam_dirname) = File::Basename::fileparse($unarchived_bam_path);
        $unarchived_bam_dirname =~ s/\/$//;
        my @unarchived_bam_dirnames = File::Spec->splitdir($unarchived_bam_dirname);
        my $allocation_bam_dirname = File::Spec->join($allocation->absolute_path,$unarchived_bam_dirnames[-1]);
        if (-e $allocation_bam_dirname) {
            $self->warning_message('Already found allocated path, skipping: %s', $allocation_bam_dirname);
            $self->warning_message('Please remove the unarchived path: %s', $unarchived_bam_dirname);
            return;
        }

        unless ($self->validate_bam_is_complete($unarchived_bam_path,$instrument_data)) {
            $self->status_message('It appears the unarchiving is still in progress for: %s', $unarchived_bam_path);
            return;
        }

        $self->debug_message('Rsync contents of %s to %s ', $unarchived_bam_dirname, $allocation_bam_dirname);
        unless (
            Genome::Sys->rsync_directory(
                source_directory => $unarchived_bam_dirname,
                target_directory => $allocation_bam_dirname,
            )
          ) {
            $self->fatal_message('Failed to rsync from %s to %s', $unarchived_bam_dirname, $allocation_bam_dirname);
        }
        unless ($ENV{UR_DBI_NO_COMMIT}) {
            Genome::Sys::CommitAction->create(
                on_commit => sub {
                    unless(Genome::Sys->remove_directory_tree($unarchived_bam_dirname)) {
                        $self->warning_message('Failed to remove directory: %s', $unarchived_bam_dirname);
                    }
                }
            );
        }
        $allocation->reallocate();
        UR::Context->commit();

        my $allocated_bam_path = File::Spec->join($allocation_bam_dirname,$unarchived_bam_filename);
        if (-e $allocated_bam_path) {
            $instrument_data->bam_path($allocated_bam_path);
            UR::Context->commit();
        } else {
            $self->fatal_message('Failed to find the new allocated BAM path: %s', $allocated_bam_path);
        }
    }
}

sub resolve_lims_unarchive_script_path {
    my $self = shift;

    my $abs_path = File::Spec->rel2abs(__FILE__);
    my ($filename,$dirname) = File::Basename::fileparse($abs_path);
    my $lims_unarchive_script_path = File::Spec->join($dirname,'generate_lims_disk_unarchive.pl');
    unless (-e $lims_unarchive_script_path) {
        $self->fatal_message('Failed to resolve the path to the LIMS script required for unarchiving.');
    }
    return $lims_unarchive_script_path;
}

sub get_or_create_allocation {
    my $self = shift;
    my $instrument_data = shift;
    my %params = (
        disk_group_name => Genome::Config::get('disk_group_alignments'),
        allocation_path => File::Spec->join('instrument_data',$instrument_data->id),
        kilobytes_requested => $instrument_data->calculate_alignment_estimated_kb_usage,
        owner_class_name => $instrument_data->class,
        owner_id => $instrument_data->id,
    );
    my ($allocation) = Genome::Disk::Allocation->get(allocation_path => $params{allocation_path});
    if ($allocation) {
        return $allocation;
    } else {
        my $create_cmd = Genome::Disk::Command::Allocation::Create->create(%params);
        unless ($create_cmd->execute) {
            $self->fatal_message('Could not create allocation for instrument data: %s', $instrument_data->id);
        }
        ($allocation) = Genome::Disk::Allocation->get(allocation_path => $params{allocation_path});
        unless ($allocation) {
            $self->fatal_message('Could not find allocation for instrument data: %s', $instrument_data->id);
        }
    }
    return $allocation;
}

sub validate_bam_is_complete {
    my $self = shift;
    my $bam_path = shift;
    my $instrument_data = shift;

    if (-e $bam_path .'.md5') {
        $self->status_message('Calculating md5 for BAM: '. $bam_path);
        my $calculated_md5 = Genome::Sys->md5sum($bam_path);
        my $md5_fh = Genome::Sys->open_file_for_reading($bam_path .'.md5');
        my ($actual_md5) = split(/\s+/, $md5_fh->getline);
        $md5_fh->close();
        if ($calculated_md5 ne $actual_md5) {
            $self->debug_message('The BAM md5s do not match. Unarchiving may still be in progress.');
            return;
        }
        $self->debug_message('md5 for BAM is a match: %s', $bam_path);
    } else {
        $self->debug_message('Running flagstat to validate BAM: %s', $bam_path);
        my $temp_flagstat = Genome::Sys->create_temp_file_path();

        my $flagstat_cmd = Genome::Model::Tools::Sam::Flagstat->create(
            bam_file => $bam_path,
            output_file => $temp_flagstat
        );
        unless ($flagstat_cmd->execute) {
            $self->fatal_message('Failed to run flagstat on BAM: %s', $bam_path);
        }

        my $calculated_flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($temp_flagstat);
        unless ($calculated_flagstat) {
            $self->fatal_message('Fail to parse flagstat data file: %s', $temp_flagstat);
        }
        if (exists $calculated_flagstat->{errors}) {
            my @errors = @{$calculated_flagstat->{errors}};
            for my $error (@errors) {
                if ($error =~ 'Truncated file') {
                    $self->debug_message('BAM file appears to be truncated: %s', $bam_path);
                    return;
                } else {
                    $self->warning_message('Error messages from flagstat: %s', $error);
                    return;
                }
            }
        }
        if (!exists($calculated_flagstat->{total_reads})) {
            $self->fatal_message('Failed to parse total_reads in BAM file: %s', $bam_path);
        }
        if (-e $bam_path .'.flagstat') {
            my $actual_flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($bam_path .'.flagstat');
            unless ($actual_flagstat) {
                $self->fatal_message('Fail to parse flagstat data file: %s.flagstat', $bam_path);
            }
            if ($calculated_flagstat->{total_reads} ne $actual_flagstat->{total_reads}) {
                $self->fatal_message('BAM file total_reads does not match expected total reads: %s', $bam_path);
            }
        } else {
            if ($calculated_flagstat->{total_reads} ne $instrument_data->read_count) {
                $self->fatal_message('BAM file total_reads does not match expected instrument data read count: %s', $bam_path);
            }
        }
        $self->debug_message('The read count in BAM is confirmed: %s'. $bam_path);
    }
    $self->status_message('Confirmed BAM file is complete: %s', $bam_path);
    return 1;
}

1;
