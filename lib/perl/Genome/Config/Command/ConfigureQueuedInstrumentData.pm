package Genome::Config::Command::ConfigureQueuedInstrumentData;

use strict;
use warnings;

use Genome;
use Genome::Sys::LockProxy qw();

use Lingua::EN::Inflect;
use List::MoreUtils qw();

class Genome::Config::Command::ConfigureQueuedInstrumentData {
    is => 'Command::V2',
    has => [
        instrument_data => {
            is          => 'Genome::InstrumentData',
            is_many     => 1,
            is_optional => 1,
            require_user_verify => 0,
            doc         => '[Re]process these instrument data.',
        },
        limit => {
            is => 'Integer',
            default => 50,
            doc => 'Maximum number of instrument data analysis project pairs to process in a single run'
        },
    ],
};

sub help_brief {
    return 'Assign instrument data with an analysis project to models';
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
EOS
}

sub execute {
    my $self = shift;

    $self->_lock();

    my @instrument_data_analysis_project_pairs = $self->_get_items_to_process();
    $self->status_message(sprintf('Found %s to process.',
            Lingua::EN::Inflect::NO('item', scalar(@instrument_data_analysis_project_pairs))));

    for my $current_pair (@instrument_data_analysis_project_pairs) {
        $self->status_message('Working on %s', $current_pair->__display_name__);

        if(my $skip_reason = $self->should_skip($current_pair)) {
            $self->_mark_pair_as_skipped($current_pair, $skip_reason);
            next;
        }

        my $current_inst_data = $current_pair->instrument_data;
        my $analysis_project = $current_pair->analysis_project;
        if (my $msg = $self->should_wait($current_inst_data, $analysis_project)) {
            $self->status_message($msg);
            next;
        }
        eval {
            $analysis_project->get_configuration_profile->process_models_for_instrument_data($current_inst_data);
        };

        my $error = $@;
        $self->error_message($error) if $error;
        $self->_mark_sync_status($current_pair, $error);
    }

    $self->_update_models_for_associated_projects(@instrument_data_analysis_project_pairs);

    return 1;
}

sub _mark_sync_status {
    my ($self, $current_pair, $error) = @_;

    if ($error) {
        $self->_mark_pair_as_failed($current_pair, $error);
    } else {
        $self->_mark_pair_as_processed($current_pair);
    }

    return 1;
}

sub should_skip {
    my ($self, $anp_instdata_bridge) = @_;

    my $inst_data = $anp_instdata_bridge->instrument_data;
    return sprintf('Instrument data (%s) has ignored flag set, skipping!', $inst_data->id) if $inst_data->ignored;
    return;
}

sub should_wait {
    my ($self, $inst_data, $analysis_project) = @_;

    if ($analysis_project->status eq 'Hold') {
        return sprintf("Analysis Project (%s) is set to status 'Hold', skipping!", $analysis_project->id);
    }
    if ($analysis_project->status eq 'Pending') {
        return sprintf("Analysis Project (%s) is still 'Pending', skipping!", $analysis_project->id);
    }
    return;
}

sub _mark_pair_as_skipped {
    my ($self, $current_pair, $skip_reason) = @_;

    $current_pair->fail_count(0);
    $current_pair->reason($skip_reason);
    $current_pair->status('skipped');

    $self->debug_message('Marking pair as skipped.  Reason: %s', $skip_reason);

    return 1;
}

sub _mark_pair_as_processed {
    my ($self, $current_pair) = @_;

    $current_pair->fail_count(0);
    $current_pair->reason(undef);
    $current_pair->status('processed');

    $self->debug_message('Marking pair as processed.');

    return 1;
}

sub _mark_pair_as_failed {
    my ($self, $current_pair, $error_message) = @_;

    my $previous_count = $current_pair->fail_count;
    $current_pair->fail_count(++$previous_count);
    $current_pair->status('failed');
    $current_pair->reason($error_message);

    $self->status_message('Marking pair as failed.  (Total failures on this pair: %s)  Reason: %s', $current_pair->fail_count, $current_pair->reason);

    return 1;
}

sub _get_items_to_process {
    my $self = shift;

    if ($self->instrument_data) {
        return Genome::Config::AnalysisProject::InstrumentDataBridge->get(
            instrument_data_id => [map { $_->id } $self->instrument_data],
            'analysis_project.status' => 'In Progress',
            -hint => ['analysis_project', 'instrument_data', 'instrument_data.sample']
        );
    } else {
        my @bridges;
        my $it = $self->_ordered_sample_id_iterator;
        while(my $sample_id = $it->()) {
            push @bridges, Genome::Config::AnalysisProject::InstrumentDataBridge->get(
                'instrument_data.sample_id' => $sample_id,
                status => ['new', 'failed'],
                'analysis_project.status' => 'In Progress',
                -hint => ['analysis_project', 'instrument_data', 'instrument_data.sample'],
            );
            last if(@bridges >= $self->limit);
        }
        return @bridges;
    }
}

sub _ordered_sample_id_iterator {
    my $self = shift;

    my $it = Genome::Config::AnalysisProject::InstrumentDataBridge->create_iterator(
        status => ['new', 'failed'],
        'analysis_project.status' => 'In Progress',
        -order => ['fail_count'],
        -hint => ['analysis_project', 'instrument_data', 'instrument_data.sample'],
    );

    my %subjects_seen;
    return sub {
        while(my $bridge = $it->next) {
            my $subject_id = $bridge->instrument_data->sample_id;
            next if $subjects_seen{$subject_id}++;
            return $subject_id;
        }
        return;
    };
}

sub _update_models_for_associated_projects {
    my $self = shift;
    my @instrument_data_analysis_project_pairs = @_;

    my @instrument_data_ids = map { $_->instrument_data->id } @instrument_data_analysis_project_pairs;

    my @projects = Genome::Project->get('parts.entity_id' => \@instrument_data_ids);

    if(@projects) {
        my $update_cmd = Genome::Project::Command::Update::Models->create(
            projects => \@projects
        );

        die $self->error_message('Failed to update models for project associated with %s', join(',', @instrument_data_ids))
            unless $update_cmd->execute();
    }

    return 1;
}

sub _lock {
    unless ($ENV{UR_DBI_NO_COMMIT}) {
        my $lock = Genome::Sys::LockProxy->new(
            resource => 'genome_config_command_configure-queued-instrument-data/lock',
            scope => 'site',
        )->lock(max_try => 1);

        die('Unable to acquire the lock! Is ConfigureQueuedInstrumentData already running or did it exit uncleanly?')
            unless $lock;

        Genome::Sys::CommitAction->create(
            on_commit => sub {
                $lock->unlock();
            }
        );
    }
}

1;
