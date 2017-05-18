use strict;
use warnings;

package Genome::Model::Command::Admin::PurgeSoftwareResultsFromAnalysisProject;

use DateTime::Format::Strptime;
use Scope::Guard;
use Genome;
use Genome::Sys::LockProxy;

class Genome::Model::Command::Admin::PurgeSoftwareResultsFromAnalysisProject {
    is => 'Command::V2',
    has => [
        days_to_retain => {
            is => 'Integer',
            default_value => 30,
            doc => 'Things older than this many days will be purged',
        },
        analysis_projects => {
            is => 'Genome::Config::AnalysisProject',
            is_many => 1,
            shell_args_position => 1,
            require_user_verify => 1,
            doc => 'List of AnalysisProjects to purge',
        },
        report => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Print a report of what would be purged, but do not actually purge any data',
        },
    ],
};

my $sql = qq(
    SELECT DISTINCT ON (candidates.result_id) candidates.result_id, candidates.profile_item_id, candidates.anp_name, candidates.anp_id, candidates.kilobytes_requested
    FROM (
        SELECT sr.id result_id, sr.class_name, bridge.analysis_project_id anp_id, anp.name anp_name, cpi.id profile_item_id, sr.outputs_path, da.id allocation_id, da.kilobytes_requested, da.status
        FROM result.software_result sr
        INNER JOIN disk.allocation da ON da.owner_id = sr.id AND da.owner_class_name = sr.class_name
        INNER JOIN result."user" sru ON sr.id = sru.software_result_id
        INNER JOIN model.build b ON b.build_id = sru.user_id
        INNER JOIN config.analysis_project_model_bridge bridge ON bridge.model_id = b.model_id
        INNER JOIN config.analysis_project anp ON bridge.analysis_project_id = anp.id
        INNER JOIN config.profile_item cpi ON cpi.id = bridge.profile_item_id
        WHERE cpi.status = 'disabled'
        AND anp.id = ?
        AND da.kilobytes_requested > 0
    ) candidates
    WHERE NOT EXISTS (
        SELECT sr.id
        FROM result.software_result sr
        INNER JOIN result."user" sru ON sr.id = sru.software_result_id
        INNER JOIN model.build b ON b.build_id = sru.user_id
        INNER JOIN config.analysis_project_model_bridge bridge ON bridge.model_id = b.model_id
        INNER JOIN config.profile_item cpi ON cpi.id = bridge.profile_item_id
        WHERE cpi.status != 'disabled'
        AND candidates.result_id = sr.id
    ) AND NOT EXISTS (
        SELECT sr.id
        FROM result.software_result sr
        INNER JOIN result."user" sru ON sr.id = sru.software_result_id
        WHERE sru.user_class_name = 'Genome::Sys::User'
        AND candidates.result_id = sr.id
    )
);

sub execute {
    my $self = shift;

    my $sth = $self->_prepare_sql_statement();
    foreach my $anp ( $self->analysis_projects ) {
        $self->purge_one_analysis_project($anp, $sth);
    }

    1;
}

sub _prepare_sql_statement {
    my $self = shift;
    my $dbh = Genome::DataSource::GMSchema->get_default_handle();
    die $self->error_message('Cannot get default database handle') unless $dbh;

    my $sth = $dbh->prepare($sql);
    die $self->error_message("preparing SQL failed : $DBI::errstr") unless $sth;
    return $sth;
}

sub analysis_project_is_old_enough_to_purge {
    my($self, $anp) = @_;

    my $strp = DateTime::Format::Strptime->new(
        pattern   => UR::Context->date_template(),
    );

    my $now_str = UR::Context->now();
    my $updated_at_str = $anp->updated_at;

    my $dt1 = $strp->parse_datetime($updated_at_str);
    my $dt2 = $strp->parse_datetime($now_str);

    my $anp_updated_duration = $dt1->delta_days($dt2);

    my $duration_to_retain = DateTime::Duration->new(
        days        => $self->days_to_retain,
    );

    return ( DateTime::Duration->compare( $anp_updated_duration, $duration_to_retain ) == -1 );
}


sub _get_lock_for_analysis_project {
    my($self, $anp) = @_;

    my $unlocker;

    if ($self->report) {
        $unlocker = sub { $self->_report_obj->print() };

    } elsif ($ENV{UR_DBI_NO_COMMIT}) {
        $unlocker = sub {};

    } else {
        my $lock = Genome::Sys::LockProxy->new(
            resource => 'genome model admin purge-software-results-from-analysis-project' . $anp->id,
            scope => 'site',
        )->lock(max_try => 1);

        die $self->error_message('Cannot get lock for analysis project '.$anp->id) unless $lock;

        $unlocker = sub {
            $lock->unlock();
        };
    }
    return $unlocker;
}

sub purge_one_analysis_project {
    my($self, $anp, $sth) = @_;

    if ($anp->is_cle) {
        die $self->error_message('Failed to process CLE analysis project: '. $anp->id);
    }

    if ( $self->analysis_project_is_old_enough_to_purge($anp) ) {
        $self->warning_message('Analysis project \''. $anp->id .'\' was updated at '. $anp->updated_at .' which is less than the '. $self->days_to_retain .' days to retain disabled results');
        return;
    }

    my $unlocker = $self->_get_lock_for_analysis_project($anp);

    my $total_kb_purged = 0;
    my $software_result_count = 0;

    my $unlock_and_print_report = Scope::Guard->new(sub {
        $unlocker->();
        my $format = $self->report
                      ? '%s in %d software results for %s'
                      : 'Removed %s from %d software results for %s';
        $self->status_message(
            $format,
            $self->_report_obj->_format_disk_size($total_kb_purged),
            $software_result_count,
            $anp->__display_name__
        );
    });

    $sth->execute($anp->id);

    while (my $data = $sth->fetchrow_hashref()) {
        my $sr = Genome::SoftwareResult->get($data->{result_id});
        die $self->error_message('Cannot get software result '. $data->{result_id}) unless $sr;

        my $reason = 'Expunge software result uniquely used by model from disabled config item ('. $data->{profile_item_id} .') for analysis project \''. $data->{anp_name} .'\' ('. $data->{anp_id} .')';

        $software_result_count++;
        $total_kb_purged += $data->{kilobytes_requested};

        if ($self->report) {
            $self->_add_to_report($sr, $data->{kilobytes_requested});
        } else {
            $self->_do_expunge($sr, $reason);
        }
    }
}

my $report_obj;
sub _report_obj {
    my $self = shift;
    $report_obj ||= Genome::Model::Command::Admin::PurgeSoftwareResultsFromAnalysisProject::PurgeReport->new($self);
}

sub _add_to_report {
    my($self, $software_result, $kb_requested) = @_;

    $self->_report_obj->add($software_result, $kb_requested);
}

sub _do_expunge {
    my($self, $sr, $reason) = @_;

    $self->status_message('Removing software result %s',
                            $sr->id);
    $sr->expunge($reason);
    UR::Context->commit() || die "commit() failed while expunging software result ".$sr->id;
}

package Genome::Model::Command::Admin::PurgeSoftwareResultsFromAnalysisProject::PurgeReport;
use List::Util qw(reduce);

sub new {
    my($class, $purger) = @_;
    my $self = {
        purger => $purger,
        data => {},
    };
    return bless($self, $class);
}

sub purger { shift->{purger} }
sub data { shift->{data} }

sub add {
    my($self, $software_result, $kb_requested) = @_;
    $self->data->{ $software_result->class }->{ $software_result->id } ||= $kb_requested;
}

sub print {
    my $self = shift;

    foreach my $type ( sort keys %{ $self->data } ) {
        my $this_class_data = $self->data->{$type};
        my $total_size = reduce { $a + $b } values %$this_class_data;
        my $count = scalar keys %$this_class_data;
        $self->purger->status_message('%s: %s in %s software results',
                                        $type,
                                        $self->_format_disk_size($total_size),
                                        $count);
        foreach my $id ( keys %$this_class_data ) {
            $self->purger->status_message("\t%s => %s",
                                            $id,
                                            $self->_format_disk_size($this_class_data->{$id}));
        }
    }
    $self->{data} = {};
}

sub _format_disk_size {
    my($self, $kb) = @_;

    if ($kb < 1024) {
        "$kb KB";
    } elsif ($kb < 1024 * 1024) {
        sprintf('%0.3f MB', $kb / 1024);
    } else {
        sprintf('%0.3f GB', $kb / 1024 / 1024);
    }
}

1;
