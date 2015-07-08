use strict;
use warnings;

package Genome::Model::Command::Admin::PurgeSoftwareResultsFromAnalysisProject;

use DateTime::Format::Strptime;
use Scope::Guard;
use Genome;
use Genome::Sys::LockProxy;

use constant KB_IN_ONE_GB => 1048576;

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
        dry_run => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Do not actually purge anything',
        },
    ],
};

my $sql = qq(
    SELECT candidates.result_id, candidates.profile_item_id, candidates.anp_name, candidates.anp_id, candidates.kilobytes_requested
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
        WHERE cpi.status = 'active'
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

    foreach my $anp ( $self->analysis_projects ) {
        $self->purge_one_analysis_project($anp);
    }
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
    if ($ENV{UR_DBI_NO_COMMIT}) {
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
    my($self, $anp) = @_;

    if ($anp->is_cle) {
        die $self->error_message('Failed to process CLE analysis project: '. $anp->id);
    }

    if ( $self->analysis_project_is_old_enough_to_purge($anp) ) {
        $self->warning_message('Analysis project \''. $anp->id .'\' was updated at '. $anp->updated_at .' which is less than the '. $self->days_to_retain .' days to retain disabled results');
        return;
    }

    my $dbh = Genome::DataSource::GMSchema->get_default_handle();
    die $self->error_message('Cannot get default database handle') unless $dbh;

    my $unlocker = $self->_get_lock_for_analysis_project($anp);

    my $total_kb_purged = 0;
    my $software_result_count = 0;

    my $unlock_and_print_report = Scope::Guard->new(sub {
        $unlocker->();
        $self->status_message('Removed %d GB from %d software results',
                                int($total_kb_purged / KB_IN_ONE_GB),
                                $software_result_count);
    });

    my $sth = $dbh->prepare($sql);
    die $self->error_message("preparing SQL failed : $DBI::errstr") unless $sth;

    $sth->execute($anp->id);

    while (my $data = $sth->fetchrow_hashref()) {
        my $sr = Genome::SoftwareResult->get($data->{result_id});
        die $self->error_message('Cannot get software result '. $data->{result_id}) unless $sr;

        my $reason = 'Expunge software result uniquely used by model from disabled config item ('. $data->{profile_item_id} .') for analysis project \''. $data->{anp_name} .'\' ('. $data->{anp_id} .')';

        $software_result_count++;
        $total_kb_purged += $data->{kilobytes_requested};

        if ($self->dry_run) {
            $self->warning_message('Dry run, not removing software result '.$sr->id);
        } else {
            $self->status_message($reason);
            $sr->expunge($reason);
            UR::Context->commit();
        }
    }
}
