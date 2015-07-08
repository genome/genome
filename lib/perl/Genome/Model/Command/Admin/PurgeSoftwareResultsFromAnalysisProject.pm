use strict;
use warnings;

use Data::Dumper;
use DateTime::Format::Strptime;
use Genome;

my $days_to_retain = 30;

my $sql = qq(
    SELECT * FROM (
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

my $anp_id = $ARGV[0];

my $anp = Genome::Config::AnalysisProject->get($anp_id);
die unless $anp;

if ($anp->is_cle) {
    die('Failed to process CLE analysis project: '. $anp->id);
}

my $strp = DateTime::Format::Strptime->new(
    pattern   => UR::Context->date_template(),
);

my $now_str = UR::Context->now();
my $updated_at_str = $anp->updated_at;

my $dt1 = $strp->parse_datetime($updated_at_str);
my $dt2 = $strp->parse_datetime($now_str);

my $anp_updated_duration = $dt1->delta_days($dt2);

my $duration_to_retain = DateTime::Duration->new(
    days        => $days_to_retain,
);

if ( DateTime::Duration->compare( $anp_updated_duration, $duration_to_retain ) == -1 ) {
    warn('Analysis project \''. $anp_id .'\' was updated at '. $updated_at_str .' which is less than the '. $days_to_retain .' days to retain disabled results');
    exit;
}

my $dbh = Genome::DataSource::GMSchema->get_default_handle();
die unless $dbh;

my $sth = $dbh->prepare($sql);
die unless $sth;

$sth->execute($anp_id);
while (my $data = $sth->fetchrow_hashref()) {
    my $sr = Genome::SoftwareResult->get($data->{result_id});
    die unless $sr;

    my $reason = 'Expunge software result uniquely used by model from disabled config item ('. $data->{profile_item_id} .') for analysis project \''. $data->{anp_name} .'\' ('. $data->{anp_id} .')';
    print $reason ."\n";
    $sr->expunge($reason);
    UR::Context->commit();
}

exit;
