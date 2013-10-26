package Genome::InstrumentData::Command::Summary1;

#REVIEW fdu 11/20/2009
#What does Summary1 stand for ?

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::InstrumentData::Command::Summary1 {
    is => 'Command',
    has => [
        sql => { 
            is_constant => 1,
            value => q{ 
                select project_name,sample_name,library_name,run_date,flow_cell_id,lanes,total_clusters,
                    round(total_bases/1000) kb, 
                    filter_clusters, 
                    pf_clusters, 
                    (case when filter_clusters_r1 > 0 then round(aligned_r1_clusters/filter_clusters_r1,2) else 0 end) pct_align_r1, 
                    (case when filter_clusters_r2 > 0 then round(aligned_r2_clusters/filter_clusters_r2,2) else 0 end) pct_align_r2, 
                    (case when filter_clusters_r1 > 1 then round(errors_r1/filter_clusters_r1,2) else 0 end) error_rate_r1,                       
                    (case when filter_clusters_r2 > 1 then round(errors_r2/filter_clusters_r2,2) else 0 end) error_rate_r2                       
                from (
                    select research_project project_name,sample_name,library_name, flow_cell_id, 
                        count(distinct lane) lanes,
                        sum(clusters*read_length) total_bases, 
                        sum(clusters) total_clusters, 
                        sum(filt_clusters) filter_clusters, 
                        (case when sum(clusters) > 0 then round(sum(filt_clusters)/sum(clusters),2) else 0 end) pf_clusters, 
                        sum(case when run_type != 'Paired End Read 2' then filt_aligned_clusters_pct*filt_clusters else 0 end) aligned_r1_clusters, 
                        sum(case when run_type = 'Paired End Read 2' then filt_aligned_clusters_pct*filt_clusters else 0 end) aligned_r2_clusters, 
                        sum(case when run_type != 'Paired End Read 2' then filt_error_rate_avg*filt_clusters else 0 end) errors_r1,  
                        sum(case when run_type = 'Paired End Read 2' then filt_error_rate_avg*filt_clusters else 0 end) errors_r2, 
                        '20' || substr(run_name,0,2) || '-' || substr(run_name,3,2) || '-' || substr(run_name,5,2) run_date, 
                        sum(case when run_type != 'Paired End Read 2' then filt_clusters else 0 end) filter_clusters_r1, 
                        sum(case when run_type = 'Paired End Read 2' then filt_clusters else 0 end) filter_clusters_r2   
                    from solexa_lane_summary
                    where sample_name like ? 
                    group by research_project,sample_name,library_name,flow_cell_id,run_name
                ) x  
                /*
                left join (
                    setup@oltp s 
                    join setup_project@oltp p on p.setup_project_id = s.setup_id
                    join contact@oltp c on c.ca_id = p.ext_con_id
                ) on s.setup_name = x.project_name
                */
                where run_date like ?
            },
        },
    ],
    has_optional => [
        sample_pattern  => { 
            is => 'Text', 
            default_value => '%',
            shell_args_position => 1, 
            doc => 'all or part of a sample name, using % as a wildcard', 
        },
        date_of_run => { 
            is => 'Text', 
            default_value => '%',
            doc => 'the run date in YYMMDD format, using % as a wildcard (defaults to all dates)', 
        },
    ],
    doc => 'summarize illumina data for a sample by library/flow cell',
};

sub help_synopsis {
    return <<EOS
# everything with a given patient ID in the name
genome instrument-data summary1 %6888%

# everything for a given year
genome instrument-data summary1 -d 2009-%%-%%
EOS
}

sub execute {
    my $self = shift;

    my $sql = $self->sql;
    my $sample_pattern = $self->sample_pattern;
    my $date = $self->date_of_run;

    if ($date ne '%') {
        $sql .= '
                order by run_date,sample_name,library_name,flow_cell_id';
    }
    else {
        $sql .= '
                order by sample_name,library_name,run_date,flow_cell_id';
    }
    
    $sample_pattern = "\%$sample_pattern\%" unless $sample_pattern eq '%';

    my $dbh = Genome::DataSource::Dwrac->get_default_handle();
    my $sth = $dbh->prepare($sql) or die "Failed to connect to the database!";
    $sth->execute($sample_pattern,$date) or die "Error querying the database!";
    UR::DBI::Report->print_formatted(sth => $sth);
    return 1;
}

1;

