package Genome::Site::TGI::Synchronize::Classes::IndexIllumina;

use strict;
use warnings;

use Genome;

=comment
# FLOW_CELL_ILLUMINA
CREATION_EVENT_ID NUMBER   (15)                    {null} NOT NULL NOT_SYNCED
FLOW_CELL_ID      VARCHAR2 (15)                    {null} NOT NULL NOT_SYNCED
GROUP_NAME        VARCHAR2 (64)                    {null} {null}   NOT_SYNCED
MACHINE_NAME      VARCHAR2 (64)                    {null} NOT NULL NOT_SYNCED
RUN_NAME          VARCHAR2 (64)                    {null} NOT NULL ok
RUN_TYPE          VARCHAR2 (25)                    {null} NOT NULL NOT_SYNCED [used in join, but value stored is calculated from reads]
TEAM_NAME         VARCHAR2 (64)                    {null} {null}   NOT_SYNCED

# INDEX_ILLUMINA
ANALYSIS_ID                    NUMBER   (15)                    {null} NOT NULL {null} ok [id]
ANALYSIS_NAME                  VARCHAR2 (15)                    {null} NOT NULL {null} NOT_SYNCED
ANALYSIS_SOFTWARE_VERSION      VARCHAR2 (32)                    {null} NOT NULL {null} ok
AVG_FILT_CLUSTERS_ACROSS_TILES NUMBER   (10)                    {null} {null}   {null} NOT_SYNCED
CLUSTERS                       NUMBER   (10)                    {null} NOT NULL {null} NOT_SYNCED
CLUSTERS_AVG                   NUMBER   (12,2)                  {null} NOT NULL {null} NOT_SYNCED
CLUSTERS_STDEV                 NUMBER   (12,2)                  {null} NOT NULL {null} NOT_SYNCED
CREATION_EVENT_ID              NUMBER   (15)                    {null} NOT NULL {null} NOT_SYNCED
EXPERIMENT_PURPOSE             VARCHAR2 (32)                    {null} {null}   {null} NOT_SYNCED
FILT_CLUSTERS                  NUMBER   (17,2)                  {null} NOT NULL {null} ok [clusters]
FILT_CLUSTERS_AVG              NUMBER   (14,2)                  {null} NOT NULL {null} NOT_SYNCED
FILT_CLUSTERS_STDEV            NUMBER   (14,2)                  {null} NOT NULL {null} NOT_SYNCED
FLOW_CELL_ID                   VARCHAR2 (15)                    {null} NOT NULL {null} ok
GERALD_DIRECTORY               VARCHAR2 (255)                   {null} NOT NULL {null} ok
INDEX_SEQUENCE                 VARCHAR2 (64)                    {null} {null}   {null} ok
INTEN_AVG                      NUMBER   (15)                    {null} {null}   {null} NOT_SYNCED
INTEN_STDEV                    NUMBER   (15)                    {null} {null}   {null} NOT_SYNCED
LANE                           NUMBER   (3)                     {null} NOT NULL {null} ok
LIBRARY_ID                     NUMBER   (20)                    {null} NOT NULL {null} ok
LIBRARY_STRATEGY               VARCHAR2 (30)                    {null} {null}   {null} NOT_SYNCED
MEDIAN_INSERT_SIZE             NUMBER   (7)                     {null} {null}   {null} ok
RESEARCH_PROJECT               VARCHAR2 (64)                    {null} {null}   {null} NOT_SYNCED
SD_ABOVE_INSERT_SIZE           NUMBER   (7)                     {null} {null}   {null} ok
SD_BELOW_INSERT_SIZE           NUMBER   (7)                     {null} {null}   {null} ok
SEQ_ID                         NUMBER   (15)                    {null} NOT NULL {null} NOT_SYNCED
SS_ID                          NUMBER   (15)                    {null} {null}   {null} NOT_SYNCED
TAGGED_DNA_ID                  NUMBER   (20)                    {null} {null}   {null} NOT_SYNCED
TARGET_REGION_SET_NAME         VARCHAR2 (512)                   {null} {null}   {null} ok

# READ_ILLUMINA
II_SEQ_ID                   NUMBER   (15)                    {null} NOT NULL NOT_SYNCED
READ_NUMBER                 NUMBER   (1)                     {null} NOT NULL NOT_SYNCED
READ_LENGTH                 NUMBER   (3)                     {null} NOT NULL ok [from r2]
INTEN_AFTER_20_CYCLES_PCT   NUMBER   (9,2)                   {null} {null}   NOT_SYNCED
INTEN_AFTER_20_CYCLES_STDEV NUMBER   (9,2)                   {null} {null}   NOT_SYNCED
FILT_ALIGNED_CLUSTERS_PCT   NUMBER   (9,2)                   {null} {null}   ok [fwd_filt_aligned_clusters_pct from r1; rev_filt_aligned_clusters_pct from r2]
FILT_ALIGNED_CLUSTERS_STDEV NUMBER   (9,2)                   {null} {null}   NOT_SYNCED
FILT_ERROR_RATE_AVG         NUMBER   (9,2)                   {null} {null}   ok [from r2; fwd_filt_error_rate_avg from r1; rev_filt_error_rate_avg from r2]
FILT_ERROR_RATE_STDEV       NUMBER   (9,2)                   {null} {null}   NOT_SYNCED
FIRST_CYCLE_INTEN_AVG       NUMBER   (17,2)                  {null} NOT NULL NOT_SYNCED
FIRST_CYCLE_INTEN_STDEV     NUMBER   (17,2)                  {null} NOT NULL NOT_SYNCED
PHASING_PCT                 NUMBER   (9,2)                   {null} NOT NULL NOT_SYNCED
PREPHASING_PCT              NUMBER   (9,2)                   {null} NOT NULL NOT_SYNCED
KILOBASES_READ              NUMBER   (12)                    {null} NOT NULL ok [fwd_kilobases_read from r2 if r1.seq_id is not null; rev_kilobases_read from r2 if r1.seq_id is not null] 
RECIPE_READ_NUMBER          NUMBER   (1)                     {null} NOT NULL NOT_SYNCED
ORIENTATION                 VARCHAR2 (50)                    {null} NOT NULL NOT_SYNCED
SLS_SEQ_ID                  NUMBER   (15)                    {null} {null}   ok [fwd_seq_id from r1; rev_seq_id from r2 if r1.sed_id is not null]
SEQ_ID                      NUMBER   (15)                    {null} NOT NULL NOT_SYNCED
N_BASE_COUNT                NUMBER   (15)                    {null} {null}   NOT_SYNCED
N_READ_COUNT                NUMBER   (15)                    {null} {null}   NOT_SYNCED
Q2_BASE_COUNT               NUMBER   (15)                    {null} {null}   NOT_SYNCED
Q2_READ_COUNT               NUMBER   (15)                    {null} {null}   NOT_SYNCED
Q20_BASE_COUNT              NUMBER   (15)                    {null} {null}   NOT_SYNCED
BASE_QUALITY_SUM            NUMBER   (15)                    {null} {null}   NOT_SYNCED
Q30_BASE_COUNT              NUMBER   (15)                    {null} {null}   NOT_SYNCED
FIRST_AVG_Q_BELOW_20        VARCHAR2 (7)                     {null} {null}   NOT_SYNCED

=cut

class Genome::Site::TGI::Synchronize::Classes::IndexIllumina {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsInstDataBase',
    table_name => <<'EOS'
        (
            select
                --Index Illumina
                to_char(i.analysis_id) id,
                i.library_id library_id,
                i.index_sequence index_sequence,
                i.flow_cell_id,
                i.lane lane,
                i.target_region_set_name,
                i.gerald_directory,
                i.median_insert_size old_median_insert_size,
                i.sd_above_insert_size old_sd_above_insert_size,
                i.sd_below_insert_size old_sd_below_insert_size,
                i.filt_clusters clusters,
                i.analysis_software_version,
                (case when i.index_sequence is null then to_char(i.lane) else to_char(i.lane) || '-' || i.index_sequence end) subset_name,

                --Constant
                0 is_external,

                --Flow Cell
                fc.run_name run_name,

                --Read Illumina #1 r1 is always reverse when defined [only if the run is paired end]
                (case when r1.seq_id is not null then 'Paired' else 'Standard' end) run_type,

                --Read Illumina #2 r2 is fwd for fragment and rev for paired end
                r2.read_length,
                r2.filt_error_rate_avg old_filt_error_rate_avg,

                --Fwd
                r1.sls_seq_id fwd_seq_id,
                (case when r1.seq_id is not null then 'Paired End Read 1' else null end) fwd_run_type,
                nvl(r1.read_length,-1) fwd_read_length,
                (case when r1.seq_id is not null then r2.kilobases_read else -1 end) fwd_kilobases_read,
                (case when r1.seq_id is not null then i.filt_clusters else null end) fwd_clusters,
                (case when r1.seq_id is not null then i.filt_clusters else null end) fwd_filt_clusters,
                r1.filt_aligned_clusters_pct old_fwd_filt_aligned_clusters_pct,
                r1.filt_error_rate_avg old_fwd_filt_error_rate_avg,

                --Rev
                (case when r1.seq_id is not null then r2.sls_seq_id else null end) rev_seq_id,
                (case when r1.seq_id is not null then 'Paired End Read 2' else null end) rev_run_type,
                (case when r1.seq_id is not null then r2.read_length else -1 end) rev_read_length,
                (case when r2.seq_id is not null then r2.kilobases_read else -1 end) rev_kilobases_read,
                (case when r1.seq_id is not null then i.filt_clusters else null end) rev_clusters,
                (case when r1.seq_id is not null then i.filt_clusters else null end) rev_filt_clusters,
                (case when r1.seq_id is not null then r2.filt_error_rate_avg else null end) old_rev_filt_error_rate_avg,
                (case when r1.seq_id is not null then r2.filt_aligned_clusters_pct else null end) old_rev_filt_aligned_clusters_pct,

                --Misc Paths
                archive2.path archive_path,
                gerald_bam.path bam_path,
                collect_gc_bias.path gc_bias_path,
                fastqc.path fastqc_path,
                '/gscmnt/sata114/info/medseq/adaptor_sequences/solexa_adaptor_pcr_primer'
                    || (case when sam.sample_type = 'rna' then '_SMART' else '' end) adaptor_path

                from index_illumina i
                    join flow_cell_illumina fc on fc.flow_cell_id = i.flow_cell_id
                    join read_illumina r2
                        on i.seq_id = r2.ii_seq_id
                        and (
                            (fc.run_type = 'Paired End' and r2.read_number = 2)
                            or
                            (fc.run_type = 'Fragment' and r2.read_number = 1)
                        )
                    left join seq_fs_path archive2 on archive2.seq_id = i.seq_id
                        and archive2.data_type = 'illumina fastq tgz'
                    left join seq_fs_path gerald_bam on gerald_bam.seq_id = i.seq_id
                        and gerald_bam.data_type = 'gerald bam'
                    left join seq_fs_path collect_gc_bias on collect_gc_bias.seq_id = i.seq_id
                        and collect_gc_bias.data_type = 'collect gc bias'
                    left join seq_fs_path fastqc on fastqc.seq_id = i.seq_id
                        and fastqc.data_type = 'fastqc'
                    left join read_illumina r1
                        on run_type = 'Paired End'
                        and r1.ii_seq_id = i.seq_id
                        and r1.read_number = 1
                    join GSC.library_summary lib on lib.library_id = i.library_id
                    join GSC.organism_sample sam on sam.organism_sample_id = lib.sample_id
        )
        index_illumina
EOS
    ,
    id_by => [
        id => { is => 'Text', },
    ],
    has => [
        library_id                      => { is => 'Number', },
        is_external                     => { },
    ],
    has_optional => [
        flow_cell_id                     => { }, # = short name
        lane                             => { },
        subset_name                      => { },
        index_sequence                   => { },
        run_name                         => { },
        run_type                         => { },
        fastqc_path                      => { },
        read_length                      => { },
        fwd_read_length                  => { },
        rev_read_length                  => { },
        fwd_kilobases_read               => { },
        rev_kilobases_read               => { },
        fwd_run_type                     => { },
        rev_run_type                     => { },
        gerald_directory                 => { },
        old_median_insert_size           => { },
        old_sd_above_insert_size         => { },
        old_sd_below_insert_size         => { },
        adaptor_path                     => { },
        archive_path                     => { },
        bam_path                         => { },
        gc_bias_path                     => { },
        analysis_software_version        => { },
        clusters                         => { },
        fwd_clusters                     => { },
        rev_clusters                     => { },
        old_fwd_filt_aligned_clusters_pct => { },
        old_rev_filt_aligned_clusters_pct => { },
        target_region_set_name           => { },
        old_filt_error_rate_avg          => { },
        fwd_seq_id                       => { },
        rev_seq_id                       => { },
        old_fwd_filt_error_rate_avg      => { },
        old_rev_filt_error_rate_avg      => { },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub entity_name { return 'instrument data solexa'; }

# TODO Require rebuild?
# library_id
# target_region_set_name
# gerald_directory
# TODO not updated. Fix?
# run_name [on flow_cell_illumina]
sub properties_to_copy {
    return ( 'id', 'library_id', properties_to_keep_updated() );
}

sub properties_to_keep_updated {
    return (qw/ 
        adaptor_path
        analysis_software_version
        bam_path
        clusters
        fastqc_path
        flow_cell_id
        fwd_clusters
        fwd_kilobases_read
        fwd_read_length
        fwd_run_type
        gerald_directory
        index_sequence
        is_external
        lane
        old_median_insert_size
        old_sd_above_insert_size
        old_sd_below_insert_size
        old_filt_error_rate_avg
        old_rev_filt_error_rate_avg
        old_fwd_filt_error_rate_avg
        old_rev_filt_aligned_clusters_pct
        old_fwd_filt_aligned_clusters_pct
        read_length
        rev_clusters
        rev_kilobases_read
        rev_read_length
        rev_run_type
        run_name
        run_type
        subset_name
        target_region_set_name
        /);
}

sub lims_property_name_to_genome_property_name {
    my ($class, $name) = @_;
    my %lims_to_genome = (
        filt_clusters => 'clusters',
        median_insert_size => 'old_median_insert_size',
        sd_above_insert_size => 'old_sd_above_insert_size',
        sd_below_insert_size => 'old_sd_below_insert_size',
        filt_error_rate_avg => 'old_filt_error_rate_avg',
        rev_filt_error_rate_avg => 'old_rev_filt_error_rate_avg',
        fwd_filt_error_rate_avg => 'old_fwd_filt_error_rate_avg',
        rev_filt_aligned_clusters_pct => 'old_rev_filt_aligned_clusters_pct',
        fwd_filt_aligned_clusters_pct => 'old_fwd_filt_aligned_clusters_pct',
    );
    return $lims_to_genome{$name} if exists $lims_to_genome{$name};
    return $name;
}

1;

