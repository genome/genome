package Genome::Site::TGI::Synchronize::Classes::IndexIllumina;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::IndexIllumina {
    table_name => <<EOS
        (
            select
                --Index Illumina
                to_char(i.analysis_id) id,
                i.library_id library_id,
                i.index_sequence,
                i.lane,
                i.target_region_set_name,
                i.gerald_directory,
                i.median_insert_size,
                i.sd_above_insert_size,
                i.sd_below_insert_size,
                i.filt_clusters clusters,
                i.analysis_software_version,
                (case when i.index_sequence is null then to_char(i.lane) else to_char(i.lane) || '-' || i.index_sequence end) subset_name,

                --Constant
                0 is_external,

                --Flow Cell
                fc.run_name run_name,
                fc.flow_cell_id,

                --Read Illumina #1
                (case when r1.seq_id is not null then 'Paired' else 'Standard' end) run_type,

                --Read Illumina #2
                r2.read_length,
                r2.filt_error_rate_avg,

                --Fwd
                r1.sls_seq_id fwd_seq_id,
                (case when r1.seq_id is not null then 'Paired End Read 1' else null end) fwd_run_type,
                nvl(r1.read_length,-1) fwd_read_length,
                (case when r1.seq_id is not null then r2.kilobases_read else -1 end) fwd_kilobases_read,
                (case when r1.seq_id is not null then i.filt_clusters else null end) fwd_clusters,
                (case when r1.seq_id is not null then i.filt_clusters else null end) fwd_filt_clusters,
                r1.filt_aligned_clusters_pct fwd_filt_aligned_clusters_pct,
                r1.filt_error_rate_avg fwd_filt_error_rate_avg,

                --Rev
                (case when r1.seq_id is not null then r2.sls_seq_id else null end) rev_seq_id,
                (case when r1.seq_id is not null then 'Paired End Read 2' else null end) rev_run_type,
                (case when r1.seq_id is not null then r2.read_length else -1 end) rev_read_length,
                (case when r2.seq_id is not null then r2.kilobases_read else -1 end) rev_kilobases_read,
                (case when r1.seq_id is not null then i.filt_clusters else null end) rev_clusters,
                (case when r1.seq_id is not null then i.filt_clusters else null end) rev_filt_clusters,
                (case when r1.seq_id is not null then r2.filt_error_rate_avg else null end) rev_filt_error_rate_avg,
                (case when r1.seq_id is not null then r2.filt_aligned_clusters_pct else null end) rev_filt_aligned_clusters_pct,

                --Misc Paths
                archive2.path archive_path,
                gerald_bam.path bam_path,
                collect_gc_bias.path gc_bias_path,
                fastqc.path fastqc_path,
                '/gscmnt/sata114/info/medseq/adaptor_sequences/solexa_adaptor_pcr_primer'
                    || (case when sam.sample_type = 'rna' then '_SMART' else '' end) adaptor_path

                from GSC.index_illumina i
                    join GSC.flow_cell_illumina fc on fc.flow_cell_id = i.flow_cell_id
                    join GSC.read_illumina r2
                        on i.seq_id = r2.ii_seq_id
                        and (
                            (fc.run_type = 'Paired End' and r2.read_number = 2)
                            or
                            (fc.run_type = 'Fragment' and r2.read_number = 1)
                        )
                    left join GSC.seq_fs_path archive2 on archive2.seq_id = i.seq_id
                        and archive2.data_type = 'illumina fastq tgz'
                    left join GSC.seq_fs_path gerald_bam on gerald_bam.seq_id = i.seq_id
                        and gerald_bam.data_type = 'gerald bam'
                    left join GSC.seq_fs_path collect_gc_bias on collect_gc_bias.seq_id = i.seq_id
                        and collect_gc_bias.data_type = 'collect gc bias'
                    left join GSC.seq_fs_path fastqc on fastqc.seq_id = i.seq_id
                        and fastqc.data_type = 'fastqc'
                    left join GSC.read_illumina r1
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
        flow_cell_id                    => { }, # = short name
        lane                            => { },
        subset_name                     => { },
        index_sequence                  => { },
        run_name                        => { },
        run_type                        => { },
        fastqc_path                     => { },
        read_length                     => { },
        fwd_read_length                 => { },
        rev_read_length                 => { },
        fwd_kilobases_read              => { },
        rev_kilobases_read              => { },
        fwd_run_type                    => { },
        rev_run_type                    => { },
        gerald_directory                => { },
        median_insert_size              => { },
        sd_above_insert_size            => { },
        sd_below_insert_size            => { },
        adaptor_path                    => { },
        archive_path                    => { },
        bam_path                        => { },
        gc_bias_path                    => { },
        analysis_software_version       => { },
        clusters                        => { },
        fwd_clusters                    => { },
        rev_clusters                    => { },
        fwd_filt_aligned_clusters_pct   => { },
        rev_filt_aligned_clusters_pct   => { },
        target_region_set_name          => { },
        filt_error_rate_avg             => { },
        fwd_seq_id                      => { },
        rev_seq_id                      => { },
        fwd_filt_error_rate_avg         => { },
        rev_filt_error_rate_avg         => { },
    ],
    data_source => 'Genome::DataSource::GMSchema',
};

1;
