package Genome::Site::TGI::Synchronize::Classes::IndexIllumina;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::IndexIllumina {
    table_name => <<EOS
        (
            select
                --to_char(s_rev.seq_id) id,
                to_char(i.analysis_id) id,

                'solexa' sequencing_platform,

                --s_rev.sample_id,
                lib.sample_id,

                i.library_id,

                --s_rev.run_name,
                fc.run_name,

                fc.flow_cell_id,
                i.lane,

                r2.read_length,

                --(case when s_rev.run_type = 'Paired End Read 2' then s_rev.filt_aligned_clusters_pct else null end) rev_filt_aligned_clusters_pct,
                (case when r1.seq_id is not null then r2.filt_aligned_clusters_pct else null end) rev_filt_aligned_clusters_pct,

                --(case when s_fwd.run_type = 'Paired End Read 1' then s_fwd.filt_aligned_clusters_pct else null end) fwd_filt_aligned_clusters_pct,
                r1.filt_aligned_clusters_pct fwd_filt_aligned_clusters_pct,

                --(case when s_rev.run_type = 'Paired End Read 2' then s_rev.read_length else null end) rev_read_length,
                --(case when r1.seq_id is not null then r2.read_length else null end) rev_read_length,
                (case when r1.seq_id is not null then r2.read_length else -1 end) rev_read_length,

                (case when r1.seq_id is not null then r2.kilobases_read else -1 end) fwd_kilobases_read,
                (case when r2.seq_id is not null then r2.kilobases_read else -1 end) rev_kilobases_read,

                --(case when s_fwd.run_type = 'Paired End Read 1' then s_fwd.read_length else null end) fwd_read_length,
                --r1.read_length fwd_read_length,
                nvl(r1.read_length,-1) fwd_read_length,

                --(case when s_rev.run_type = 'Paired End Read 2' then s_rev.run_type else null end) rev_run_type,
                (case when r1.seq_id is not null then 'Paired End Read 2' else null end) rev_run_type,

                --(case when s_fwd.run_type = 'Paired End Read 1' then s_fwd.run_type else null end) fwd_run_type,
                (case when r1.seq_id is not null then 'Paired End Read 1' else null end) fwd_run_type,

                --(case when s_rev.run_type = 'Paired End Read 2' then 'Paired' else 'Standard' end) run_type,
                (case when r1.seq_id is not null then 'Paired' else 'Standard' end) run_type,

                --s_rev.gerald_directory,
                i.gerald_directory,

                --s_rev.median_insert_size,
                i.median_insert_size,

                --s_rev.sd_above_insert_size,
                i.sd_above_insert_size,
                
                --s_rev.sd_below_insert_size,
                i.sd_below_insert_size,

                --s_rev.is_external,
                0 is_external,

                --archive.path archive_path,
                archive2.path archive_path,
                gerald_bam.path gerald_bam_path,

                fastqc.path fastqc_path,

                --adaptor.path adaptor_path,
                --adaptor2.path adaptor_path,
                '/gscmnt/sata114/info/medseq/adaptor_sequences/solexa_adaptor_pcr_primer'
                    || (case when sample_type = 'rna' then '_SMART' else '' end) adaptor_path,

                --(case when s_fwd.run_type = 'Paired End Read 1' then s_fwd.FILT_CLUSTERS else null end) fwd_filt_clusters,
                (case when r1.seq_id is not null then i.filt_clusters else null end) fwd_filt_clusters,

                --(case when s_rev.run_type = 'Paired End Read 2' then s_rev.FILT_CLUSTERS else null end) rev_filt_clusters,
                (case when r1.seq_id is not null then i.filt_clusters else null end) rev_filt_clusters,

                --(nvl(s_fwd.FILT_CLUSTERS,0) + s_rev.FILT_CLUSTERS) filt_clusters, 	
                -- s_rev.FILT_CLUSTERS is still the expected value for fragment reads
                i.filt_clusters,
                    
                --s_rev.analysis_software_version,
                i.analysis_software_version,

                i.index_sequence

                --from GSC.solexa_lane_summary s_rev
                --join read_illumina r2 on r2.sls_seq_id = s_rev.seq_id --and r1.read_number = 1
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
        solexa_detail
EOS
    ,
    id_by => [
        id => { is => 'Text', },
    ],
    has_optional => [
        flow_cell_id                    => { }, # = short name
        flow_cell                       => { is => 'Genome::InstrumentData::FlowCell', id_by => 'flow_cell_id' },
        lane                            => { },
        fastqc_path                     => { },
        index_sequence                  => { },
        read_length                     => { },
        fwd_read_length                 => { },
        rev_read_length                 => { },
        fwd_kilobases_read              => { },
        rev_kilobases_read              => { },
        #TODO These three columns will point to "read_length" or whatever name is decided
        #(see also https://gscweb.gsc.wustl.edu/wiki/Software_Development:Illumina_Indexed_Runs_Warehouse_Schema)
        run_type                        => { },
        fwd_run_type                    => { },
        rev_run_type                    => { },
        gerald_directory                => { },
        median_insert_size              => { },
        sd_above_insert_size            => { },
        sd_below_insert_size            => { },
        is_external                     => { },
        adaptor_path                    => { },
        archive_path                    => { },
        bam_path                        => { column_name => 'gerald_bam_path'},
        analysis_software_version       => { },
        clusters                        => { column_name => 'filt_clusters' },
        fwd_clusters                    => { column_name => 'fwd_filt_clusters' },
        rev_clusters                    => { column_name => 'rev_filt_clusters' },
        fwd_filt_aligned_clusters_pct   => { },
        rev_filt_aligned_clusters_pct   => { },
    ],
    data_source => 'Genome::DataSource::GMSchema',
};

1;
