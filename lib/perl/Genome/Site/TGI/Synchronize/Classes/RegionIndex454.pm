package Genome::Site::TGI::Synchronize::Classes::RegionIndex454;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::RegionIndex454 {
    table_name => <<'SQL'
        (
            select 
                --ri454
                to_char(ri454.seq_id) id,
                ri454.index_sequence index_sequence,
                ri454.library_id,
                ri454.num_bases total_bases_read,
                ri454.num_reads total_reads,
                --rr454
                rr454.paired_end is_paired_end,
                rr454.region_id,
                rr454.region_number,
                rr454.run_name run_name,
                ( 
                 case when ri454.index_sequence is null
                  then to_char(rr454.region_number)
                  else to_char(rr454.region_number) || '-' || ri454.index_sequence 
                 end
                ) subset_name,
                --constants
                '454' sequencing_platform
            from GSC.region_index_454 ri454
            join GSC.run_region_454 rr454 on rr454.region_id = ri454.region_id
        ) region_index_454
SQL
    ,
    id_by => [
        id => { is => 'Text', },
    ],
    has_optional => [
        index_sequence => { is => 'Text', },
        is_paired_end => { is => 'Text', },
        library_id => { is =>'Text', },
        region_id => { is => 'Text', },
        region_number => { is => 'Text', },
        run_name => { is => 'Text', },
        sequencing_platform => { is => 'Text', },
        subset_name => { is => 'Text', },
        total_reads => { is => 'Text', },
        total_bases_read => { is => 'Text', },
    ],
    data_source => 'Genome::DataSource::GMSchema',
};

