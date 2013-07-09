package Genome::Site::TGI::Synchronize::Classes::Genotyping; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::Genotyping { # EXTERNAL 2877138689 ILLUMINA (INTERNAL) 2883765759
    table_name => <<SQL
    (
        select g.seq_id id, g.status status,
         g.organism_sample_id sample_id, s.full_name sample_name,
	     p.name platform_name, p.chip_type chip_name, p.version version,
         'external' import_source_name
	    from external_genotyping g
        join genotyping_platform p on p.genotyping_platform_id = g.genotyping_platform_id
        join organism_sample s on s.organism_sample_id = g.organism_sample_id
        union all
        select g.seq_id id, g.status status,
         g.organism_sample_id sample_id, s.full_name sample_name,
         'infinium' platform_name, 'HumanOmniExpress' chip_name, '12v1_A' version,
         'wugc' import_source_name
        from illumina_genotyping g
        join organism_sample s on s.organism_sample_id = g.organism_sample_id
        where g.status = 'pass'
    ) genotyping
SQL
    ,
    id_by => [
        id => { is => 'Text', },
    ],
    has => [
        status => { is => 'Text', },
        sample_id => { is => 'Text', },
        sample_name => { is => 'Text', },
        platform_name => { is => 'Text', },
        chip_name => { is => 'Text', },
        version => { is => 'Text', },
        import_source_name => { is => 'Text', },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

1;

