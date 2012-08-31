package Genome::Site::TGI::Sample::Genomic;

use strict;
use warnings;
use Genome;

class Genome::Site::TGI::Sample::Genomic {
    table_name => "
        (select dna.*,o.taxon_id 
         from dna\@oltp
         left join (
	            dna_resource\@oltp dr 
	            join entity_attribute_value\@oltp eav
		    on eav.entity_id = dr.dr_id
		        and eav.type_name = 'dna'
		        and eav.attribute_name = 'org id'	
	            join GSC.organism_taxon o 
		        on o.legacy_org_id = eav.value
                   ) on dr.dna_resource_prefix = substr(dna_name,0,4)
         where dna_type = 'genomic dna'
         and taxon_id is not null
) genomic_sample",
    id_by => [
              id => { is => 'Number', column_name => 'DNA_ID' },
    ],
    has => [
            name => { is => 'Text', len => 64, column_name => 'DNA_NAME' },
            taxon                       => { is => 'Genome::Site::TGI::Taxon', id_by => 'taxon_id' },
            species_name                => { via => 'taxon' },
        ],
    data_source => 'Genome::DataSource::GMSchema',
};
