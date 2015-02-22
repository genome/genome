package Genome::InstrumentData::Solexa::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Solexa::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                {
                    name => 'flow_cell',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'flow_cell_id',
                        'machine_name',
                        'run_name',
                        'run_type',
                        {
                            name => 'lanes',
                            perspective => 'default',
                            toolkit => 'xml',
                            aspects => [
                                'id',
                            ],
                        },
                    ],
                },
                {
                    name => 'library',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'library_id',
                        'name',
                    ],
                },
                {
                    name => 'sample',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name',
                        'common_name',
                        'extraction_label',
                        'extraction_type',
                        'extraction_desc',
                        'tissue_label',
                        'tissue_desc',
                        'organ_name',
                        {
                            name => 'taxon',
                            aspects => ['id', 'species_name'],
                            perspective => 'default',
                            toolkit => 'xml'
                        },
                        {
                            name => 'source',
                            aspects => ['id', 'name', 'common_name'],
                            perspective => 'default',
                            toolkit => 'xml'
                        }
                    ],
                },
                {
                    name => 'taxon',
                    aspects => [
                        'id',
                        'domain',
                        'species_name',
                        'domain',
                        'species_latin_name',
                        'strain_name',
                        'ncbi_taxon_id',
                    ],
                    perspective => 'default',
                    toolkit => 'xml'
                },
                'archive_path',
                'fwd_clusters', 'rev_clusters', 'clusters',
                'fwd_read_length', 'rev_read_length', 'read_length',
                'subset_name',
                'target_region_set_name',
                'index_sequence',
                'project_name',
            ]
        }
    ]
};

1;

