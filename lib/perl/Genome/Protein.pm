package Genome::Protein;
#:adukes short term: move data directory into id_by, but this has to be done in parallel w/ rewriting all file-based data sources.  It might be better to wait until long term: custom datasource that incorporates data_dir, possibly species/source/version, eliminating the need for these properties in the id, and repeated multiple times in the files

use strict;
use warnings;

use Genome;

class Genome::Protein {
    type_name => 'genome protein',
    table_name => 'PROTEIN',
    id_by => [
        protein_id => { is => 'Text' },
        species => { is => 'varchar',
            is_optional => 1,
        },
        source => { is => 'VARCHAR',
            is_optional => 1,
        },
        version => { is => 'VARCHAR',
            is_optional => 1,
        },
    ],
    has => [
        protein_name => { 
            is => 'String' 
        },
        transcript_id => { 
            is => 'Text' 
        },
         reference_build_id => {
            is => 'Text',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'reference_build_id',
        },
        amino_acid_seq => { 
            is => 'String' 
        },
        transcript => {
            is => 'Genome::Transcript', 
            calculate_from => ['data_directory', 'reference_build_id', 'transcript_id'],
            calculate => q/
               return Genome::Transcript->get(transcript_id => $transcript_id, data_directory => $data_directory, reference_build_id => $reference_build_id);
            /,
        },
        data_directory => {
            is => "Path",
        },
    ],
    schema_name => 'files',
    data_source => 'Genome::DataSource::Proteins',
};

1;

#TODO
=pod
=cut
