package Genome::SequenceDiffPart;

#:eclark 11/17/2009 Code review.

# 44 Rows in the database
# Used by Genome/Model/Tools/Htest/Diff/Define.pm and Define/Lst.pm

use strict;
use warnings;

use Genome;
class Genome::SequenceDiffPart {
    type_name => 'sequence diff part',
    table_name => 'SEQUENCE_DIFF_PART',
    id_by => [
        orig_position => { is => 'NUMBER', len => 10, sql => 'delete_position' },
        diff_id         => { is => 'NUMBER', len => 10 },
        refseq_path     => { is => 'VARCHAR2', len => 256 },
    ],
    has => [
        confidence_value => { is => 'NUMBER', len => 10, is_optional => 1 },
        orig_length    => { is => 'NUMBER', len => 10, is_optional => 1, sql => 'delete_length' },
        orig_sequence  => { is => 'VARCHAR2', len => 4000, is_optional => 1, sql => 'delete_sequence' },
        patched_length    => { is => 'NUMBER', len => 10, is_optional => 1, sql => 'insert_length' },
        patched_position  => { is => 'NUMBER', len => 10, is_optional => 1, sql => 'insert_position' },
        patched_sequence  => { is => 'VARCHAR2', len => 4000, is_optional => 1, sql => 'insert_sequence' },
        sequence_diff   => { is => 'Genome::SequenceDiff', id_by => 'diff_id', constraint_name => 'SDP_FK' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
