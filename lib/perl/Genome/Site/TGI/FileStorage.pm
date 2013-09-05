package Genome::Site::TGI::FileStorage;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::FileStorage {
    is => 'UR::Object',
    table_name => '(SELECT file_storage_id, file_name FROM gsc.file_storage) file_storage',
    id_by => [
        file_storage_id => {
            is => 'Number', len => '20',
        },
    ],
    has => [
        file_name => {
            is => 'VARCHAR2', len => '255',
        },
        #cannot db-link LOB columns (see ORA-22992), so this must go through the GSC object or use an OLTP datasource
        content => {
            is => 'BLOB',
            calculate_from => ['_gsc_file_storage'],
            calculate => q{
                $_gsc_file_storage->content,
            },
        },
        _gsc_file_storage => {
            is => 'GSC::FileStorage',
            calculate_from => ['file_storage_id'],
            calculate => q{
                GSC::FileStorage->get(file_storage_id => $file_storage_id);
            },
        },
    ],
    data_source => 'Genome::DataSource::Oltp',
};

1;
