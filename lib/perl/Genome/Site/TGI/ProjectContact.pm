package Genome::Site::TGI::ProjectContact; 
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::ProjectContact {
    table_name =>   "(select * from contact\@oltp) contact",
    id_properties => [
        con_id      => { is => 'Number', len => 10 },
    ],
    has => [
        email       => { is => 'Email', column_name => 'CONTACT_EMAIL' },
        name        => { is => 'Text',  column_name => 'CONTACT_NAME'  },
        type        => { is => 'Text',  column_name => 'CONTACT_TYPE'  },
    ],
    data_source => 'Genome::DataSource::GMSchema',
};

1;

