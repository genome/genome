package Genome::Site::TGI::Synchronize::Classes::PopulationGroupMember; 

use strict;
use warnings;

=pod
MEMBER_ID NUMBER   (10)                    {null} NOT NULL ok
PG_ID     NUMBER   (10)                    {null} NOT NULL ok [population_group_id]
properties 2, copied NA, updated NA
=cut

class Genome::Site::TGI::Synchronize::Classes::PopulationGroupMember {
    is => 'UR::Object',
    table_name => 'POPULATION_GROUP_MEMBER',
    id_by => [
        population_group  => { is => 'Genome::Site::TGI::Synchronize::Classes::PopulationGroup', id_by => 'pg_id' },
        member_id => { is => 'Number', },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

1;

