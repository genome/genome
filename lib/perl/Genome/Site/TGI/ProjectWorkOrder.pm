package Genome::Site::TGI::ProjectWorkOrder; 

use strict;
use warnings;

class Genome::Site::TGI::ProjectWorkOrder {
    id_by => [
        project_id => { is => 'Number', len => 10, column_name => 'PROJECT_ID' },
        setup_wo_id => { is => 'Number', column_name => 'SETUP_WO_ID' }
    ],
    has => [
        creation_event_id => { is => 'Number', column_name => 'CREATION_EVENT_ID' },
        administration_project => { is => 'Genome::Site::TGI::AdministrationProject', id_by => 'project_id' },
    ],
    doc => 'LIMS bridge table',
    data_source => 'Genome::DataSource::Oltp',
    table_name => 'project_work_order',
};


1;


#PROJECT_WORK_ORDER  GSC::ProjectWorkOrder   oltp    production
#
#    CREATION_EVENT_ID creation_event_id NUMBER(10)  (fk)    
#    PROJECT_ID        project_id        NUMBER(10)  (pk)(fk)
#    SETUP_WO_ID       setup_wo_id       NUMBER(10)  (pk)(fk)



