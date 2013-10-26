package Genome::Site::TGI::AdministrationProject; 

use strict;
use warnings;

class Genome::Site::TGI::AdministrationProject {
    id_by => [
        id => { is => 'Number', len => 10, column_name => 'PROJECT_ID' },
    ],
    has => [
        creation_event_id => { is => 'Number', column_name => 'CREATION_EVENT_ID' },
        priority => { is => 'Text', column_name => 'PRIORITY' },
        project_name => { is => 'Text', column_name => 'PROJECT_NAME', doc => 'LIMS project_name' },
        project_work_orders => {
            is => 'Genome::Site::TGI::ProjectWorkOrder',
            reverse_as => 'administration_project',
            is_many => 1,
        },
        analysis_project_ids => {
            via => 'project_work_orders',
            to => 'setup_wo_id',
            is_many => 1,
        },
        analysis_projects => {
            calculate => q {
                return Genome::Project->get(id => [$self->analysis_project_ids]);
            }
        }
    ],
    has_optional => [
        parent_project_id => { is => 'Number', column_name => 'PARENT_PROJECT_ID' },
        status => { is => 'Text', column_name => 'STATUS' }
    ],
    doc => 'yet another table about projects',
    data_source => 'Genome::DataSource::Oltp',
    table_name => 'administration_project',
};


1;


#ADMINISTRATION_PROJECT  GSC::AdministrationProject  oltp    production
#
#    BACKGROUND        background        BLOB(2147483647) NULLABLE         
#    CREATION_EVENT_ID creation_event_id NUMBER(10)                (fk)    
#    OVERVIEW          overview          BLOB(2147483647) NULLABLE         
#    PARENT_PROJECT_ID parent_project_id NUMBER(10)       NULLABLE (fk)    
#    PRIORITY          priority          VARCHAR2(16)                      
#    PROJECT_ID        project_id        NUMBER(10)                (pk)    
#    PROJECT_NAME      project_name      VARCHAR2(100)             (unique)
#    PROJECT_UPDATE    project_update    BLOB(2147483647) NULLABLE         
#    STATUS            status            VARCHAR2(32)     NULLABLE  


