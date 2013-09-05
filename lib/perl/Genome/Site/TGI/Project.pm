package Genome::Site::TGI::Project; 
use strict;
use warnings;
use Genome;

# NOTE: this replaces the old Genome::Project
# The new Genome::Project is built from this plus analyst-created projects

class Genome::Site::TGI::Project {
    is => 'Genome::Searchable',
    id_properties => ['setup_project_id'],
    table_name =>   "(SELECT * FROM setup_project p JOIN setup s ON setup_id = setup_project_id WHERE project_type != 'setup project finishing' AND setup_status != 'abandoned') project ",
    has => [
        
        # get these directly, since we can't join through any app objects
        name                => { is => 'Text', len => 64, column_name => 'SETUP_NAME' },
        status              => { is => 'Text', column_name => 'SETUP_STATUS' },
        description         => { is => 'Text', len => 256, column_name => 'SETUP_DESCRIPTION' },
        project_type        => { is => 'Text', len => 32 },
        mailing_list        => { is => 'Text', column_name => 'MAILING_LIST' }, 
        
        # not used, but available
        acct_id             => { is => 'Number', len => 10 },
        het_testing         => { is => 'Boolean' },
        is_submitted        => { is => 'Boolean' },
        parent_project_id   => { is => 'Number', len => 10 },
        setup_project_id    => { is => 'Number', len => 10 },
    ],
    has_optional => [ 
        work_orders         => { is => 'Genome::WorkOrder', reverse_as => 'project', is_many => 1 },
        samples             => { is => 'Genome::Site::TGI::Sample', via => 'work_orders', to => 'samples' },
        models              => { is => 'Genome::Model', via => 'samples', to => 'models' },

        external_contact        => { is => 'Genome::Site::TGI::ProjectContact', id_by => 'ext_con_id' },
        external_contact_name   => { is => 'Text', via => 'external_contact', to => 'name' }, 
        external_contact_email  => { is => 'Text', via => 'external_contact', to => 'email' }, 
       
        internal_contact        => { is => 'Genome::Site::TGI::ProjectContact', id_by => 'internal_con_id' },
        internal_contact_name   => { is => 'Text', via => 'internal_contact', to => 'name' }, 
        internal_contact_email  => { is => 'Text', via => 'internal_contact', to => 'email' }, 
    ],
    data_source => 'Genome::DataSource::Oltp',
};

1;


