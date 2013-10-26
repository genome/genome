package Genome::WorkOrder;
use strict;
use warnings;
use Genome;

class Genome::WorkOrder {
    is => 'Genome::Searchable',
    table_name => '(SELECT swo.*, 
                           s.*, 
                           s.setup_name as name,
                           s.setup_description as description
                      FROM setup_work_order swo
                      JOIN setup s
                        ON s.setup_id = swo.setup_wo_id
                    ) work_order',
    id_by => [
        id => {
            is => 'Integer',
            len => 10,
            column_name => 'SETUP_WO_ID',
        },
    ],    
    has => [
        name => {
             is => 'Text',
            len => 32,
        },
        description => {
             is => 'Text',
            len => 255,
        },
        acct_id => {
            is => 'Integer',
            len => 10,
        },
        wo_facilitator => {
            is => 'Integer',
            len => 10,
        },
        deadline => {
            is => 'Date',
        },
        project_tracking_number => {
            is => 'Text',
            len => 32,
        },
        pipeline => {
            is => 'Text',
            len => 256,
        },
        project_id => {
            is => 'Text',
            len => 10,
        },
        project => { 
            is => 'Genome::Site::TGI::Project', 
            id_by => 'project_id' 
        },
        project_name => {
            via => 'project',
            to => 'name',
            is_many => 0
        },
        collaborator => {
            via => 'project',
            to => 'external_contact_name',
            is_many => 0
        },
        estimate_id => {
            is => 'Text',
            len => 32,
        },
        is_test => {
            is => 'Integer',
            len => 1,
        },
        #setup_ss_id => {
        #    is => 'Integer',
        #    len => 10,
        #},
        barcode => {
            is => 'Text',
            len => 16,
        },
        requester_con_id => {
                is => 'Integer',
                len => 10,
        },
        items => {
            is => 'Genome::WorkOrderItem',
            is_many => 1,
            reverse_as => 'work_order',
        },
        samples => {
            is => 'Genome::Sample',
            via => 'items',
            to => 'sample',
            is_many => 1,
        },
        models => {
            is => 'Genome::Model',
            via => 'items',
        },
        builds => {
            is => 'Genome::Model::Build',
            via => 'models',
        },
    ],
    data_source => 'Genome::DataSource::Oltp',
};

sub xitems {
    return Genome::WorkOrderItem->get(setup_wo_id => $_[0]->id);
}

sub sample_description {
    my ($self) = @_;
    my @samples = $self->samples();
    my @names = map {$_->sample} @samples;
    return join(',',@names);
}

1;

