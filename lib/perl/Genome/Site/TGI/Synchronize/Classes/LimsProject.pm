package Genome::Site::TGI::Synchronize::Classes::LimsProject; 

use strict;
use warnings;
use Genome;

class Genome::Site::TGI::Synchronize::Classes::LimsProject {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    table_name => <<SQL
    (
		--Administration Project
		select ap.project_id id, ap.project_name name, 'na' pipeline 
		from administration_project ap
		where ap.project_id > 2570000
         and ap.status != 'abandoned'
		union all
		--Setup Work Order
		select wo.setup_wo_id id, s.setup_name name, wo.pipeline pipeline
		from setup s
		join setup_work_order wo on wo.setup_wo_id = s.setup_id
		where s.setup_id > 2570000
		 and s.setup_status != 'abandoned'
    ) lims_project
SQL
    ,
    id_by => [
        id => { is => 'Text', },
    ],
    has => [
        name => { is => 'Text', },
        pipeline => { is => 'Text', },
    ],
    data_source => 'Genome::DataSource::Oltp',
};

sub entity_name { return 'project'; }

sub params_for_create_in_genome {
    my $self = shift;
    return ( map { $_ => $self->$_ } (qw/ id name /) );
}

sub __display_name__ {
    return $_[0]->name.' ('.$_[0]->id.')';
}

1;

