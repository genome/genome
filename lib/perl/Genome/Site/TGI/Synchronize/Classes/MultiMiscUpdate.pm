package Genome::Site::TGI::Synchronize::Classes::MultiMiscUpdate;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::MultiMiscUpdate {
    is => 'UR::Object',
    id_by => [
        subject_class_name => { is => 'Text', },
        subject_id => { is => 'Text', },
        description => { is => 'Text' },
        edit_date => { is => 'Text', },
    ],
    has_optional => [
        result => { is => 'Text', },
        #status => { is => 'Text', },
    ],
    has_many => [
        misc_updates => { is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate', },
    ],
};

sub perform_update {
    my $self = shift;
    my $method = '_'.lc($self->description);
    return $self->$method;
}

sub genome_entity_params {
    my $self = shift;

    my @misc_updates = $self->misc_updates;
    if ( not @misc_updates ) {
        $self->error_message('No misc updates to get genome entity params!');
        return;
    }

    my (%params, $lims_table_name);
    for my $misc_update ( $self->misc_updates ) {
        $params{ $misc_update->subject_property_name } = $misc_update->new_value;
        $lims_table_name = $misc_update->lims_table_name;
    }

    if ( $lims_table_name eq 'population_group_member' ) {
        $params{subject_id} = delete $params{pg_id};
        $params{attribute_label} = 'member';
        $params{attribute_value} = delete $params{member_id};
        $params{nomenclature} = 'WUGC',
    }
    elsif ( $lims_table_name eq 'sample_attribute' ) {
        $params{subject_id} = delete $params{organism_sample_id};
    }
    else {
        $self->error_message('Unsupported lims table name => '.$lims_table_name);
        return;
    }

    for my $required_key (qw/ subject_id attribute_label attribute_value nomenclature /) {
        next if defined $params{$required_key};
        $self->error_message("Missing required key ($required_key) in genome entity params!");
        return;
    }

    return %params;
}

sub _insert {
    my $self = shift;

    my %params = $self->genome_entity_params;
    return $self->_failure if not %params;

    my $genome_entity = Genome::SubjectAttribute->get(%params);
    if ( not $genome_entity ) {
        $genome_entity = Genome::SubjectAttribute->create(%params);
        if ( not $genome_entity ) {
            $self->error_message('Failed to create genome entity!');
            return $self->_failure;
        }
    }

    return $self->_success;
}

sub _delete {
    my $self = shift;

    my %params = $self->genome_entity_params;
    return $self->_failure if not %params;

    my $genome_entity = Genome::SubjectAttribute->get(%params);
    if ( $genome_entity ) {
        if ( not $genome_entity->delete ) {
            $self->error_message('Failed to delete genome entity!');
            return $self->_failure;
        }
    }

    return $self->_success;
}

sub _success {
    my $self = shift;
    for my $misc_update ( $self->misc_updates ) {
        $misc_update->success;
    }
    $self->result( $self->description );
    return 1;
}

sub _failure {
    my $self = shift;
    my $error_message = $self->error_message // 'NO ERROR SET!';
    for my $misc_update ( $self->misc_updates ) {
        $misc_update->error_message($error_message);
        $misc_update->failure;
    }
    $self->result('FAILED');
    return; 
}

1;

