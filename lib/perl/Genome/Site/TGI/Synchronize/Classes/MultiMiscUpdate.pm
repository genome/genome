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
    has => [
        priority => { is => 'Number', },
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

sub genome_enitity_params {
    my $self = shift;

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
        Carp::confess('Unsupported lims table name => '.$lims_table_name);
    }

    return %params;
}

sub _insert {
    my $self = shift;

    my %params = $self->genome_enitity_params;
    my $genome_enitity = Genome::SubjectAttribute->get(%params);
    
    if ( not $genome_enitity ) {
        $genome_enitity = Genome::SubjectAttribute->create(%params);
        if ( not $genome_enitity ) {
            return $self->_failure('failed', 'Failed to create genome entity!');
        }
    }

    return $self->_success('inserted');
}

sub _delete {
    my $self = shift;

    my %params = $self->genome_enitity_params;
    my $genome_enitity = Genome::SubjectAttribute->get(%params);
    if ( $genome_enitity ) {
        if ( not $genome_enitity->delete ) {
            return $self->_failure('failed', 'Failed to delete genome entity!');
        }
    }

    return $self->_success('deleted');
}

sub _success {
    my $self = shift;
    for my $misc_update ( $self->misc_updates ) {
        $misc_update->success;
    }
    return 1;
}

sub _failure {
    my ($self, $error) = @_;
    for my $misc_update ( $self->misc_updates ) {
        $misc_update->failure($error);
    }
    return 1;
}

1;

