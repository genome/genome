package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::SubjectAttribute;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::SubjectAttribute {
    is => 'UR::Object',
    id_by => [
        subject_class_name => { is => 'Text', },
        subject_id => { is => 'Text', },
        description => { is => 'Text', },
        edit_date => { is => 'Text', },
    ],
    has_optional => [
        result => { is => 'Text', },
        status => { is => 'Text', },
    ],
    has_many => [
        misc_updates => { is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate', },
    ],
    has_constant_calculated => [
        value_method => {
            calculate => q|
                my %descriptions_and_methods = (
                    INSERT => 'new_value',
                    DELETE => 'old_value',
                );
                return $descriptions_and_methods{ $self->description };
            |,
        },
    ],
};

sub get_or_create_from_misc_updates {
    my ($class, @misc_updates) = @_;

    if ( not @misc_updates ) {
        $class->error_message('No misc updates to get or create multi misc update!');
        return;
    }

    my %params = map { $_ => $misc_updates[0]->$_ } $class->__meta__->id_property_names;
    my $self = $class->get(%params);
    if ( not $self ) {
        $self = $class->create(%params);
    }

    for my $misc_update ( @misc_updates ) {
        return if not $self->add_misc_update($misc_update);
    }

    return $self;
}

sub add_misc_update {
    my ($self, $misc_update) = @_;

    if ( not $misc_update ) {
        $self->error_message('No misc update given to add!');
        return;
    }

    for my $attr ( $self->__meta__->id_property_names ) {
        next if defined $self->$attr and $self->$attr eq $misc_update->$attr;
        $self->error_message("Mismatch in id property ($attr) when adding misc update! Values do not match! ".$self->$attr." vs ".$misc_update->$attr);
        return;
    }

    return $self->SUPER::add_misc_update($misc_update);
}

sub perform_update {
    my $self = shift;
    my $method = '_'.lc($self->description);
    return $self->$method;
}

sub _insert {
    my $self = shift;

    my %params = $self->_resolve_genome_entity_params;
    return $self->_failure if not %params;

    my $genome_subject = Genome::Subject->get(id => $params{subject_id});
    if ( not $genome_subject ) {
        $self->status_message('No genome subject for id! '. $params{subject_id});
        return $self->_skip;
    }

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

    my %params = $self->_resolve_genome_entity_params;
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

sub _resolve_genome_entity_params {
    my $self = shift;

    my @misc_updates = $self->misc_updates;
    if ( not @misc_updates ) {
        $self->error_message('No misc updates set to get genome entity params!');
        return;
    }

    my $value_method = $self->value_method;
    if ( not $value_method ) {
        $self->error_message('Failed to get value method!');
        return;
    }

    my (%params);
    for my $misc_update ( $self->misc_updates ) {
        $params{ $misc_update->subject_property_name } = $misc_update->$value_method;
    }

    my $lims_table_name = $misc_updates[0]->lims_table_name;
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

sub _set_result {
    my ($self, $result) = @_;

    $self->result($result);

    my $value_method = $self->value_method;
    $self->status(
        join(
            "\t", 
            $result, 
            map({ $self->$_; } (qw/ description subject_class_name subject_id /)),
            map({ $_->$value_method } $self->misc_updates),
        )
    );
}

sub _success {
    my $self = shift;
    for my $misc_update ( $self->misc_updates ) {
        $misc_update->success;
    }
    $self->_set_result('PASS');
    return 1;
}

sub _skip {
    my $self = shift;
    for my $misc_update ( $self->misc_updates ) {
        $misc_update->skip;
    }
    $self->_set_result('SKIP');
    return;
}

sub _failure {
    my $self = shift;
    my $error_message = $self->error_message // 'NO ERROR MSG SET!';
    for my $misc_update ( $self->misc_updates ) {
        $misc_update->error_message($error_message);
        $misc_update->failure;
    }
    $self->_set_result('FAIL');
    return; 
}

1;

