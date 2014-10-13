package Genome::Site::TGI::Synchronize::CheckForNameUpdates;

use strict;
use warnings;

use Genome;

use constant TABLES_TO_PROCESS => [qw(gsc.organism_sample gsc.library_summary)];

class Genome::Site::TGI::Synchronize::CheckForNameUpdates {
    is => 'Command::V2',
    has_optional => [
        tables_to_process => {
            is => 'Text',
            is_many => 1,
            valid_values => TABLES_TO_PROCESS,
            default => TABLES_TO_PROCESS,
            doc => 'which LIMS tables to consider updates to',
        },
        limit => {
            is => 'Number',
            doc => 'How many updates to consider',
        },
        updates => {
            is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate',
            is_many => 1,
            doc => 'Process these specific updates, bypassing the normal query',
        },
    ]
};

sub execute {
    my $self = shift;

    my @updates = $self->updates;
    unless(@updates) {
        @updates = Genome::Site::TGI::Synchronize::Classes::MiscUpdate->get(
            subject_class_name => [$self->tables_to_process],
            subject_property_name => 'full_name',
            description => 'UPDATE',
            is_reconciled => 0,
            ($self->limit? (-limit => $self->limit) : ()),
        );
    }

    if(@updates) {
        $self->status_message('Found %s unreconciled updates.', scalar(@updates));
    } else {
        $self->status_message('No unreconciled updates found.');
        return 1;
    }

    $self->_preload_entities(@updates);

    for my $update (@updates) {
        $self->process_update($update);
    }

    my @failures = grep { $_->has_failed } @updates;
    if(@failures) {
        $self->status_message('Found %s updates that need attention.', scalar(@failures));
        for my $failure (@failures) {
            $self->status_message('Need %s renamed to %s', $failure->old_value, $failure->new_value);
        }
    } else {
        $self->status_message('No updates need attention.');
    }

    return 1;
}

sub process_update {
    my $self = shift;
    my $update = shift;

    unless($self->update_is_current($update)) {
        $self->debug_message('Found superseded update for %s of %s', $update->subject_id, $update->subject_class_name);
        return $update->success(); #don't need to process an update which has been updated again already
    }

    if($self->update_has_occurred($update)) {
        $self->debug_message('Found completed update for %s of %s', $update->subject_id, $update->subject_class_name);
        return $update->success(); #update is already handled
    } else {
        $self->debug_message('Found pending update for %s of %s', $update->subject_id, $update->subject_class_name);
        return $update->failure(); #update requires intervention
    }
}

sub update_is_current {
    my $self = shift;
    my $update = shift;

    my $foreign_entity = $update->site_tgi_class_name->get($update->subject_id);
    return unless $foreign_entity;

    my $property_name = $update->genome_property_name;
    return ($foreign_entity->$property_name eq $update->new_value);
}

sub update_has_occurred {
    my $self = shift;
    my $update = shift;

    my $genome_entity = $update->genome_entity;
    return unless $genome_entity;

    my $property_name = $update->genome_property_name;
    return ($genome_entity->$property_name eq $update->new_value);
}

sub _preload_entities {
    my $self = shift;
    my @updates = @_;

    my %classes_and_ids;

    for my $update (@updates) {
        $classes_and_ids{$update->site_tgi_class_name} ||= [];
        push @{$classes_and_ids{$update->site_tgi_class_name}}, $update->subject_id;
        $classes_and_ids{$update->genome_class_name} ||= [];
        push @{$classes_and_ids{$update->genome_class_name}}, $update->subject_id;
    }

    for my $class (keys %classes_and_ids) {
        $class->get($classes_and_ids{$class});
    }

    return 1;
}

1;
