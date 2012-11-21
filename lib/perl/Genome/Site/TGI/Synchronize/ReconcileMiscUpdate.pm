package Genome::Site::TGI::Synchronize::ReconcileMiscUpdate;

use strict;
use warnings;

use Genome;

use Mail::Sender;

class Genome::Site::TGI::Synchronize::ReconcileMiscUpdate {
    is => 'Command::V2',
    has => [
    ],
};

sub table_to_class_mapping {
    return (
        "gsc.organism_taxon" => 'Genome::Taxon',
        "gsc.organism_individual" => 'Genome::Individual',
        "gsc.organism_population_group" => 'Genome::PopulationGroup',
        "gsc.organism_sample" => 'Genome::Sample',
        "gsc.sample_attribute" => 'Genome::SubjectAttribute',
        "gsc.population_group_member" => 'Genome::SubjectAttribute',
    );
}

sub class_for_table {
    my ($self, $table) = @_;
    my %m = $self->table_to_class_mapping;
    return $m{$table};
}

sub execute {
    my $self = shift;

    my @msg;

    push @msg, $self->sync_properties_that_map_directly();

#    push @msg, $self->sync_work_order_project();
    my $sender = Mail::Sender->new({ smtp => 'gscsmtp.wustl.edu', from => 'jlolofie@genome.wustl.edu', });

    $sender->MailMsg({
        to => 'apipe-builder@genome.wustl.edu',
        subject => 'lims -> apipe sync',
        msg => join("\n",@msg)
    });

}

sub sync_work_order_project {

    my ($self) = @_;
    # TODO: here
}

sub sync_properties_that_map_directly {

    my ($self) = @_;

    my $table_for_class = $self->table_to_class_mapping();
    my $msg;

    for my $table (keys %$table_for_class) {
        my @updates = Genome::Site::TGI::MiscUpdate->get(is_reconciled => 0, table_name => $table);
        my $class = $table_for_class->{$table};
        if (!$class) {
            $self->warning_message("Skipping table $table (no class mapping)");
            next;
        }
        $self->status_message("Found " . scalar @updates . " unreconciled updates for $table ($class)");
        my $updates_to_review = $self->process_updates($table, $class, \@updates);

        $msg .= "Table: $table Class: $class\n";
        $msg .= Data::Dumper::Dumper $updates_to_review;
    }

    return $msg;
}

sub process_updates {

    my ($self, $table, $class, $updates) = @_;

    my $updates_to_review = {};
    my @completed_updates;
    my $update_list = {};

    for my $update (@$updates) {

        my $column = $update->column_name;
        unless ($column) {
            $self->warning_message("Update " . $update->__display_name__ . " has no column name, cannot apply!");
            $update->is_reconciled(1);
            next;
        }

        my $lims_old_value = $update->old_value;
        my $lims_new_value = $update->new_value;

        if ($lims_old_value && $lims_new_value) {
            next if $lims_old_value eq $lims_new_value;
        }

        my $object = $class->get($update->updated_row_primary_key);
        unless ($object) {
            $self->warning_message("Found no object of class $class with id " . $update->updated_row_primary_key);
#            push @{ $updates_to_review->{'could_not_find_genome_object'} }, join("\t",$class,$update->updated_row_primary_key);
            $updates_to_review->{'could_not_find_genome_object'}++;
            $update->is_reconciled(1);
            next;
        }

        if (! grep { $_->property_name eq $column } $class->__meta__->properties() ) {
#            $self->warning_message("Column '$column' doesnt exist for class '$class'");
#            push @{ $updates_to_review->{'property_doesnt_exist_for_column'} }, join("\t",$class,$column,$update->id);
            $updates_to_review->{'property_doesnt_exist_for_column'}++;
            $update->is_reconciled(1);
            next;
        }

        my $apipe_current_value = $object->$column;
        if ($apipe_current_value and $apipe_current_value eq $lims_new_value) {
            $update->is_reconciled(1);
            next;
        }

        if ( ! $lims_old_value or ($apipe_current_value ne $lims_old_value) ) {
            $self->warning_message("Object " . $object->__display_name__ . " property $column " . 
                "currently has value '$apipe_current_value', but update claims value should be '$lims_old_value'. " . 
                "Cannot safely apply update, marking for later review and skipping");
            $updates_to_review->{'current_value_doesnt_match_old'}++;
#            push @{ $updates_to_review->{'current_value_doesnt_match_old'} }, join("\t",$class,$column,$update->id);
            next;
        }

#        print "JTAL: would have changed $column to '$new_value' from '$current_value'\n";

        $update_list->{$column}++;
        printf("%s\t%s\t%s\t%s\t%s\t%s\n", $apipe_current_value, $lims_old_value, $lims_new_value, $class, $object->id, $column);
#        $object->$column($new_value);
#        $update->is_reconciled(1);
        push @completed_updates, $update;
    }

    print $table . " -> " . $class . "\n";
    print Data::Dumper::Dumper $update_list;
    return $updates_to_review; 
}

1;

