package Genome::Site::TGI::Synchronize::Expunge;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Expunge {
    is => 'Genome::Command::Base',
    has => [
        report => {
            is => 'Hashref',
            is_input => 1,
            is_optional => 0,
            doc => 'Hashref containing objects of interest from Genome::Site::TGI::Synchronze::UpdateApipeClasses',
        },
    ],
};

sub execute {
    my $self = shift;
    my %report = %{$self->report};

    my $dictionary = Genome::Site::TGI::Synchronize::Classes::Dictionary->get;
    for my $entity_name (keys %report){
        my $lims_class = $dictionary->lims_class_for_entity_name($entity_name);
        my $class = $lims_class->genome_class_for_create;
        next unless $class =~ m/Genome::InstrumentData/; #only remove instrument data for now
        next if $class eq 'Genome::InstrumentData::Imported'; #imported instrument data doesn't come from LIMS, so skip it
        my @ids = @{$report{$class}->{missing}} if $report{$class}->{missing};
        next if not @ids;
        printf("DELETING %s %s\n", $class, join(' ', @ids));
        my @deleted;
        for my $id (@ids){
            my $successfully_deleted = $self->_remove_expunged_object($class, $id);
            push @deleted, $successfully_deleted;
        }
        $report{$class}->{deleted} = \@deleted;
    }

    return 1;
}

sub _remove_expunged_object {
    my $self = shift;
    my $class = shift;
    my $id = shift;

    my $object = $class->get($id);

    $object->delete;

    return $id;
}

1;

