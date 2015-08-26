package Genome::Role::Comparable::Command::DiffBlessed;

use strict;
use warnings;
use UR::Role;

use YAML qw();

role Genome::Role::Comparable::Command::DiffBlessed {
    requires => ['new_object'],
};

sub bless_message {
    my $self = shift;
    my $rel_db_file = $self->rel_db_file();
    my $new_object_id = $self->new_object->id;
    my $m = sprintf('If you want to bless this object (%s) update and commit the DB file (%s).', $new_object_id, $rel_db_file);
}

sub blessed_id {
    my ($self, $key) = @_;
    my ($blessed_ids) = YAML::LoadFile($self->db_file);
    my $blessed_id = $blessed_ids->{$key};
    unless(defined $blessed_id) {
        die "Undefined id returned for $key\n";
    }
    return $blessed_id;
}

sub db_file {
    my $self = shift;
    return $self->__meta__->module_path.".YAML";
}

sub rel_db_file {
    my $self = shift;
    my $db_file = $self->db_file;
    my $ns_base_dir = Genome->get_base_directory_name;
    my $rel_db_file = File::Spec->abs2rel($db_file, $ns_base_dir);
    return $rel_db_file;
}

1;

