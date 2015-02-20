package Genome::Interfaces::Comparable::Command::DiffBlessed;

use strict;
use warnings;
use Genome;

use YAML;

class Genome::Interfaces::Comparable::Command::DiffBlessed {
    is => 'Genome::Interfaces::Comparable::Command::Diff',
    is_abstract => 1,
};

sub bless_message {
    my $self = shift;
    my $rel_db_file = $self->rel_db_file();
    my $new_object_id = $self->new_object->id;
    my $m = sprintf('If you want to bless this object (%s) update and commit the DB file (%s).', $new_object_id, $rel_db_file);
}

sub diffs_message {
    my $self = shift;
    my $diff_string = $self->SUPER::diffs_message(@_);
    my $bless_msg = $self->bless_message();
    $diff_string = join("\n", $diff_string, $bless_msg);
    return $diff_string;
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

