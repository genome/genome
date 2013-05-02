package Genome::Model::Build::Command::DiffBlessed;

use strict;
use warnings;

use Genome;
use File::Spec;
use YAML;

class Genome::Model::Build::Command::DiffBlessed {
    is => 'Genome::Model::Build::Command::Diff::Base',
    has => [
        perl_version => {
            is => 'Text',
            valid_values => ['5.8', '5.10'],
            default_value => '5.10',
        },
    ],
};

sub db_file {
    return __FILE__.".YAML";
}

sub rel_db_file {
    my $self = shift;
    my $db_file = $self->db_file;
    my $ns_base_dir = Genome->get_base_directory_name;
    my $rel_db_file = File::Spec->abs2rel($db_file, $ns_base_dir);
    return $rel_db_file;
}

sub blessed_build {
    my $self = shift;
    my $model_name = $self->new_build->model_name;
    my $perl_version = $self->perl_version;
    my $db_file = $self->db_file;
    return retrieve_blessed_build($model_name, $perl_version, $db_file);
}

sub retrieve_blessed_build {
    my ($model_name, $perl_version, $db_file) = @_;
    $db_file ||= db_file();
    my ($blessed_ids) = YAML::LoadFile($db_file);
    my $blessed_id = $blessed_ids->{$model_name};
    unless(defined $blessed_id) {
        die "Undefined id returned for $model_name\n";
    }
    my $blessed_build = Genome::Model::Build->get(
        id => $blessed_id,);
    return $blessed_build;
}

sub diffs_message {
    my $self = shift;
    my $diff_string = $self->SUPER::diffs_message(@_);
    
    $diff_string = join("\n", $diff_string, sprintf('If you want to bless this build, update the file %s.', $self->db_file));

    return $diff_string;
}
