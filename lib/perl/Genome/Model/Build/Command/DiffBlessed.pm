package Genome::Model::Build::Command::DiffBlessed;

use strict;
use warnings;

use Genome;
use YAML;

class Genome::Model::Build::Command::DiffBlessed {
    is => 'Genome::Model::Build::Command::Diff::Base',
    has => [
        perl_version => {
            is => 'Text',
            valid_values => ['5.8', '5.10'],
            default_value => '5.10',
        },
        blessed_build => {
            is_optional => 1,
            calculate_from => ['new_build', 'perl_version', 'db_file'],
            calculate => q{ blessed_build($new_build->model->name, $perl_version, $db_file) },
        },
        db_file => {
            is_optional => 1,
            default_value => default_db_file(),
        },
    ],
};

sub default_db_file {
    return __FILE__.".YAML";
}

sub blessed_build {
    my ($model_name, $perl_version, $db_file) = @_;
    $db_file ||= default_db_file();
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
