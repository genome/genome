package Genome::Model::Build::Command::DiffBlessed;

use strict;
use warnings;

use Genome;
use File::Spec;
use YAML;

class Genome::Model::Build::Command::DiffBlessed {
    is => 'Command::V2',
    roles => [ 'Genome::Role::Comparable::Command::Diff',
               'Genome::Role::Comparable::Command::DiffBlessed',
            ],
    has => [
        new_build => {
            is => 'Genome::Model::Build',
        },
        perl_version => {
            is => 'Text',
            valid_values => ['5.8', '5.10'],
            default_value => '5.10',
        },
    ],
};

sub help_brief : Overrides(Genome::Role::Comparable::Command::Diff) {
    &Genome::Role::Comparable::Command::Diff::help_brief;
}

sub help_detail : Overrides(Genome::Role::Comparable::Command::Diff) {
    &Genome::Role::Comparable::Command::Diff::help_detail;
}

sub diffs_message : Overrides(Genome::Role::Comparable::Command::Diff) {
    my $self = shift;
    my $diff_string = $self->Genome::Role::Comparable::Command::Diff::diffs_message(@_);
    my $bless_msg = $self->Genome::Role::Comparable::Command::DiffBlessed::bless_message(@_);
    return join("\n", $diff_string, $bless_msg);
}

sub blessed_object {
    my $self = shift;
    return $self->blessed_build;
}

sub new_object {
    my $self = shift;
    return $self->new_build;
}

sub blessed_build {
    my $self = shift;
    my $model_name = $self->new_build->model_name;
    return $self->retrieve_blessed_build($model_name);
}

sub retrieve_blessed_build {
    my ($self, $model_name) = @_;
    my $blessed_id = $self->blessed_id($model_name);
    my $blessed_build = Genome::Model::Build->get(
        id => $blessed_id,);
    return $blessed_build;
}

