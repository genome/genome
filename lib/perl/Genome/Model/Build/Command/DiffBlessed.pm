package Genome::Model::Build::Command::DiffBlessed;

use strict;
use warnings;

use Genome;

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
            calculate_from => ['new_build', 'perl_version'],
            calculate => q{
                $DB::single = 1;
                my $model_id = $new_build->model->id;
                my $blessed_build_raw = qx(/gsc/scripts/opt/genome/bin/list-blessed-builds -m $model_id -p $perl_version);
                my $commit = (split(/\s+/, $blessed_build_raw))[2];
                my $blessed_build = Genome::Model::Build->get(
                    model_id => $model_id,
                    software_revision => "$perl_version-$commit");
                return $blessed_build;
            },
        },
    ],
};
