package Genome::Model::MetagenomicCompositionShotgun::Command::ImportedDataArchives;

use Genome;
use strict;
use warnings;

class Genome::Model::MetagenomicCompositionShotgun::Command::ImportedDataArchives {
    is => 'Genome::Model::MetagenomicCompositionShotgun::Command',
    doc => 'print path to imported data archives',
    has => [
        model_id => {
            is => 'Integer',
        },
    ],
};

sub execute {
    my $self = shift;

    my $mcs_model = Genome::Model->get($self->model_id);
    my ($meta_model) = $mcs_model->_metagenomic_alignment_models;
    my @imported_data = $meta_model->instrument_data;
    print "WARNING: These archives are in an unverified quality format do not assume it is in Sanger quality format.\n";
    for my $imported_data (@imported_data) {
        print $imported_data->archive_path . "\n";
    }
}

1;
