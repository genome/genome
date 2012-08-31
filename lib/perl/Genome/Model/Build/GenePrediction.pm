package Genome::Model::Build::GenePrediction;

use strict;
use warnings;

use Genome;
use Carp 'confess';

class Genome::Model::Build::GenePrediction {
    is => 'Genome::Model::Build',
    is_abstract => 1,
    subclassify_by => 'subclass_name',
    has => [
        subclass_name => {
            calculate_from => 'model_id',
            calculate => sub {
                my $model_id = shift;
                confess "Not given a model id with which to subclass the gene prediction build!" unless defined $model_id;
                my $model = Genome::Model->get($model_id);
                confess "Could not get model $model_id!" unless $model;
                my $subclass = $model->subclass_name;
                $subclass =~ s/Model/Model::Build/;
                return $subclass;
            },
        },
    ],
};

1;

