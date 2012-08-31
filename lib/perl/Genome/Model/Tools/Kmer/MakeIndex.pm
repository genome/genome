package Genome::Model::Tools::Kmer::MakeIndex;

use strict;
use warnings;

use Genome;
class Genome::Model::Tools::Kmer::MakeIndex {
    is => 'Genome::Model::Tools::Kmer',
    has => [
        index_name => {
            is => 'Text',
        },
        mer_size => {
            is => 'Integer',
            default_value => 20,
        },
        ],
    has_optional => [
            min_occ => {
            is => 'Integer',
             },
             
            max_occ => {
            is => 'Integer',
             }, 
    ],
};

sub execute {
    my $self = shift;

    my $gt_path = $self->genometools_path;
    my $options = '-mersize '. $self->mer_size;
    if ($self->min_occ) {
        $options .= ' -minocc '. $self->min_occ;
    }
    if ($self->max_occ) {
        $options .= ' -maxocc '. $self->max_occ;
    }
    my $cmd = $gt_path .' tallymer mkindex '. $options .' -esa '. $self->index_name;
    Genome::Sys->shellcmd(
        cmd => $cmd,
    );
    return 1;
}
