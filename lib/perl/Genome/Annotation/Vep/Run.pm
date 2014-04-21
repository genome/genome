package Genome::Annotation::Vep::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Vep::Run {
    is => 'Genome::Annotation::CommandBase',
    has_input => [
        ensembl_version => {
            is => 'String',
        },
        feature_list_ids_and_tags => {
            is => 'String',
            is_many => 1,
            doc => 'List of feature lists to be annotated in the 
                    INFO field, along with the tag to be used
                    e.g. 12345:SEGDUP,58676:ROI
                    The id and tag should be separated by a colon',
        },
        reference_build => {is => 'Genome::Model::Build::ReferenceSequence'},
        species => { is => 'Text', },
        polyphen => { is => 'String', },
        sift => { is => 'String', },
        condel => { is => 'String', },
        terms => {is => 'String',},
        regulatory => {is => 'Boolean',},
        gene => {is => 'Boolean',},
        most_severe => {is => 'Boolean',},
        per_gene => {is => 'Boolean',},
        hgnc => {is => 'Boolean',},
        coding_only => {is => 'Boolean',},
        canonical => {is => 'Boolean',},
        plugins => {is => 'String',
                    is_many => 1,
                    is_optional => 1},
        plugins_version => {is => 'String',},
        hgvs => {is => 'Boolean', },
    ],
};

sub execute {
    my $self = shift;

    $self->output_result(Genome::Annotation::Vep::RunResult->get_or_create($self->input_hash));
    return 1;
}
