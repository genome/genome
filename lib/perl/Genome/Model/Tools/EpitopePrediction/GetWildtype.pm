package Genome::Model::Tools::EpitopePrediction::GetWildtype;

use strict;
use warnings;

use Genome;
use Workflow;
use Carp;

class Genome::Model::Tools::EpitopePrediction::GetWildtype {
    is => ['Genome::Model::Tools::EpitopePrediction::Base'],
    doc => "Get the Wildtype protein sequence from the specified Annotation Database for the variant proteins which have been annotated",
    has_input => [
        input_tsv_file => {
            is => 'Text',
            doc => 'A tab separated input file from the annotator',
        },
        output_directory => {
            is => 'Text',
            doc => 'Location of the output',
        },
        anno_db => {
            is => 'Text',
            is_optional=> 1,
            doc => 'The name of the annotation database.  Example: NCBI-human.combined-annotation',
        },
        anno_db_version => {
            is => 'Text',
            is_optional=> 1,
            doc => 'The version of the annotation database. Example: 54_36p_v2',
        },
    ],
    has_output => {
        output_tsv_file => {
            is => 'Text',
            doc => 'A tab separated output file with the amino acid sequences both wildtype and mutant',
            calculate_from => ['output_directory'],
            calculate => q| return File::Spec->join($output_directory, "snvs_wildtype.tsv"); |,
        },
    },
};

sub execute {
    my $self = shift;
    my $input = $self->input_tsv_file;
    my $output = $self->output_tsv_file;

    my $result = Genome::Model::Tools::Annotate::VariantProtein->execute(
        input_tsv_file  => $input,
        output_tsv_file => $output,
        anno_db         => $self->anno_db,
        version => $self->anno_db_version,
    );
    unless ($result) {
        confess $self->error_message("Couldn't execute Genome::Model::Tools::Annotate::VariantProtein $!");
    }

    return 1;
}

1;

__END__
gmt annotate variant-protein 
--input-tsv-file=Shared-Somatic-Tier1-Missense-d42m1.fullAnnotation_withHeader.tsv 
--output-tsv-file=Shared-Somatic-Tier1-Missense-d42m1.fullAnnotation_withHeader_WT.tsv 
--anno-db=NCBI-mouse.combined-annotation --anno-db-version=58_37k_v2
