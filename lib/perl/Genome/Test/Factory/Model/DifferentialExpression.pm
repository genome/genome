package Genome::Test::Factory::Model::DifferentialExpression;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::Build;
use Genome::Test::Factory::ProcessingProfile::DifferentialExpression;
use Genome::Test::Factory::ProcessingProfile::RnaSeq;
use Genome::Test::Factory::Model::ImportedAnnotation;
use Genome::Test::Factory::Model::ImportedReferenceSequence;
use Genome::Test::Factory::Model::RnaSeq;

our @required_params = qw(condition_labels_string condition_model_ids_string reference_sequence_build annotation_build);

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::DifferentialExpression->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::DifferentialExpression->setup_object();
    return $p->id;
}

sub setup_differential_expression_model {
    my $self = shift;
    my %params = @_;
    my $normal_model;
    my $tumor_model;
    my $normal_model_id = delete $params{normal_model_id};
    my $tumor_model_id = delete $params{tumor_model_id};
    if ($normal_model_id and $tumor_model_id) {
        $normal_model = Genome::Model->get($normal_model_id);
        $tumor_model = Genome::Model->get($tumor_model_id);
    }
    else {
        my $pp = Genome::Test::Factory::ProcessingProfile::RnaSeq->setup_object();
        my $reference = Genome::Test::Factory::Model::ImportedReferenceSequence->setup_object();
        my $reference_build = Genome::Test::Factory::Build->setup_object(model_id => $reference->id,
            status => "Succeeded");
        my $annot_model = Genome::Test::Factory::Model::ImportedAnnotation->setup_object();
        my $annot_build = Genome::Test::Factory::Build->setup_object(model_id => $annot_model->id,
            status => "Succeeded", reference_sequence => $reference_build);
        $normal_model = Genome::Test::Factory::Model::RnaSeq->setup_object(processing_profile_id => $pp->id,
            reference_sequence_build => $reference_build, annotation_build => $annot_build);
        my $normal_build = Genome::Test::Factory::Build->setup_object(model_id => $normal_model->id, status => "Succeeded");
        $tumor_model = Genome::Test::Factory::Model::RnaSeq->setup_object(
            processing_profile_id => $pp->id,
            subject => $normal_model->subject,
            reference_sequence_build => $reference_build,
            annotation_build => $annot_build,
        );
        my $tumor_build = Genome::Test::Factory::Build->setup_object(model_id => $tumor_model->id, status => "Succeeded");
        $normal_model_id = $normal_model->id;
        $tumor_model_id = $tumor_model->id;
    }
    my $condition_labels_string = "normal,tumor";
    my $condition_model_ids_string = $normal_model_id." ".$tumor_model_id;

    my $de_model = Genome::Test::Factory::Model::DifferentialExpression->setup_object(
        condition_labels_string => $condition_labels_string,
        condition_model_ids_string => $condition_model_ids_string,
        reference_sequence_build => $normal_model->reference_sequence_build,
        annotation_build => $normal_model->annotation_build,
        %params,
    );
    print STDERR "Created de model\n";
    return $de_model;
}

1;
