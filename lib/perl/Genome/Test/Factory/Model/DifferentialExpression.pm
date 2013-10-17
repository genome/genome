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
    my $condition_labels_string = "normal,tumor";
    my $pp = Genome::Test::Factory::ProcessingProfile::RnaSeq->setup_object();
    my $reference = Genome::Test::Factory::Model::ImportedReferenceSequence->setup_object();
    my $reference_build = Genome::Test::Factory::Build->setup_object(model_id => $reference->id,
                                                                     status => "Succeeded");
    print STDERR "Created reference build\n";
    my $m = Genome::Test::Factory::Model::ImportedAnnotation->setup_object();
    my $annot_build = Genome::Test::Factory::Build->setup_object(model_id => $m->id,
                                                      status => "Succeeded", reference_sequence => $reference_build);
    print STDERR "Created annotation build\n";
    my $normal = Genome::Test::Factory::Model::RnaSeq->setup_object(processing_profile_id => $pp->id,
                                                                    reference_sequence_build => $reference_build);
    my $normal_build = Genome::Test::Factory::Build->setup_object(model_id => $normal->id, status => "Succeeded");
    print STDERR "Created normal rnaseq build\n";
    my $tumor = Genome::Test::Factory::Model::RnaSeq->setup_object(
                                                                   processing_profile_id => $pp->id, 
                                                                   subject => $normal->subject,
                                                                   reference_sequence_build => $reference_build,
                                                                  );
    my $tumor_build = Genome::Test::Factory::Build->setup_object(model_id => $tumor->id, status => "Succeeded");

    print STDERR "Created tumor rnaseq build\n";
    my $condition_model_ids_string = $normal->id." ".$tumor->id;

    my $de_model = Genome::Test::Factory::Model::DifferentialExpression->setup_object(
        condition_labels_string => $condition_labels_string,
        condition_model_ids_string => $condition_model_ids_string,
        reference_sequence_build => $reference_build,
        annotation_build => $annot_build,
        %params,
    );
    print STDERR "Created de model\n";
    return $de_model;
}

1;
