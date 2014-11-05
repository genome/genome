package Genome::Model::Tools::MetagenomicClassifier::Rdp::ClassifierBase;

use strict;
use warnings;

use Genome::InlineConfig;

class Genome::Model::Tools::MetagenomicClassifier::Rdp::ClassifierBase {
    is_abstract => 1,
    has => [
        training_set => {
            is => 'Text',
            valid_values => [ valid_training_sets() ],
            doc => 'Training set to use.'
        },
    ],
    has_optional => [
        _factory => { is_transient => 1, },
        _classifier => { is_transient => 1, },
    ],
};

$ENV{PERL_INLINE_JAVA_JNI} = 1;

sub valid_training_sets {
    return (qw/ 4 6 9 10 broad /);
}

sub get_training_path {
    my $class = shift;
    my $training_set = shift;

    $training_set |= '';
    return "/gsc/scripts/share/rdp/$training_set";
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $classifier_properties_path = eval{ $self->classifier_properties_path_for_set( $self->training_set ); };
    if ( not $classifier_properties_path ) {
        $self->error_message($@);
        return;
    }

    my $factory = new FactoryInstance($classifier_properties_path);
    $self->_factory($factory);
    $self->_classifier( $factory->createClassifier );

    return $self;
}

sub new {
    my $class = shift;
    warn "Called method 'new' instantiate $class, but it is deprecated. Please use 'create.'";
    return $class->create(@_);
}

sub base_training_path {
    return "/gsc/scripts/share/rdp/";
}

sub training_path_for_set {
    my ($self, $training_set) = @_;

    die 'No training set given to get training path!' if not $training_set;
    die 'Invalid training set given to get training path!' if not grep { $training_set eq $_ } $self->valid_training_sets;

    my $training_path = $self->base_training_path.'/'.$training_set;
    if ( not -d $training_path ) {
        $self->error_message("Training path does not exist: $training_path");
        return;
    }

    return $training_path;
}

sub classifier_properties_path_for_set {
    my ($self, $training_set) = @_;

    my $training_path = $self->training_path_for_set($training_set);
    return if not $training_path;

    my $classifier_properties_path = $training_path.'/rRNAClassifier.properties';
    if ( not -s $classifier_properties_path ) {
        $self->error_message('No rRNAClassifier.properties in training path! '.$training_path);
        return;
    }

    return $classifier_properties_path;
}

sub classify {
    my ($self, $seq) = @_;

    return if not $self->verify_seq($seq);

    my $parsed_seq = eval{ $self->create_classifier_sequence($seq); };
    if ( not $parsed_seq ) {
        $self->error_message("Can't classify sequence (".$seq->{id}."). Can't create rdp parsed sequence.");
        return;
    }

    # Try to classify 2X - per kathie 2009mar3
    my $classification_result = eval{ $self->_classifier->classify($parsed_seq); };
    if ( not $classification_result ) {
        $classification_result = eval{ $self->_classifier->classify($parsed_seq); };
        if ( not $classification_result ) {
            $self->error_message("Can't classify sequence (".$seq->{seq}."). No classification result was returned from the classifier.");
            return;
        }
    }

    my $complemented = $self->_is_reversed(
        parsed_seq => $parsed_seq,
        classification_result => $classification_result,
    );

    my @assignments = @{$classification_result->getAssignments()->toArray()};
    my %taxa = (
        id => $seq->{id},
        complemented => $complemented,
        classifier => 'rdp',
        root => {
            id => 'Root',
            confidence => $assignments[0]->getConfidence,
        },
    );
    for my $assignment ( @assignments[1..$#assignments] ) {
        # print Dumper([map{ $_->getName } @{$assignment->getClass->getMethods}]);
        # Methods are: getConfidence, getTaxid, getName, getRank	
        my $id = $assignment->getName;
        $id =~ s/\s+/_/g;
        $id =~ s/['"]//g;
        $taxa{ $assignment->getRank || 'root' } = {
            id => $id,
            confidence => $assignment->getConfidence,
        };
    }

    return \%taxa;
}

sub create_classifier_sequence {
    return new edu::msu::cme::rdp::classifier::readseqwrapper::ParsedSequence($_[1]->{id}, $_[1]->{seq});
}

sub verify_seq {
    my ($self, $seq) = @_;

    if ( not $seq ) {
        Carp::confess("No sequence given to classify");
    }

    if ( not $seq->{id} ) {
        Carp::confess('Seq does not have an id: '.Data::Dumper::Dumper($seq));
    }

    if ( not $seq->{seq} or length $seq->{seq} < 50) {
        return;
    }

    return 1;
}

1;

