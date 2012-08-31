# review gsanders
# Is there a reason this is not using the UR framework for standard stuff like accessors/mutators and the like?

package Genome::Utility::MetagenomicClassifier::ChimeraClassifier;

use strict;
use warnings;

require Bio::Seq;
require Bio::Taxon;
require Genome::Utility::MetagenomicClassifier::ChimeraClassification;

sub create {
    my ($class, %params) = @_;
    
    my $self = bless \%params, $class;

    my $classifier = $self->classifier;
    unless ($classifier) {
        if ($self->{training_set}) {
            my $training_set = $self->{training_set};
            $classifier = new Genome::Utility::MetagenomicClassifier::Rdp::Version2x1(training_set => $training_set);
        }
        elsif ($self->{training_path}) {
            $classifier = new Genome::Utility::MetagenomicClassifier::Rdp::Version2x1(training_path => $self->{training_path});
        }
        else {
            $classifier = new Genome::Utility::MetagenomicClassifier::Rdp::Version2x1();
        }
        $self->classifier($classifier);
    }

    my $probe_pairs = $self->probe_pairs;
    unless ($probe_pairs) {
        $self->probe_pairs([[0,400],[400,400],[800,400], [-400,400]]);
    }

    return $self;
}

sub probe_pairs {
    my $self = shift;
    my $probe_pairs = shift;
    if ($probe_pairs) {
        $self->{probe_pairs} = $probe_pairs;
    }
    else {
        $probe_pairs = $self->{probe_pairs};
    }
    return $probe_pairs;
}
sub classifier {
    my $self = shift;
    my $classifier = shift;
    if ($classifier) {
        $self->{classifier} = $classifier;
    }
    else {
        $classifier = $self->{classifier};
    }
    return $classifier;
}

sub classify {
    my ($self, $seq) = @_;

    unless ( $seq ) {
        warn ("No sequence to classify");
        #$self->error_message("No sequence to classify");
        return;
    }

    my $sequence = $seq->seq;
#    if (length $sequence < 1000) {
#        print "sequence is too short\n";
#        return;
#    }
    my $classification = $self->classifier->classify($seq);

    unless ($classification) {
        warn "Failed to classify ".$seq->display_name. " with rdp";
        return;
    }

    my @probe_sequences = $self->generate_probe_sequences($seq);


    my @probe_classifications;
    foreach my $probe (@probe_sequences) {
        my $probe_classification = $self->classifier->classify($probe);
        if ($probe_classification) {
            push @probe_classifications, $probe_classification;
        }
    }
    return Genome::Utility::MetagenomicClassifier::ChimeraClassification->new(
        name => $seq->display_name,
        classification => $classification,
        probe_classifications => \@probe_classifications,
    );
}

sub generate_probe_sequences {
    my $self = shift;
    my $seq = shift;

    my $seq_name = $seq->display_name;
    my $sequence = $seq->seq;
    my $full_length = $seq->length;

    my @probe_pairs = @{$self->probe_pairs};

    my @probes;
    my $probe_count = 1;
    foreach my $probe_pair (@probe_pairs) {
        my ($start, $length) = @$probe_pair;
        if ($start > $full_length) {
            warn "$seq_name is too short";
            next;
        }

        if ($start < 0) {
            $start = $full_length + $start;
        }
        my $stop = $start + $length;
        my $probe_sequence = substr($sequence, $start, $length);
        my $probe_name = $seq_name."-$start-$stop";
        my $probe_seq = Bio::Seq->new(
            -id => $probe_name,
            -seq => $probe_sequence,
        );
        push @probes, $probe_seq;
    }
    return @probes;
}

1;

#$HeadURL: $
#$Id: $
