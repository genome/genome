package Genome::Model::ClinSeq::Command::Converge::Snvs;
use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::Converge::Snvs {
    is => 'Genome::Model::ClinSeq::Command::Converge::Base',
    doc => 'converge SNV results from mutiple clinseq builds'
};

sub resolve_labels {
    my $self = shift;
  
    # tweak what gets passed in
    
    my %results = $self->SUPER::resolve_labels(@_);
    
    # tweak what came out

    return %results;
}

sub execute {
    my $self = shift;

    my @stuff;
    my @labels;
    my %labeled_stuff = $self->resolve_labels(\@stuff,\@labels);   

    return 1;
};

1;

__END__
@ISA = ('Genome::Model::ClinSeq::Command::Converge::Base');

