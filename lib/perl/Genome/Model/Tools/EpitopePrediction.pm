package Genome::Model::Tools::EpitopePrediction;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::EpitopePrediction {
    is  => 'Command::Tree',
    doc => 'Pipeline for Immune Epitope Prediction',
};

sub help_brief {
    "Different pipeline steps for Immune Epitope Prediction for MHC Class-I binding",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt small-rna ...    
EOS
}
sub help_detail {
    "These commands are setup to run perl5.12.1";
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    unless ($] > 5.012) {
        die 'Bio::DB::Sam requires perl 5.12! Consider using gmt5.12.1 instead of gmt for all epitope-prediction commands!';
    }
    require Bio::DB::Sam;
    return $self;
}



1;