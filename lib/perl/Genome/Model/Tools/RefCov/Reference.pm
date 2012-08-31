package Genome::Model::Tools::RefCov::Reference;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RefCov::Reference {
    has => [
        coverage => {},
        id => {},
        fasta => {
            is => 'String',
            doc => 'This should be deprecated for passing the fai directly',
        },
        chromosome => {},
        start => {},
        stop => {},
    ],
    has_optional => {
        fai => {},
        reflen => {
            is_calculated => 1,
            calculate_from => ['start','stop'],
            calculate => sub {
                my ($start,$stop) = @_;
                return (($stop - $start ) + 1);
            },
        },
        covlen => {
            is_calculated => 1,
            calculate_from => ['coverage'],
            calculate => sub {
                my $coverage = shift;
                return scalar(@{$coverage});
            },
        },
        sequence => {
            is => 'String',
        },
    }
};

sub create {
    my $class = shift;
    my %params = @_;
    unless ($] > 5.010) {
        die "Bio::DB::Sam requires perl 5.10 or greater!";
    }
    require Bio::DB::Sam;
    my $coverage = delete($params{coverage});
    my $self = $class->SUPER::create(%params);
    $self->coverage($coverage);
    
    # TODO: We shouldn't have to load the FAI file, the loaded object should be passed in
    # This could cause a lot more IO than is really necessary loading the file every time
    my $fai = Bio::DB::Sam::Fai->load( $self->fasta );
    unless ($fai) {
        die('Failed to load FASTA index for FASTA: '. $self->fasta);
    }
    $self->fai($fai);

    if ($self->covlen != $self->reflen) {
        die(__PACKAGE__ .' coverage length '. $self->covlen .' does not match reference length '. $self->reflen .'. Investigate.');
    }

    # Load the sequence string.
    my $fetch_string   = $self->chromosome() .':'. $self->start() .'-'. $self->end();
    $self->sequence($self->fai->fetch( $fetch_string ));

    return $self;
}

1;  # End of package.
