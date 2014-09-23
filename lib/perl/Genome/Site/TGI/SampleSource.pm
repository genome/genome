package Genome::Site::TGI::SampleSource;

use strict;
use warnings;

use Genome;
require Carp;

class Genome::Site::TGI::SampleSource {
    is_abstract => 1,
    subclassify_by => '_subclass_by_subject_type',
    has_optional => [
        taxon           => { is => 'Genome::Taxon', id_by => 'taxon_id' },
        species_name    => { via => 'taxon' },
    ],
    has_many_optional => [
        samples         => { is => 'Genome::Sample', reverse_id_by => 'source' },
        sample_names    => { via => 'samples', to => 'name' },
    ],
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'Base class for individual and population group',
};

sub _subclass_by_subject_type {
    my ($clas, $sample_source) = @_;

    my $subject_type = $sample_source->subject_type;
    if ( $subject_type eq 'organism individual' or $subject_type eq 'organism_individual' ) {
        return 'Genome::Individual';
    }
    elsif ( $subject_type eq 'population group' or $subject_type eq 'population_group' ) {
        return 'Genome::PopulationGroup';
    }
    else {
        Carp::confess("Unknown subject type ($subject_type). Cannot determine subclass.");
    }
 
}

1;

