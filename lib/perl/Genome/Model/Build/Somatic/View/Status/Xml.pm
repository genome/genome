package Genome::Model::Build::Somatic::View::Status::Xml;

use strict;
use warnings;

use Genome;


class Genome::Model::Build::Somatic::View::Status::Xml {
    is => 'Genome::Model::Build::View::Status::Xml',
};

#Separate method to make it easy to override in subclasses
sub get_summary_report_location {
    my $self = shift;
    my $build = $self->subject;
    
    #A default value equivalent to what was previously hardcoded in the XSL.
    return $build->data_directory . '/cancer_report.html';
}

1;
