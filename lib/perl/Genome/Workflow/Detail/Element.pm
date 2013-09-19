package Genome::Workflow::Detail::Element;

use strict;
use warnings;

use Genome;


class Genome::Workflow::Detail::Element {
    is_abstract => 1,
};


# ------------------------------------------------------------------------------
# Abstract methods
# ------------------------------------------------------------------------------

sub get_xml {}
sub get_xml_element {}
sub validate {}


1;
