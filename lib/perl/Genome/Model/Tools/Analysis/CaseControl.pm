package Genome::Model::Tools::Analysis::CaseControl;

use strict;
use warnings;

use Genome;

use XML::Simple;
use File::Basename;

class Genome::Model::Tools::Analysis::CaseControl {
    is => ['Genome::Model::Tools::Analysis'],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools for analysis of CaseControl disease datasets.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis mendelian --help ...
EOS
}

sub help_detail {
    return <<EOS
EOS
}

sub newbler_bin {
    my $self = shift;

    my $bin_path = $self->bin_path;

    return $bin_path;
}

sub full_bin_path {
    my $self = shift;
    my $cmd = shift;

    return $self->newbler_bin .'/'. $cmd;
}

sub get_newbler_version_from_xml_file {
    my $class = shift;
    my $xml_file = shift;

    if (ref($class)) {
        $class = ref($class);
    }
    my $xml = XML::Simple->new( Forcearray => [ qw( book ) ], keyattr => { book => 'isbn' } );
    my $xml_data = $xml->XMLin($xml_file);
    return unless exists $xml_data->{ProjectInformation}->{SoftwareVersion};
    return $xml_data->{ProjectInformation}->{SoftwareVersion};
}


1;

