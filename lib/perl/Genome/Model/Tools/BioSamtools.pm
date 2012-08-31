package Genome::Model::Tools::BioSamtools;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools {
    is => ['Command'],
    is_abstract => 1,
};

sub help_detail {
    "These commands are setup to run perl5.10.0 scripts that use Bio-Samtools and require bioperl v1.6.0.  Most require 64-bit architecture except those that simply work with output files from other Bio-Samtools commands.";
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    unless ($] > 5.010) {
        die 'Bio::DB::Sam requires perl 5.10 or greater! Consider using \'/usr/bin/perl `which gmt`\' instead of \'gmt\' for all bio-samtools commands!';
    }
    require Bio::DB::Sam;
    return $self;
}

sub perl_path {
    die('Remove perl-5.10.0');
}

sub bin_path {
    die('Remove /gsc/var/tmp/Bio-SamTools');
}

sub execute_path {
    die('Do not call execute_path in '. __PACKAGE__);
}

sub bioperl_path {
    die('Do not define bioperl version');
}

1;
