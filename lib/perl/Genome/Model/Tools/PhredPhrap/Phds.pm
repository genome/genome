package Genome::Model::Tools::PhredPhrap::Phds;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::PhredPhrap::Phds {
    is => 'Genome::Model::Tools::PhredPhrap::Base',
    has => [
    ],
};

require Genome::Model::Tools::PhredPhrap::PhdToFasta;
use Data::Dumper;

sub help_brief {
    return 'Phrap starting with PHDs in a project\'s phd_dir';
}

sub help_detail {
    return '';
}

sub _files_to_remove {
    return (qw/ default_fasta_file default_qual_file /);
}

sub _handle_input {
    my $self = shift;

    $self->status_message("Verifying PHDs");
    my $phd_file = $self->_verify_phds;

    $self->status_message("PHD to FASTA and Quality");
    $self->_phd2fnq($phd_file);

    return 1;
}

sub _verify_phds {
    my $self = shift;

    my $phd_dir = $self->_directory->phd_dir;
    my $dh = IO::Dir->new($phd_dir)
        or $self->fatal_msg( sprintf('Can\'t open dir (%s)', $phd_dir) );
    my $phd_file = $self->default_phd_file;
    unlink $phd_file if -e $phd_file;
    my $phd_fh = IO::File->new("> $phd_file")
        or $self->fatal_msg("Can\'t open phd file ($phd_file) for writing");

    while ( my $phd_name = $dh->read ) {
        next unless $phd_name =~ m#\.phd\.\d+$#;
        # TODO Exclude
        $phd_fh->print("$phd_name\n");
    }

    $dh->close;
    $phd_fh->close;

    $self->fatal_msg("No phds found in directory ($phd_dir)") unless -s $phd_file;

    return $phd_file;
}

sub _phd2fnq {
    my ($self, $phd_file) = @_;

    return Genome::Model::Tools::PhredPhrap::PhdToFasta->execute(
        phd_file => $phd_file,
        phd_dir => $self->_directory->phd_dir,
        fasta_file => $self->default_fasta_file,
    );
}

1;

#$HeadURL$
#$Id$
