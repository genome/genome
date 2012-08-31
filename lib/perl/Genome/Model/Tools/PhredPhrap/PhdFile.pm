package Genome::Model::Tools::PhredPhrap::PhdFile;

use strict;
use warnings;

use Genome;

require Cwd;
use Data::Dumper;
require IO::File;

class Genome::Model::Tools::PhredPhrap::PhdFile {
    is => 'Genome::Model::Tools::PhredPhrap::Phds',
    has => [
    phd_file => {
        is => 'String', #file_r
        doc => "File of PHDs",
    },
    ],
};

sub help_brief {
    return 'Phrap starting with PHDs in a file';
}

sub help_detail {
    return $_->[0]->help_brief;
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    $self->phd_file( Cwd::abs_path( $self->phd_file ) );
    
    return $self;
}

sub _verify_phds {
    my $self = shift;

    my $fh = IO::File->new('<' . $self->phd_file)
        or $self->fatal_msg( sprintf('Can\'t open include PHD file (%s): %s', $self->phd_file, $!) );
    chdir $self->phd_dir
        or $self->fatal_msg('Can\'t access directory (%s): %s', $self->phd_dir, $!);

    while ( my $phd_name = $fh->getline ) {
        chomp $phd_name;

        # Verify
        $self->fatal_msg( 
            sprintf('Can\'t find phd (%s) in directory (%s)', $phd_name, $self->phd_dir)
        ) unless -s $phd_name;
    }

    $fh->close;

    chdir $self->_cwd;
    
    return $self->phd_file;
}

1;

#$HeadURL$
#$Id$
