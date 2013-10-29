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

    my $fh = IO::File->new('<' . $self->phd_file);
    if ( not $fh ) {
        $self->error_message("$!") if $!;
        $self->error_message('Can\'t open PHD file! '.$self->phd_file);
        return;
    }

    my $chdir_ok = chdir $self->phd_dir;
    if ( not $chdir_ok ) {
        $self->error_message("$!") if $!;
        $self->error_message('Can\'t access directory! '.$self->phd_dir);
        return
    }

    while ( my $phd_name = $fh->getline ) {
        chomp $phd_name;

        # Verify
        if ( not -s $phd_name ) {
            $self->error_message( sprintf('Can\'t find phd (%s) in directory (%s)', $phd_name, $self->phd_dir) );
            return;
        }
    }

    $fh->close;

    chdir $self->_cwd;

    return $self->phd_file;
}

1;

