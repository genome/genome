package Genome::Model::Tools::PhredPhrap::ScfFile;

use strict;
use warnings;

use Genome;

use Data::Dumper;
require IO::File;

class Genome::Model::Tools::PhredPhrap::ScfFile {
    is => 'Genome::Model::Tools::PhredPhrap::Scfs',
    has => [
    scf_file => {
        is => 'String', #file_r
        is_optional => 0,
        doc => "File of SCFs",
    },
    ],
};

sub help_brief {
    return 'Phrap starting with a file of SCFs';
}

sub help_detail {
    return '';
}

sub _verify_scfs {
    my $self = shift;

    my $fh = IO::File->new('<' . $self->scf_file);
    unless ( $fh ) {
        $self->error_message( sprintf('Can\'t open include scf file (%s) for reading', $self->scf_file) );
        return;
    }
    
    unless ( chdir $self->_directory->chromat_dir ) {
        $self->error_message('Can\'t access directory (%s): %s', $self->_directory->chromat_dir, $!);
        return;
    }

    while ( my $scf_name = $fh->getline ) {
        chomp $scf_name;
        $scf_name =~ s/\s+//g;
        $scf_name =~ s/\.gz//;
        $scf_name =~ s/\.exp//;
        $scf_name =~ s/\.phd\.\d+//;
        next if $scf_name eq '';

        # Verify
        unless ( -s $scf_name or -s "$scf_name.gz" ) {
            $self->error_message( 
                sprintf('Can\'t find scf (%s) in directory (%s)', $scf_name, $self->_directory->chromat_dir)
            );
            return;
        }
    }

    chdir $self->_cwd;
    $fh->close;

    return $self->scf_file;
}

1;

#$HeadURL$
#$Id$
