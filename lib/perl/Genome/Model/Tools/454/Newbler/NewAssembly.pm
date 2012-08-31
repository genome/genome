package Genome::Model::Tools::454::Newbler::NewAssembly;

use strict;
use warnings;


class Genome::Model::Tools::454::Newbler::NewAssembly {
    is => 'Genome::Model::Tools::454::Newbler',
    has => [
            dir => {
		    is => 'String',
                    doc => 'pathname of the output directory',
                   },
        ],

};

sub help_brief {
    "create a new newbler project for denovo assembly"
}

sub help_detail {
    return <<"EOS"
This runs the "createProject" newbler command to make a new assembly.
EOS
}

sub execute {
    my $self = shift;

    $DB::single = $DB::stopper;
    my $cmd = $self->full_bin_path('createProject') .' -t asm '. $self->dir;
    my $rv = system($cmd);
    unless ($rv == 0) {
        $self->error_message("non-zero return status from command '$cmd'");
        return
    }
    return 1;
}

1;

