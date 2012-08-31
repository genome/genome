package Genome::Model::Tools::454::Newbler::RemoveRun;

use strict;
use warnings;

class Genome::Model::Tools::454::Newbler::RemoveRun {
    is => 'Genome::Model::Tools::454::Newbler',
    has => [
            dir => {
                    is => 'String',
                    doc => 'pathname of the output directory for project',
                },
            runs => {
                       is => 'array',
                       doc => 'a list of sff files, directories, or fasta files',
                   },
        ],

};

sub help_brief {
    "remove a 454 region to a newbler project"
}

sub help_detail {
    return <<"EOS"

EOS
}

sub execute {
    my $self = shift;

    $DB::single = $DB::stopper;
    my $cmd = $self->full_bin_path('removeRun') .' '. $self->dir .' '. join(' ',@{$self->runs});
    my $rv = system($cmd);
    unless ($rv == 0) {
        $self->error_message("non-zero return status from command '$cmd'");
        return
    }
    return 1;
}

1;

