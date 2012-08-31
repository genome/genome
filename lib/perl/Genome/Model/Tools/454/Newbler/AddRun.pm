package Genome::Model::Tools::454::Newbler::AddRun;

use strict;
use warnings;

class Genome::Model::Tools::454::Newbler::AddRun {
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
            is_paired_end => {
                              is => 'boolean',
                              doc => 'a flag for paired end status',
                              default_value => 0,
                          },
        ],

};


sub help_brief {
    "add a 454 region to a newbler project"
}

sub help_detail {
    return <<"EOS"

EOS
}

sub execute {
    my $self = shift;

    $DB::single = $DB::stopper;
    my $options = $self->is_paired_end ? ' -p ' : ' ';
    my $cmd = $self->full_bin_path('addRun') . $options . $self->dir .' '. join(' ',@{$self->runs});

    my $rv = system($cmd);
    unless ($rv == 0) {
        $self->error_message("non-zero return status from command '$cmd'");
        return;
    }
    return 1;
}

1;

