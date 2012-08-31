package Genome::Model::Tools::454::Newbler::SetRef;

use strict;
use warnings;

class Genome::Model::Tools::454::Newbler::SetRef {
    is => 'Genome::Model::Tools::454::Newbler',
    has => [
            dir => {
                    is => 'String',
                    doc => 'pathname of the output directory for project',
                },
            reference_fasta_files => {
                                      is => 'array',
                                      doc => 'file path of the reference sequence(s)',
                                  },
        ],

};

sub help_brief {
"genome-model tools newbler set-ref --dir=DIR --reference-fasta-files='FileA FileB'";
}

sub help_detail {
    return <<"EOS"

EOS
}

sub execute {
    my $self = shift;

    $DB::single = $DB::stopper;
    my $cmd = $self->full_bin_path('setRef') .' '. $self->dir .' '. join(' ',@{$self->reference_fasta_files});
    my $rv = system($cmd);
    unless ($rv == 0) {
        $self->error_message("non-zero return status from command '$cmd'");
        return
    }
    return 1;
}

1;

