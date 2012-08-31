package Genome::Model::Tools::454::Seqclean;

use strict;
use warnings;

use Genome;

use File::Basename;
use Cwd;

class Genome::Model::Tools::454::Seqclean {
    is => ['Command'],
    has => [
            in_fasta_file => {
                              is_input => 1,
                              is => 'string',
                              doc => 'a fasta file path to run seqclean on',
                          },
        ],
    has_optional => [
                     seqclean_params => {
                                         is_param => 1,
                                         is => 'string',
                                         doc => 'a set of params to use with seqclean, run \'seqclean --help\' for available options',
                                     },
                     adaptor_sequence_db => {
                                             is_input => 1,
                                             is => 'string',
                                             doc => 'a path to the blastable database of adaptor sequences, equivalent to -v seqclean option',
                                         },
                     contaminant_db => {
                                        is_input => 1,
                                        is => 'string',
                                        doc => 'a path to the blastable database of contaminant sequences, equivalent to -s seqclean option',
                                    },
                     out_fasta_file => {
                                        is_output => 1,
                                        is => 'string',
                                        doc => 'a file path for the fasta seqcleaned reads, equivalent to -o seqclean option',
                                    },
                     seqclean_report => {
                                         is_output => 1,
                                         is => 'string',
                                         doc => 'a file path for the seqclean report output, equivalent to -r seqclean option',
                                     },
                 ],
};

sub help_brief {
    "a genome-model tool to run seqclean",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt 454 seqclean ...
EOS
}

sub help_detail {
    return <<EOS
This tool runs the seqclean command on the given fasta file.
The default seqclean report goes to <in fasta file>.cln
and the default fasta output goes to <in fasta file>.clean.
The blastable databases must be formatted first by 'formatdb'.
'seqclean' creates directories and several supporting files.
All 'seqclean' output is dumped to the output fasta directory.
For additional parameters and options see 'seqclean --help'.
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    unless (-e $self->in_fasta_file) {
        die 'Input fasta file '. $self->in_fasta_file .' does not exist';
    }

    return $self;
}

sub execute {
    my $self = shift;

    my $cwd = getcwd;

    my $dirname;
    if ($self->out_fasta_file) {
        $dirname = dirname($self->out_fasta_file);
    } else {
        $dirname = dirname($self->in_fasta_file);
    }
    chdir ($dirname) || die ("Failed to change directory to '$dirname':  $!");

    my $cmd = 'seqclean '. $self->in_fasta_file;
    $cmd .= ' '. $self->seqclean_params if $self->seqclean_params;
    $cmd .= ' -o '. $self->out_fasta_file if $self->out_fasta_file;
    $cmd .= ' -r '. $self->seqclean_report if $self->seqclean_report;
    $cmd .= ' -v '. $self->adaptor_sequence_db if $self->adaptor_sequence_db;
    $cmd .= ' -s '. $self->contaminant_db if $self->contaminant_db;

    $self->status_message('Running: '. $cmd);

    my $rv = system($cmd);

    chdir ($cwd) || die ("Failed to change directory to '$cwd':  $!");;

    unless ($rv == 0) {
        $self->error_message("non-zero exit code($rv) returned from command '$cmd'");
        return;
    }

    return 1;
}


1;

