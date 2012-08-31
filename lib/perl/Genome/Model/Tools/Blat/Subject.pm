package Genome::Model::Tools::Blat::Subject;

use strict;
use warnings;

use Genome;
use Workflow;
use File::Basename;

class Genome::Model::Tools::Blat::Subject {
    is => 'Command',
    has => [
        query_file => {
            doc => 'the query file to blat against subject(either .fa, .nib, or .2bit)',
            is => 'String',
            is_input => 1,
        },
        subject_file => {
            doc => 'the subject file to align against(either .fa, .nib, or .2bit)',
            is => 'String',
            is_input => 1,
        },
    ],
    has_optional => [
        output_directory => {
            doc => 'The directory where output alignment files are written',
            is_input => 1,
        },
        blat_params => {
            doc => 'additional command line parameters passed to blat execution',
            is => 'String',
            is_param => 1,
            is_input => 1,
        },
        alignment_file => {
            doc => 'the output file to store blat alignments',
            is => 'String',
            is_output => 1,
        },
        aligner_output_file => {
            doc => 'a file path to store blat aligner output',
            is => 'String',
            is_output => 1,
        },
    ],
};

sub help_brief {
    "run one instance of a blat alignment",
}

sub help_synopsis {
    return <<"EOS"
gmt blat subject --subject_file --query_file {--alignment_file, --blat_params]
EOS
}

sub help_detail {
    return <<EOS
This tool runs a single instance of a blat alignment using
one subject/database file(--subject-file) and one query file(--query-file).
The subject ans query files must be in .fa, .nib, or .2bit format.
The output of the blat alignment is stored in the supplied file path(--alignment-file).
Any command line parameters for blat should be supplied with --blat-params.
For blat parameters see the blat help/man documentation.
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    my @suffix_list = (qr{\.fa.*}, qr{\.nib}, qr{\.2bit});
    my ($query_basename,$query_directory,$query_suffix) = fileparse($self->query_file,@suffix_list);
    my ($subject_basename,$subject_directory,$subject_suffix) = fileparse($self->subject_file,@suffix_list);

    unless ($query_suffix =~ /^.(fa.*|nib|2bit)/) {
        $self->error_message("Expected a .fa, .nib, or .2bit format file and found '$query_suffix' for query file.");
        die;
    }
    unless ($subject_suffix =~ /^.(fa|nib|2bit)/) {
        $self->error_message("Expected a .fa, .nib, or .2bit format file and found '$subject_suffix' for subject file.");
        die;
    }
    my $output_directory;
    if ($self->output_directory) {
        $output_directory = $self->output_directory;
    } else {
        $output_directory = $query_directory;
    }
    unless ($self->alignment_file) {
        $self->alignment_file($output_directory .'/'. $query_basename .'_'. $subject_basename .'.psl');
    }
    unless ($self->aligner_output_file) {
        $self->aligner_output_file($output_directory .'/'. $query_basename .'_'. $subject_basename .'.out');
    }
    return $self;
}

sub execute {
    my $self = shift;

    my $blat_param_string = $self->blat_params || '';

    my $cmd = 'blat '. $self->subject_file .' '. $self->query_file .' '. $blat_param_string .' '. $self->alignment_file .' > '. $self->aligner_output_file .' 2>&1';
    $self->status_message('Running: '. $cmd);
    my $rv = system($cmd);
    unless ($rv == 0) {
        $self->error_message("non-zero return value($rv) from command $cmd");
        return;
    }
    return 1;
}


1;
