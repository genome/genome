package Genome::Model::Tools::Crossmatch::Run;

use strict;
use warnings;

use Genome;
use Workflow;
use File::Basename;

class Genome::Model::Tools::Crossmatch::Run {
    is => 'Command',
    has => [
            query_file => {
                           doc => 'the query file to cross_match against subject(fasta format)',
                           is => 'String',
                           is_input => 1,
                       },
            subject_file => {
                             doc => 'the subject file to align against(fasta format)',
                             is => 'String',
                             is_input => 1,
                         },
        ],
    has_optional => [
                     cm_params => {
                                   doc => 'additional command line parameters passed to cross_match execution',
                                   is => 'String',
                                   is_input => 1,
                                   is_param => 1,
                               },
                     alignment_file => {
                                        doc => 'the output file to store cross_match alignments',
                                        is => 'String',
                                        is_output => 1,
                                    },
                     aligner_output_file => {
                                             doc => 'the file to store cross_match aligner output',
                                             is => 'String',
                                             is_output => 1,
                                         },
                 ],
};

sub help_brief {
    "run one instance of a cross_match alignment",
}

sub help_synopsis {
    return <<"EOS"
gmt cross-match run --subject_file --query_file {--alignment_file, --cm_params]
EOS
}

sub help_detail {
    return <<EOS
This tool runs a single instance of a cross_match alignment using
one subject/database file(--subject-file) and one query file(--query-file).
The subject and query files must be in fasta format with a suffix like .fa*.
The output of the cross_match alignment is stored in the supplied file path(--alignment-file).
Any command line parameters for cross_match should be supplied with --cm-params.
For cross_match parameters see the cross_match documentation.(http://www.phrap.org/phredphrap/phrap.html)
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    my @suffix_list = (qr{\.fa.*});
    my ($query_basename,$query_directory,$query_suffix) = fileparse($self->query_file,@suffix_list);
    my ($subject_basename,$subject_directory,$subject_suffix) = fileparse($self->subject_file,@suffix_list);

    unless ($query_suffix =~ /^\.fa/) {
        $self->error_message("Expected a fasta format file named like .fa* for query file but found '$query_suffix'");
        die;
    }
    unless ($subject_suffix =~ /^\.fa/) {
        $self->error_message("Expected a fasta format file named like .fa* for subject file but found '$subject_suffix'");
        die;
    }
    unless ($self->alignment_file) {
        $self->alignment_file($query_directory .'/'. $query_basename .'_'. $subject_basename .'.cm');
    }
    unless ($self->aligner_output_file) {
        $self->aligner_output_file($query_directory .'/'. $query_basename .'_'. $subject_basename .'.out');
    }
    return $self;
}

sub execute {
    my $self = shift;

    my $cm_param_string = $self->cm_params || '';

    my $cmd = '(cross_match.test '. $self->subject_file .' '. $self->query_file .' '. $cm_param_string .' > '. $self->alignment_file .') > '. $self->aligner_output_file . ' 2>&1';
    $self->status_message('Running: '. $cmd);
    my $rv = system($cmd);
    unless ($rv == 0) {
        $self->error_message("non-zero return value($rv) from command $cmd");
        return;
    }
    return 1;
}


1;
