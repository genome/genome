package Genome::Model::Tools::Fasta::Sanitize;

use strict;
use warnings;

use Genome;

require Bio::SearchIO;
require Bio::Seq;
require Bio::SeqIO;
require Bio::Seq::PrimaryQual;
use Data::Dumper;
require File::Temp;
require Genome::Model::Tools::WuBlast::Blastn;
require Genome::Model::Tools::WuBlast::Xdformat::Create;
require IO::File;
use Regexp::Common;

my @SANITATION_PROPERTIES = ( # used arrow format for convenience
    min_threshold => {
        is => 'Integer',
        doc => 'Sequecne minimum threshold, includes with an length greater than or equal to this value',
    },
    max_threshold => {
        is => 'Integer',
        doc => 'Contig maximum threshold, includes contigs with an unpadded length less than or equal to this value',
    }, 
    replace_xs_with_ns => {
        is => 'Boolean',
        doc => 'Replace all Xs (case insensitive) with Ns.  This is an alias to "replace" with parameter "x:N"',
    },
    replace => {
        is => 'Text',
        doc => 'Replace a case insensitive string (or regular expression) of a characters with a new string.  Separate the match and new string by a colon. From the command line, enclose the replace string w/in single quotes (as shown in the examples) to ensure that it is not interperetted by the operating system.  Examples: Replace 15 As with 10 Ns => \'AAAAAAAAAAAAAAA:NNNNNNNNNN\', or as a reg exp => \'A{15}:NNNNNNNNNN\'; Replace 15 or more As w/ Ns (must use reg exp) => \'A{15,}:NNNNNNNNNN\'',
    },
    capitalize => {
        is => 'Boolean',
        doc => 'Captilize all characters',
    },
    rm_descs => {
        is => 'Boolean',
        doc => 'Removes the description in each fasta header',
    },
);
class Genome::Model::Tools::Fasta::Sanitize {
    is  => 'Genome::Model::Tools::Fasta',
    has_optional => [
    sanitized_fasta_file => {
        is => 'Text',
        doc => 'The name of the sanitized fasta file (output).  Default will have ".sanitized" between the file base and the fasta extension',
    },
    @SANITATION_PROPERTIES
    ],
};

sub sanitation_methods {
    return grep { not ref } @SANITATION_PROPERTIES;
}

sub default_sanitized_fasta_file {
    die "'default_sanitized_fasta_file' is an object method, but called was called as class method\n" unless ref $_[0];
    return $_[0]->fasta_file_with_new_suffix('sanitized');
}

sub help_brief {
    return 'Cleans FASTA (and Quality) files';
}

sub help_detail { 
    return help_brief();
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    # Make sure we're sanitizing...
    if ( not grep { defined $self->$_ } sanitation_methods() ) {
        $self->error_message("No sanitation methods indicated");
        return;
    }

    # Verify the replace param 
    if ( $self->replace ) {
        my ($find, $replace) = split(':', $self->replace);
        unless ( defined $find and $find ne '' ) {
            $self->error_message("No find string (regexp) for param 'replace'.  Use format 'FIND:REPLACE'");
            return;
        }
        unless ( defined $replace and $replace =~ /^\w+$/i ) { # allow replacing w/ nothing?
            $self->error_message("No or invalid replace string for param 'replace'.  Use format 'FIND:REPLACE'");
            return;
        }
    }

    return $self;
}

sub execute {
    my $self = shift;

    #< In ># - create own parser?? otherwise can't remove illegal characters
    my $bioseq_in = $self->get_fasta_reader( $self->fasta_file )
        or return;
    
    #< OUT >#
    unless ( $self->sanitized_fasta_file ) {
        $self->sanitized_fasta_file( $self->default_sanitized_fasta_file );
        #unlink??
    }
    my $bioseq_out = $self->get_fasta_writer( $self->sanitized_fasta_file )
        or return;

    #< Sanitize >#
    my @methods = map { '_'.$_ } grep { defined $self->$_ } sanitation_methods();
    SEQ: while ( my $bioseq = $bioseq_in->next_seq ) {
        for my $method ( @methods ) {
            next SEQ unless $self->$method($bioseq);
        }
        $bioseq_out->write_seq($bioseq);
    }

    $self->debug_message('Wrote sanitized sequences to '.$self->sanitized_fasta_file);

    return 1;
}

sub _max_threshold {
    my ($self, $bioseq) = @_;

    return $bioseq->length <= $self->max_threshold;
}

sub _min_threshold {
    my ($self, $bioseq) = @_;

    return $bioseq->length >= $self->min_threshold;
}

sub _capitalize {
    my ($self, $bioseq) = @_;

    return $bioseq->seq( uc $bioseq->seq );
}

sub _find_and_replace {
    my ($self, $bioseq, $find, $replace) = @_;

    my $seq = $bioseq->seq;
    $seq =~ s/$find/$replace/ig;
    $bioseq->seq($seq);

    return 1;
}

sub _replace {
    my ($self, $bioseq) = @_;

    return $self->_find_and_replace( $bioseq, split(':', $self->replace) );
}

sub _replace_xs_with_ns {
    my ($self, $bioseq) = @_;

    return $self->_find_and_replace($bioseq, 'x', 'N');
}

sub _rm_descs {
    my ($self, $bioseq) = @_;

    $bioseq->desc(undef);
    
    return 1;
}

1;

#$HeadURL$
#$Id$
