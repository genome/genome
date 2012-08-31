package Genome::Model::Tools::Fasta::Trim;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fasta::Trim {
    is => 'Genome::Model::Tools::Fasta',
    is_abstract => 1,
    has_optional => [
    output_fasta_file => {
        type => 'Text',
        doc => 'Output file.  Defaults to inserting the command name between the fasta file base name and the fasta extension.  The resulting qual file will be named as above, but with a ".qual" extension',
    },
    ],
};

#< Command Interface >#
sub help_brief {
    my $class = shift;

    return sprintf(
        "Trim vector and/or low quality from FASTA and Quality files using %s",
        ( $class eq __PACKAGE__ ? '' : 'using '.$class->executable ),
    );
}

sub help_detail { 
    return <<EOS 
    Requires a quality file.
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    unless ( $self->have_qual_file ) {
        $self->error_message("A qual file (fasta file name w/ .qual) is required");
        $self->delete;
        return;
    }

    unless ( $self->output_fasta_file ) {
        $self->output_fasta_file( $self->_default_output_fasta_file );
        unlink $self->output_fasta_file if -e $self->output_fasta_file;
        unlink $self->output_qual_file if -e $self->output_qual_file;
    }
    
    return $self;
}

#< Output Fasta and Qual >#
sub _default_output_fasta_file {
    return $_[0]->fasta_file_with_new_suffix( $_[0]->executable );
}

sub output_qual_file {
    return $_[0]->output_fasta_file.'.qual';
}

1;

#$HeadURL$
#$Id$
