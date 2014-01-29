package Genome::Model::Tools::WuBlast::Xdformat;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::WuBlast::Xdformat {
    is => 'Command',
    has => [
            database => {
                         is => 'String',
                         is_input => 1,
                         doc => 'The path to a xdformat database',
                     },
    ],
    has_optional => [
                     db_type => {
                                 is => 'String',
                                 is_input => 1,
                                 default => 'n',
                                 doc => 'Database type (n)compare nucleotide sequence, (p) compare protein sequence, (x) compare peptide sequence queries to nucleotide sequence databases dynamically translated in all 6 reading frames',
                             },
    ],
};

sub help_brief {
    "Wrapper for running xdformat",
}

sub execute {
    my $self = shift;

    for my $method ( $self->_pre_execute_methods  ) {
        $self->$method
            or return;
    }

    my $cmd = sprintf(
        'xdformat -%s -%s %s %s',
        $self->db_type,
        $self->_operation_character,
        $self->database,
        ( $self->can('fasta_files') ? join(' ', $self->fasta_files) : '' ),
    );
    $self->debug_message('Running: '.$cmd);

    my $rv = system($cmd);
    unless ( $rv == 0 ) {
        $self->error_message( 
            sprintf('Could not %s database (%s).  Return value: %s.', $self->_operation_name, $self->database, $rv) 
        );
        return;
    }

    for my $method ( $self->_post_execute_methods  ) {
        $self->$method
            or return;
    }
   
    return 1;
}

#< Pre and Post Execute Methods >#
sub _pre_execute_methods {
    return;
}

sub _post_execute_methods {
    return;
}

1;

#$HeadURL$
#$Id$
