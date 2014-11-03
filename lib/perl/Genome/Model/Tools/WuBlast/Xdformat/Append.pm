package Genome::Model::Tools::WuBlast::Xdformat::Append;

use strict;
use warnings;

use Genome;

require Genome::Model::Tools::WuBlast::Xdformat::Verify;
    
class Genome::Model::Tools::WuBlast::Xdformat::Append {
    is => 'Genome::Model::Tools::WuBlast::Xdformat',
    has => [
    database => {
        is => 'String',
        is_input => 1,
        doc => 'the path to a new or existing database',
    },
    ],
    has_many => [
    fasta_files => {
        is => 'String',
        is_input => 1,
        doc => 'a list of paths to fasta sequence files',
    },
    ],
};

#< Standard command methods >#
sub help_brief {
    return "Appends to an xdformat database";
}

sub help_detail {
    return help_brief();
    return <<EOS
EOS
}

#< Pre and Post Execute Methods >#
sub _pre_execute_methods {
    return (qw/ _verify_db /);
}

sub _verify_db {
    my $self = shift;

    return Genome::Model::Tools::WuBlast::Xdformat::Verify->execute(
        database => $self->database,
        db_type => $self->db_type,
    )->result;
}

#< Operation >#
sub _operation_name {
    return 'append';
}

sub _operation_character {
    return 'a';
}

1;

#$HeadURL$
#$Id$
