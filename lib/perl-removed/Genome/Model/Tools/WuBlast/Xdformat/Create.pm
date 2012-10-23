package Genome::Model::Tools::WuBlast::Xdformat::Create;

use strict;
use warnings;

use Genome;

use Data::Dumper;
require Genome::Model::Tools::WuBlast::Xdformat::Verify;

class Genome::Model::Tools::WuBlast::Xdformat::Create {
    is => 'Genome::Model::Tools::WuBlast::Xdformat',
    has_optional => [
    overwrite_db => {
        is => 'Boolean',
        default => 0,
        doc => 'If existing, remove the database files to be made by xdformat',
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
    return "Creates an xdformat database";
}

sub help_detail {
    return help_brief();
    return <<EOS
EOS
}

#< Pre and Post Execute Methods >#
sub _pre_execute_methods {
    return (qw/ _handle_db _verify_fastas /);
}

sub _handle_db {
    my $self = shift;

    for my $file ( $self->xdformat_files ) {
        if ( $self->overwrite_db ) {
            unlink $file
        }
        else {
            $self->error_message(
                sprintf(
                    'Files for database (%s) already exist.  Remove, or indicate to overwrite the db.',
                    $self->database,
                )
            );
            return;
        }
    }

    return 1;
}

sub _verify_fastas {
    my $self = shift;

    unless ( $self->fasta_files ) {
        $self->error_message("No fasta files to create database from");
        return;
    }

    my @missing_fastas;
    for my $fasta_file ( $self->fasta_files ) {
        push @missing_fastas, $fasta_file unless -e $fasta_file;
    }

    if ( @missing_fastas ) {
        $self->error_message(
            sprintf(
                'FASTA files (%s) do not exist',
                join(', ', @missing_fastas),
            )
        );
        return;
    }

    return 1;
}

sub _post_execute_methods {
    return (qw/ _verify_db /);
}

sub _verify_db {
    my $self = shift;

    return Genome::Model::Tools::WuBlast::Xdformat::Verify->execute(
        database => $self->database,
        db_type => $self->db_type,
    );
}

#< Operation >#
sub _operation_name {
    return 'create';
}

sub _operation_character {
    return 'o';
}

#< Db files >#
sub xdformat_files {
    my $self = shift;
    
    #print Dumper([ glob( sprintf('%s.x%s[a-z]', $self->database, $self->db_type) )]);
    return glob( sprintf('%s.x%s[a-z]', $self->database, $self->db_type) );
}

1;

#$HeadURL$
#$Id$
