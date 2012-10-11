package Genome::DruggableGene::Command::DumpDatabase;

use strict;
use warnings;
use Genome;

class Genome::DruggableGene::Command::DumpDatabase {
    is => 'Genome::Command::Base',
    has => [
        output_file => {
            is => 'Path',
            doc => 'Path for the database dump',
        },
        postgres_host => {
            is => 'text',
            doc => 'The hostname of the postgres box from which data will be dumped',
            default => 'postgres',
        },
    ],
};

sub help_brief { 'Dump the dgidb database using pg_dump' }

sub help_synopsis { help_brief() }

sub help_detail { help_brief() }

sub execute {
    my $self = shift;
    my $postgres_host = $self->postgres_host;
    my $output_file = $self->output_file;
    system("pg_dump -a -f $output_file -h $postgres_host -U genome -n dgidb");
    1;
}

1;
