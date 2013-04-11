package Genome::Db::Command::List;
use strict;
use warnings;

class Genome::Db::Command::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name => { is_constant => 1, value => 'Genome::Db' }, 
    ],
    doc => 'list databases',
};

sub help_synopsis {
return <<EOS
    GENOME_DB=\$HOME/mydb:/gsc/scripts/opt/genome/db genome db list --filter "source_dir like '%mydb%'"
EOS
}

1;

