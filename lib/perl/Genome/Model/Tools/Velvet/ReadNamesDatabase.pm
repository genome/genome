package Genome::Model::Tools::Velvet::ReadNamesDatabase;

use strict;
use warnings;

use DBI;
use Genome;

use Bio::SeqIO;

class Genome::Model::Tools::Velvet::ReadNamesDatabase {
    is => 'Genome::Model::Tools::Velvet::Base',
    has => [
	sequences_file => {
	    is => 'Text',
	    doc => 'Velvet generated Sequences file',
	    is_optional => 1,
	},
	directory => {
	    is => 'Text',
	    doc => 'Velvet assembly directory',
#	    is_optional => 1,
	},
	create_new => {
	    is => 'Boolean',
	    doc => 'Overwrite existing seqlite db file',
	    is_optional => 1,
	},
    ],
};

sub help_brief {
    'Tool to create seqlite database of input read names indexed by read number in afg file',
}

sub help_synopsis {
    return <<EOS
gmt velvet read-names-database --seq-file /gscmnt/111/velvet_assembly/Sequences --directory /gscmnt/111/velvet_assembly
EOS
}

sub help_detail {
    return <<EOS
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    unless ($self->sequences_file and -s $self->sequences_file) {
	my $sequences_file = $self->directory.'/Sequences';
	unless (-s $sequences_file) {
	    $self->error_message("Can't find velvet sequences file: ".$self->directory.'/Sequences'.
				 "\n\tYou must supply it");
	    return;
	}
	$self->sequences_file($sequences_file);
    }

    return $self;
}

sub execute {
    my $self = shift;

    #force make a new one
    if ($self->create_new) {
	#unlink $self->directory.'/velvet_reads.sqlite' if -s $self->directory.'/velvet_reads.sqlite';
	unlink $self->read_names_sqlite;
	unless ($self->create_db) {
	    $self->error_message("Failed to create db");
	    return;
	}
	return 1;
    }
    #if it doesn't exist make one
    #unless (-s $self->directory.'/velvet_reads.sqlite') {
    unless (-s $self->read_names_sqlite) {
	unless ($self->create_db) {
	    $self->error_message("Failed to create db");
	    return;
	}
	return 1;
    }
    #use existing one
    $self->status_message("velvet_reads.sqlite exists\n\t".$self->read_names_sqlite);

    return 1;
}

sub create_db {
    my $self = shift;

    my $dbh = $self->_get_dbh();

    my $sth = $dbh->prepare('create table read_info (id integer, name string, position integer)');
    $sth->execute or return $self->error_message('Failed to create table read_info : '.$DBI::errstr);
    $dbh->commit;

    my $sth_insert_read = $dbh->prepare("insert into read_info (id, name, position) values (?,?,?)");
    die $dbh->errstr unless $sth_insert_read; #TODO ??

    my $sth_get_read_name_and_position = $dbh->prepare( qq (select name, position from read_info where id = ?) );
    die $dbh->errstr unless $sth_get_read_name_and_position; #??

    #load from sequences file
    my $seq_fh  = Genome::Sys->open_file_for_reading($self->sequences_file) or return;
    my $seek_pos = $seq_fh->tell;
    my $io = Bio::SeqIO->new(-format => 'fasta', -fh => $seq_fh);

    while (my $seq = $io->next_seq) {
	my $read_name = $seq->primary_id;
	my ($read_index) = $seq->desc =~ /^(\d+)\s+\d+$/;
	unless ($read_index) {
	    $self->error_message("Failed to get read index number from seq->desc: ".$seq->desc);
	    return;
	}
	$sth_insert_read->execute($read_index, $read_name, $seek_pos) or
	    return $self->error_message("Failed to insert for $read_name : ".$DBI::errstr);

	$seek_pos = $seq_fh->tell;
    }
    $seq_fh->close;

    my $ssth = $dbh->prepare('create unique index ids on read_info(id)');
    $ssth->execute or return $self->error_message('Failed to create index ids : '.$DBI::errstr);
    $dbh->commit;

    return 1;
}

sub get_read_name_from_afg_index {
    my ($self, $read_id) = @_;

    unless (-s $self->read_names_sqlite) {
	$self->error_message("Can't find reads sqlite database: ".$self->read_names_sqlite);
	return;
    }

    my $dbh = $self->_get_dbh();

    my $sth = $dbh->prepare( qq(select name, position from read_info where id = '$read_id'));
    $sth->execute or return $self->error_message("Failed to select for read: $read_id ".$DBI::errstr);

    my @out = $sth->fetchrow_array;
	return $self->error_message("Got nothing from sqlite select for read: $read_id") unless @out;

    my ($read_name, $pos) = @out;

    return $read_name, $pos;
}

sub _get_dbh {
    my $self = shift;

    my $db_file = $self->read_names_sqlite; #name it with a file ending??

    my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", '', '', { AutoCommit => 0, RaiseError => 1 })
        or return $self->error_message("Failed to connect to db ($db_file): " . $DBI::errstr);

    return $dbh;
}

1;
