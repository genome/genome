package Finishing::Assembly::Phd::CreateSqliteFromBall;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use IO::File;

my %ball :name(ball:r) :isa(file_rw) :desc('phd.ball file path');
my %db_file :name(db_file:o) :isa(string) :desc('desired output sqlite db file path') :default('phd_ball.sqlite');
my %dbh :name(dbh:p);


sub START {
    my $self = shift;

    my $dbfile = $self->db_file;
    unlink $dbfile if -e $dbfile;

    my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', { AutoCommit => 0, RaiseError => 1 })
        or $self->fatal_msg("Can't connect to db ($dbfile): " . $DBI::errstr);

    my $sth = $dbh->prepare('create table phd_index (name string, position integer)');
    $sth->execute or $self->fatal_msg($DBI::errstr);
    $dbh->commit;

    $sth = $dbh->prepare('create unique index names on phd_index(name)');
    $sth->execute or $self->fatal_msg($DBI::errstr);
    $dbh->commit;
    
    $self->dbh($dbh);
    return 1;
}
    

sub execute {
    my $self = shift;

    my $fh = IO::File->new('<' . $self->ball)
        or $self->fatal_msg( sprintf( 'Can not open file (%s): %s', $self->ball, $!) );

    my @values = ();
    
	while (my $line = $fh->getline) {
        next unless $line =~ /^BEGIN_SEQUENCE\s+(\S+)/;
        push @{$values[0]}, $1;
        push @{$values[1]}, $fh->tell - length $line;
        
        if (@{$values[0]} == 1000) {
            $self->_insert_values(\@values);
            @values = ();
        }
    }

    $self->_insert_values(\@values) if @values;
    return $self->db_file;
}


sub _insert_values {
    my ($self, $values) = @_;
    my $dbh = $self->dbh;

    my $sth = $dbh->prepare("insert into phd_index (name, position) values (?, ?)");
    
    $sth->execute_array({}, @$values) or $self->fatal_msg($DBI::errstr);
    $dbh->commit;

    return 1;
}


    
1;		

=pod

=head1 NAME

 Finishing::Assembly::Phd::CreateSqliteFromBall
 

=head1 SYNOPSIS

 my $create = Finishing::Assembly::Phd::CreateSqliteFromBall->connect(
    ball => $ball_file,
 );

 my $db_file = $create->execute;

    
=head1 DESCRIPTION

This module create simple sqlite DB from phd.ball. DB only contains two column name and file index position.
This can speed up getting phd from phd.ball.

=head1 METHODS

=head2 execute

=head1 Author(s)

=cut

#$HeadURL$
#$Id$
