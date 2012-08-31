package Finishing::Assembly::Phd::Ball;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;
use File::Temp;

use Finishing::Assembly::Phd::FileReader;
use Finishing::Assembly::Phd::CreateSqliteFromBall;
use IO::File;

my %ball :name(ball:r) :isa(file_rw);
my %db_file :name(db_file:o) :isa(file_rw);
my %reader :name(reader:p) :isa(object);
my %dbh :name(dbh:p);


sub connect {
    my ($class, %param) = @_;
    return $class->new(%param);
}


sub START {
    my $self = shift;
    my $db_file = $self->db_file if $self->db_file;
    
    unless ($db_file) {
        $db_file = File::Temp->new(
            UNLINK   => 1, 
            DIR      => '/tmp', 
            SUFFIX   => '.sqlite',
            TEMPLATE => 'phd_ball_XXXXXX',
        )->filename;

        Finishing::Assembly::Phd::CreateSqliteFromBall->new(
            ball    => $self->ball,
            db_file => $db_file,
        )->execute;
    }

    $self->fatal_msg("sqlite dbfile $db_file not exist") unless -s $db_file;
    
    my $dbh = DBI->connect(
        'dbi:SQLite:dbname=' . $db_file, '', '', { AutoCommit => 0, RaiseError => 1},
    )
        or $self->fatal_msg("Can't connect: " . $DBI::errstr);

    $self->dbh($dbh);
    $self->reader(Finishing::Assembly::Phd::FileReader->instance);

    return 1;
}


sub get_phd {
	my ($self, $phd_name) = @_;
    
    my ($name, $iter);
    my $msg = 'This only works on .phd.1 version now';
    
    if ($phd_name =~ /\.phd\./) {
        ($name, $iter) = split('.phd.', $phd_name);
        unless ($iter == 1) {
            $self->warn_msg($msg);
            return;
        }
    }
    else {
        #$self->warn_msg($msg);
        $name = $phd_name;
    }
	
    my $fh = IO::File->new('<' . $self->ball)
        or $self->fatal_msg( sprintf( 'Can\'t open ball (%s): %s', $self->ball, $!) );

    my $sth = $self->dbh->prepare(
        qq(
            select position
            from phd_index
            where name = '$name'
        )
    );
    $sth->execute;
    
    my @out = $sth->fetchrow_array;
    return unless @out;
    
    $fh->seek($out[0], 0);
    my $phd = $self->reader->execute($fh);
    $fh->close;

    return $phd;
}


sub DEMOLISH {
    return shift->dbh->disconnect;
}


1;		

=pod

=head1 NAME

 Finishing::Assembly::Phd::Ball
 
  > Object oriented phd/phd.ball file reader

=head1 SYNOPSIS

 my $phd_object = Finishing::Assembly::Phd::Ball->connect(
    ball    => 'phd.ball',
    db_file => 'phd_ball.sqlite',
 );

 my $phd = $phd_object->get_phd('E3434mfdfdjfdfdbw');

    
=head1 DESCRIPTION

Finishing::Assembly::Phd:Ball takes phd.ball (required) file and sqlite db_file (optional). If without sqlite db_file, this module will create one in /tmp. sqlite store read name and their corresponding file index position. 

=head1 METHODS

=head2 get_phd

=head1 Disclaimer

=head1 Author(s)

=cut

#$HeadURL$
#$Id$
