package Finishing::Assembly::Phd::FastaAndQualDB;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;
use File::Basename;

use Bio::SeqIO;
use Data::Dumper;
use Finishing::Assembly::Phd::Exporter;
use Finishing::Assembly::Factory;
#use Finishing::Assembly::Source::AssembledRead;
use Finishing::Assembly::Source::Tags;

my %file :name(file:r) :isa(file_rw);
my %dbh :name(_dbh:p);
my %dir :name(_dir:p);

sub connect
{
    my ($class, $file) = @_;

    return $class->new
    (
        file => $file,
    );
}

sub START
{
    my $self = shift;
    my $file = $self->file;
    
    $self->fatal_msg("sqlite DB file: $file not exist") unless -s $file;
    
    my $dbh = DBI->connect
    (
        'dbi:SQLite:dbname=' . $file, '', '', { AutoCommit => 0, RaiseError => 1},
    )
        or $self->fatal_msg("Can't connect: " . $DBI::errstr);

    #my $sth = $dbh->prepare("select value from general_info where attribute = 'directory'");
    #$sth->execute or $self->fatal_msg($DBI::errstr);
    #my @ary = $sth->fetchrow_array;
    #$self->fatal_msg("No directory in database") unless @ary;
    #set fasta and qual files in the same dir as sqlite DB file.
    
    $self->_dir(dirname($file));
    $self->_dbh($dbh);

    return 1;
}

sub disconnect
{
    my ($self) = @_;

    $self->_dbh->disconnect if $self->_dbh;

    return 1;
}

sub DEMOLISH
{
    my $self = shift;

    return $self->disconnect;
}

sub _execute_st : PRIVATE
{
    my ($self, $st) = @_;

    my $sth = $self->_dbh->prepare($st)
        or $self->fatal_msg("Can't prepare statement:\n$st\nERROR:\n$DBI::errstr");
    $sth->execute
        or $self->fatal_msg("Can't prepare statement:\n$st\nERROR:\n$DBI::errstr");

    return $sth;
}

#- GENERAL PHD -#
sub phd_names
{
    my $self = shift;

    my $sth = $self->_execute_st('select name from read_index');

    my $aryref = $sth->fetchall_arrayref;

    return unless $aryref and @$aryref;
    
    return map { @$_ } @$aryref;
}

sub phd_exists
{
    my ($self, $name) = @_;

    return $self->_execute_st( sprintf('select * from read_index where name = \'%s\'', $name) )->fetchrow_array;
}

sub get_phd
{
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

    my $sth = $self->_dbh->prepare(
        qq(
            select f1.name, r.fasta_position, f2.name, r.qual_position
            from read_index r
            join file_id f1 on f1.id = r.fasta_file_id
            join file_id f2 on f2.id = r.qual_file_id
            where r.name = '$name'
        )
    );
 
    $sth->execute;
    my @ary = $sth->fetchrow_array;

    return unless @ary;

    my $fasta_file = sprintf('%s/%s', $self->_dir, $ary[0]);
    my $fasta_fh = IO::File->new("< $fasta_file");
    $fasta_fh->seek($ary[1], 0);
    my $fasta_bio = Bio::SeqIO->new(-fh => $fasta_fh, -format => 'Fasta');
    my $fasta = $fasta_bio->next_seq;

    my $qual_file = sprintf('%s/%s', $self->_dir, $ary[2]);
    my $qual_fh = IO::File->new("< $qual_file");
    $qual_fh->seek($ary[3], 0);
    my $qual_bio = Bio::SeqIO->new(-fh => $qual_fh, -format => 'qual');
    my $qual = $qual_bio->next_seq;

    main->fatal_msg
    (
        sprintf("Name in db does not match fasta and qual", $fasta->id, $qual->id)
    ) unless $name eq $fasta->id and $fasta->id eq $qual->id;

    my %phd = 
    (
        name => $name,
        phd_name => $phd_name,
        base_string => $fasta->seq,
        qualities => $qual->qual,
    );

    my $ds = ' ' . $fasta->desc;
    $ds =~ s/ (\w+): /|$1|/g;
    my @tokens = split(/\|/, $ds);
    shift @tokens;
    my %desc = @tokens;
    my %comments;
    foreach my $key ( keys %desc )
    {
        $comments{ lc($key) } = $desc{$key};
    }
    $phd{comments} = \%comments;

    undef $fasta_bio;
    $fasta_fh->close;
    undef $qual_bio;
    $qual_fh->close;

    #print Dumper(\%phd);
    #return Finishing::Assembly::Source::AssembledRead->new(%phd);
    my $fac = Finishing::Assembly::Factory->connect('source');
    return $fac->create_assembled_read(%phd);
}

################
#- LATEST PHD -#
################

sub latest_phd_name
{
    my ($self, $name) = @_;

    #TODO 
    my $latest_iteration = 0;
    foreach my $phd_name ( glob( sprintf('%s/%s.phd.*', $self->_dir, $name) ) )
    {
        my ($iteration) = $phd_name =~ /\.phd\.(\d+)$/;
        $latest_iteration = $iteration if $iteration > $latest_iteration;
    }

	return sprintf('%s.phd.%s', $name, $latest_iteration);
}

sub latest_phd_file
{
    my ($self, $name) = @_;

    my $phd_name = $self->latest_phd_name($name);

    return unless $phd_name;

    return $self->phd_file($phd_name);
}

sub get_latest_phd
{
    my ($self, $name) = @_;

    my $phd_name = $self->latest_phd_name($name);

    return unless $phd_name;
    
    return $self->get_phd($phd_name);
}

1;		

=pod

=head1 Name

 Finishing::Assembly::Phd
 
  > Object oriented phd/phd.ball file reader/writer

=head1 Synopsis

 my $phd_object = Finishing::Assembly::Phd->new
 (
    input_directory => "inputdirname",
 );

 my @phd_names = $phd_object->get_phd_names();
 my $phd = $phd_object->get_phd("vef07");

    
=head1 Description

Finishing::Assembly::Phd takes either a Phd file, and allows the user to get Contig objects from the ace file, edit them, and write the file back to the hard disk when finished.

=head1 METHODS

=head2 get_phd_names

=head2 get_phd

=head2 get_latest_phd

=head2 get_latest_phd_name

=head2 get_next_phd_name

=head2 write_phd

=head2 write_phd_ball

=head2 add_phd

=head1 Disclaimer

=head1 Author(s)

=cut

#$HeadURL$
#$Id$
