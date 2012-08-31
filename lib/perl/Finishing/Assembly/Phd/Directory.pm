package Finishing::Assembly::Phd::Directory;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use Finishing::Assembly::Phd::Exporter;
use Finishing::Assembly::Phd::FileReader;
use IO::Dir;
use IO::File;

require Cwd;
require Date::Format;
require File::Copy;

my %directory :name(directory:r) :isa(dir_rw);

sub connect
{
    my ($class, $directory) = @_;

    return $class->new
    (
        #directory => Cwd::abs_path($directory), # causes problems
        directory => $directory,
    );
}

sub disconnect
{
    return 1;
}

sub reader
{
    return Finishing::Assembly::Phd::FileReader->instance;
}

#- GENERAL -#
sub base_names {
    my $self = shift;

    my %names;
    for my $phd ( $self->phd_names ) {
        $phd =~ s/\.phd\.\d+$//;
        $names{$phd} = 1;
    }

    return keys %names;
}

sub phd_names {
    my $self = shift;

    my $directory = $self->directory;
    my @phd_files = glob("$directory/*");

    return grep { s#$directory/##; m#\.phd\.# } @phd_files;
}

sub phd_file {
    my ($self, $name) = @_;

    return sprintf('%s/%s', $self->directory, $name);
}

sub phd_file_name {
    return phd_file(@_);
}

sub phd_exists
{
	my ($self, $name) = @_;

    my $file = $self->phd_file_name($name);

	return unless -s $file;

    return $file;
}

sub get_phd
{
	my ($self, $phd_name) = @_;
	
    my $file = $self->phd_exists($phd_name)
        or return;

    my $fh = IO::File->new("< $file")
        or $self->fatal_msg( sprintf( 'Can\'t open file (%s): %s', $file, $!) );
    my $read = $self->reader->execute($fh);
    $fh->close;

    $read->phd_file($phd_name);
    
    return $read;
}

#- LATEST PHD -#
sub latest_phd_iteration
{
    my ($self, $name) = @_;

    my $latest_iteration = 0;
    foreach my $phd_name ( glob( sprintf('%s/%s.phd.*', $self->directory, $name) ) )
    {
        my ($iteration) = $phd_name =~ /\.phd\.(\d+)$/;
        $latest_iteration = $iteration if $iteration > $latest_iteration;
    }

    return $latest_iteration;
}

sub latest_phd_name
{
    my ($self, $name) = @_;

    my $latest_iteration = $self->latest_phd_iteration($name);
    
    return if $latest_iteration == 0;
    
	return sprintf('%s.phd.%s', $name, $latest_iteration);
}

sub latest_phd_file
{
    my ($self, $name) = @_;

    my $phd_name = $self->latest_phd_name($name);

    return unless $phd_name;

    return $self->phd_file_name($phd_name);
}

sub latest_phd
{
    my ($self, $name) = @_;

    my $phd_name = $self->latest_phd_name($name);

    return unless $phd_name;
    
    return $self->get_phd($phd_name);
}

sub get_latest_phd
{
    my ($self, $name) = @_;

    my $phd_name = $self->latest_phd_name($name);

    return unless $phd_name;
    
    return $self->get_phd($phd_name);
}

#- NEXT PHD -#
sub next_phd_name
{
    my ($self, $name) = @_;

    return sprintf('%s.phd.%d', $name, $self->latest_phd_iteration($name) + 1); 
}

sub next_phd_file
{
    my ($self, $name) = @_;

    return $self->phd_file_name( $self->next_phd_name($name) );
}

#- CONVENIENCE METHODS -#
sub add_tag_to_phd
{
    my ($self, $phd_name, $tags) = @_;

    my $read = $self->latest_phd( $phd_name );

    $self->error_msg("Can't find phd for name ($phd_name)") unless $read;

    $self->_set_defaults_for_tags($read, $tags);
    push @{ $read->tags }, @$tags;

    my $tmp_file = sprintf('%s/tmp', $self->directory);
    unlink $tmp_file if -e $tmp_file;

    my $exporter = Finishing::Assembly::Phd::Exporter->new
    (
        'read' => $read,
        file => $tmp_file,
    );
    $exporter->execute;

    unlink $read->phd_file;
    File::Copy::copy($tmp_file, $read->phd_file);
    unlink $tmp_file;

    return 1;
}

sub _set_defaults_for_tags : PRIVATE
{
    my ($self, $phd, $tags) = @_;

    foreach my $tag ( @$tags )
    {
        $tag->type('comment') unless $tag->type;
        $tag->source('phd_schema') unless $tag->source;
        $tag->date( Date::Format::time2string('%y/%d/%m %H:%M:%S', time) ) unless $tag->date;
        $tag->start(1) unless $tag->start;
        $tag->stop( $phd->length ) unless $tag->stop;
    }

    return 1;
}

1;		

=pod

=head1 Name

 Finishing::Assembly::Phd
 
  > Object oriented phd/phd.ball file reader/writer

=head1 Synopsis

=head1 Usage

 my $phd_schema = Finishing::Assembly::Phd::Directory->new
 (
    directory => './phd_dir',
 );

 my @phd_names = $phd_schema->phd_names;

 my $phd = $phd_schema->get_phd("vef07");
    
=head1 Methods

=head1 Disclaimer

=head1 Author(s)

=cut

#$HeadURL$
#$Id$
