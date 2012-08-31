package Genome::Model::Tools::Pcap::Config;

our $VERSION = 0.01;

use strict;
use warnings;
use Carp;
use base qw(Class::Accessor::Fast);

use IO::File;
use Storable;

#Genome::Model::Tools::Pcap::Config->mk_accessors(qw (

sub new 
{
	croak("__PKG__:new:no class given, quitting") if @_ < 1;
    my ($caller, %params) = @_;
    my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = {};
    bless ($self, $class);
	if (!exists $params{config_file})
	{
		$self->{config}{config_file} = "$ENV{HOME}/.acerc";
	}
	# my $fh = IO::File->new($self->config_file);
    if(! -e $self->{config}{config_file})
	{
		$self->_init;
	}
	else
	{	
		$self->load;
	}
	return $self;
}

sub _init
{
	my ($self) = @_;
	my $config = {};
	$config->{db_type} = 'SQLite';
#this is the dbi dsn for all db indexes
#and is the filesystem path for storable indexes
	$config->{user_name} = 'finfo_admin';
	$config->{password} = 'finfo_admin';
	$config->{mysql_def_dsn} = 'dbi:mysql:database=finfo;host=mysql1';
	$config->{sqlite_def_dsn} = 'dbi:SQLite:/gscmnt/936/info/jschindl/sqlite/assembly.db';
	$config->{default_assembly} = "assembly." . Genome::Sys->username;
	$self->{config} = $config;

}

sub save
{
	my ($self) = @_;
	Storable::nstore($self->{config},"$ENV{HOME}/.acerc");
}

sub load
{
	my ($self) = @_;
	$self->{config} = Storable::retrieve($self->{config}{config_file});
}

1;
