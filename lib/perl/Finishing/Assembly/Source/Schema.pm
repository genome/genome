package Finishing::Assembly::Source::Schema;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finfo::Singleton';

use Data::Dumper;
use Finishing::Assembly::Cache;
#use Finishing::Assembly::Project::Utils;

# models
use Finishing::Assembly::Source::AssembledRead;
use Finishing::Assembly::Source::Assembly;
use Finishing::Assembly::Source::Contig;
use Finishing::Assembly::Source::Organism;
use Finishing::Assembly::Source::Project;
use Finishing::Assembly::Source::Tags;

my %cache :name(_cache:p) :isa(object) :default( Finishing::Assembly::Cache->new() );
my %ctg_cache :name(_contig_cache:p) :isa(object) :default( Finishing::Assembly::Cache->new() );
my %proj_cache :name(_project_cache:p) :isa(object) :default( Finishing::Assembly::Cache->new() );

sub connect
{
    return __PACKAGE__->instance;
}

#- CONTIG -#
sub create_contig
{
    my ($self, %p) = @_;
 
    return Finishing::Assembly::Source::Contig->new(%p);
}

sub copy_contig
{
    my ($self, $obj) = @_;
    my @accessors = Finishing::Assembly::Source::Contig->accessors;
    my %p;
    foreach ( @accessors){
        $p{$_} = $obj->$_;
    }
    return Finishing::Assembly::Source::Contig->new(%p);
}

#- ASSEMBLED READ -#
sub create_assembled_read
{
	my ($self, %p) = @_;

	return Finishing::Assembly::Source::AssembledRead->new(%p);
}

sub copy_assembled_read
{
    my ($self, $obj) = @_;

    my @accessors = Finishing::Assembly::Source::AssembledRead->accessors;
    my %p;
    foreach ( @accessors )
    {
        $p{$_} = $obj->$_;
    }
    
    return Finishing::Assembly::Source::AssembledRead->new(%p);
}

#- CONSENSUS TAG -#
sub create_consensus_tag {
    my ($self, %p) = @_;
    return Finishing::Assembly::Source::ConsensusTag->new(%p);
}

sub copy_consensus_tag
{
    my ($self, $obj) = @_;
    my @accessors = Finishing::Assembly::Source::ConsensusTag->accessors;
    my %p;
    foreach (@accessors){
        $p{$_} = $obj->$_;
    }
    return Finishing::Assembly::Source::ConsensusTag->new(%p);
}

sub create_read_tag { 
    my ($self, %p) = @_;
    return Finishing::Assembly::Source::ReadTag->new(%p);
}

sub copy_read_tag
{
    my ($self, $obj) = @_;
    my @accessors = Finishing::Assembly::Source::ReadTag->accessors;
    my %p;
    foreach (@accessors){
        $p{$_} = $obj->$_;
    }
    return Finishing::Assembly::Source::ReadTag->new(%p);
}

sub create_assembly_tag { 
    my ($self, %p) = @_;
    return Finishing::Assembly::Source::AssemblyTag->new(%p);
}

sub copy_assembly_tag
{
    my ($self, $obj) = @_;
    my @accessors = Finishing::Assembly::Source::AssemblyTag->accessors;
    my %p;
    foreach (@accessors){
        $p{$_} = $obj->$_;
    }
    return Finishing::Assembly::Source::AssemblyTag->new(%p);
}

#- PROJECT -#
sub project_utils
{
    return Finishing::Assembly::Project::Utils->instance;
}

sub get_or_create_project
{
    my $self=shift;

    return $self->create_project(@_);
}

sub get_project
{
    return;
}

sub create_project
{
    my ($self, %p) = @_;

    my $name = delete $p{name};
    Finfo::Validate->validate
    (
        attr => 'project name',
        value => $name,
        isa => 'string',
        msg => 'fatal',
    );

    my $base_directory = delete $p{base_directory};
    unless ( $base_directory )
    {
        $base_directory = $self->project_utils->source_projects_base_directory;
    }

    my $directory = sprintf('%s/%s', $base_directory, $name);

    unless ( -d $directory )
    {
        mkdir $directory
            or $self->fatal_msg("Can't create directory ($directory) for $name\: $!");
    }

    return Finishing::Assembly::Source::Project->new
    (
        name => $name,
        directory => $directory,
    );
}

1;

=pod

=head1 Name

Finishing::Assembly::GSCSchema

=head1 Synopsis

=head1 Usage

=head1 Methods

=head1 See Also

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/TmpSchema.pm $
#$Id: TmpSchema.pm 31127 2007-12-18 20:39:53Z ebelter $

