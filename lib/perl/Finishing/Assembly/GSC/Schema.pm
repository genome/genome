package Finishing::Assembly::GSC::Schema;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finfo::Singleton';

use Data::Dumper;

my %auto_commit :name(_auto_commit:p);

sub connect
{
    my ($self, $dbi, $user, $pw, $db_params) = @_;

    $self = $self->instance unless ref($self);
    
    require GSCApp; 
    unless ( App::Init->initialized )
    {
        App::DB->db_access_level('rw');
        App->init;
    }

    $self->_auto_commit( $db_params->{AutoCommit} ) if exists $db_params->{AutoCommit};

    return $self;
}

sub commit
{
    my $self = shift;

    App::DB->sync_database
        or $self->fatal_msg("Can't sync");

    App::DB->commit
        or $self->fatal_msg("Can't commit");

    return 1;
}

sub DEMOLISH
{
    return 1;
    my $self = shift;

    $self->commit if $self->_auto_commit;

    return 1;
}

sub get_project
{
    my ($self, $name) = @_;

    Finfo::Validate->validate
    (
        attr => 'gsc project name',
        value => $name,
        isa => 'string',
        msg => 'fatal',
    );

    return GSC::Project->get(name => $name);
}

sub create_project
{
    my ($self, %p) = @_;

    my $name = delete $p{name};
    my $oraganism_name = delete $p{organism};

    my $project = $self->get_project($name);
    return $project if $project;

    my $directory;
    if ( my $base_dir = delete $p{base_directory} )
    {
        $directory = sprintf('%s/%s', $base_dir, $name);
    }
    else
    {
        $directory = Finishing::Assembly::Project::Utils->instance->determine_and_create_projects_directory
        (
            $name,
            $oraganism_name || 'unknown',
        );
    }

    my $ps = GSC::ProcessStep->get
    (
        process_to => 'new project',
        process => 'new project',
    );

    $self->fatal_msg("Can't get new project process step") unless $ps;

    my $pse = $ps->execute
    (
        name => $name,
        consensus_directory => $directory,
        project_status => 'prefinish_done',
        target => 0,
        purpose => 'finishing',
        group_name => 'crick',
        priority => 0,
    );

    $self->fatal_msg
    (
        "Can't execute new project process step for $name"
    ) unless $pse;

    $project = GSC::Project->get(name => $name);

    $self->fatal_msg("Executed pse, but can't project from db") unless $project;

    $self->commit;
    
    return $project;
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

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/GSCSchema.pm $
#$Id: GSCSchema.pm 31312 2007-12-26 23:11:14Z ebelter $

