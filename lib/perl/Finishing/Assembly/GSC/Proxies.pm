package Finishing::Assembly::GSC::Proxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::Proxy';

sub _get_only : RESTRICTED
{
    my ($self, $method) = @_;

    return $self->source->$method;
}

##################################################################################

package Finishing::Assembly::GSC::ProjectProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::GSC::Proxy';

my %dir :name(_dir:p);

sub methods_for_source_method
{
    return GSC::Project->property_names;
}

sub methods : CUMULATIVE
{
    return (qw/ get_only /);
}

sub _methods_for_get_only
{
    return 
    (qw/
        organism_name organism_nickname 
        finisher_unix_login finisher_group claim_date
        prefinisher_unix_login
        get_accession_numbers
        get_neighbor_info_from_tilepath
        get_species_name get_chromosome
        get_projects_submission_pses get_projects_last_submission_pse
        /);
}

sub owner
{
    return shift->source->finisher_unix_login;
}

sub _set_dir : PRIVATE
{
    my ($self, $dir) = @_;

    $self->source->consensus_directory($dir);

    return $self->_dir($dir);
}

sub directory
{
    my ($self, $new_dir) = @_;
    
    if ( $new_dir )
    {
        return $self->set_dir($new_dir);
    }

    if ( $self->_dir )
    {
        return $self->_dir;
    }
    
    my $con_dir = $self->source->consensus_directory;

    return $con_dir if $con_dir; 
    
    my $abs_dir = $self->source->consensus_abs_path;
    
    return $self->_set_dir($abs_dir);
}

sub create_wg_clone_link
{
    my $self = shift;

    return 1 if GSC::CloneProject->get(project_id => $self->source->project_id);

    my ($clone) = GSC::Clone->get
    (
        sql =>
        "select * from clones where ct_clone_type = 'genome' and cs_clone_status = 'active' and clone_name like 'C\\_%' escape '\\'"
    );

    $self->error_msg("Can't get wg clone")
        and return unless $clone;

    my $new_cp = GSC::CloneProject->create
    (
        project_id => $self->source->project_id,
        clo_id => $clone->clo_id
    );

    $self->fatal_msg("Can't create clone proj link for " . $self->source->name)
        and return unless $new_cp;

    return $new_cp;
}

##################################################################################

package Finishing::Assembly::GSC::ScaffoldProxy;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::GSC::Proxy';

sub methods_for_source_method
{
    return
    {
        get_contig_iterator => 'contig_iterator',
    };
}

sub name
{
    my $self = shift;

    return $self->source->sequence_item_name;
}

sub padded_length
{
    my $self = shift;

    $self->source->seq_length(@_);
}

sub unpadded_length
{
    my $self = shift;

    return $self->source->get_unpadded_position( $self->source->seq_length );
}

##################################################################################

1;

=pod

=head1 Name

Finishing::Project::Proxy

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

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/GSProxies.pm $
#$Id: GSCProxies.pm 31311 2007-12-26 23:10:43Z ebelter $

