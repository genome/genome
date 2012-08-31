package Finishing::Assembly::Commands::AnalyzeJoin;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::Commands::Base';

use Data::Dumper;
use Finishing::Assembly::ContigTools;
use IO::File;

my %organism_name :name(organism_name:r)
    :desc('Organism name');
my %assembly_name :name(assembly_name:r)
    :desc('Assembly name');
my %phd; #TODO
my %file :name(stats_file:o)
    :isa(file_rw)
    :desc('Stats file to print join info to');
my %lc_name :name(left_contig_name:r)
    :desc('Left contig name');
my %rc_name :name(right_contig_name:r)
    :desc('Right contig name');

my %phd_obj; #TODO
my %lc :name(_left_contig:p);
my %rc :name(_right_contig:p);
my %ctp :name(_contig_tools_params:p)
    :ds(hashref);

sub START
{
    my $self = shift;

    my $assembly = $self->_factory->get_assembly
    (
        name => $self->Assembly_name,
        organism_name => $self->orgasnism_name,
    );

    $self->fatal_msg
    (
        sprintf
        (
            'Can\'t get assembly (assembly name: %s, organism name: %s)', 
            $self->assembly_name, 
            $self->organism_name,
        )
    ) unless $assembly;

    my $left_contig = $assembly->get_contig( $self->left_contig_name );
    $self->fatal_msg
    (
        sprintf('Can\'t get left contig (%s)', $self->left_contig_name)
    ) unless $left_contig;

    my $right_contig = $assembly->get_contig( $self->right_contig_name );
    $self->fatal_msg
    (
        sprintf('Can\'t get right contig (%s)', $self->right_contig_name)
    ) unless $right_contig;

    my $contig_tools = Finishing::Assembly::ContigTools->new;
    my $po;# = GSC::IO::Assembly::PhdDB->new;

    my %params = 
    (
        cutoffs =>
        {
            hqlength => 1000, # $self->hq_length,
            hq_percent_identity => 90, # $self->hq_pecent_id,
        },
    );

    if ( my $file = $self->stats_file )
    {
        my $fh = IO::File->new(">> $file");
        $self->fatal_msg("Can't open file ($file): $!") unless $fh;
        $params{statsfh} = $fh;
    }
    
    $self->_contig_tools_params(\%params);
    
    return 1;
}

sub execute
{
    my $self = shift;

    my $contig_tools = Finishing::Assembly::ContigTools->new;
    $contig_tools->merge
    (
        $self->_left_contig,
        $self->_right_contig,
        $self->_phd, #TODO
        %{ $self->_contig_tools_params },
    );
    
    return 1;
}

1;

=pod

=head1 Name

Finishing::Assembly::Commands::JoinContigs

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

#$HeadURL$
#$Id$
