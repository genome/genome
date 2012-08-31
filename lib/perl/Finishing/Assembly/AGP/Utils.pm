package Finishing::Assembly::AGP::Utils;

use strict;
use warnings;

use base 'Finfo::Singleton';

use Data::Dumper;
use Finishing::Assembly::AGP::Info;

sub agp_info : PRIVATE
{
    return Finishing::Assembly::AGP::Info->instance;
}

sub string_to_agp_hashref
{
    my ($self, $string) = @_;

    $self->_enforce_instance;
    
    return unless Finfo::Validate->validate
    (
        attr => 'string to convert to agp',
        value => $string,
        msg => 'fatal',
    );
    
    my (@values) = split(/\s+/, $string);

    $self->error_msg("Can't split string\n$string")
        and return unless @values;

    my $comp_pos = $self->agp_info->component_position;

    my @attrs = $self->agp_info->attributes_for_component_type( $values[$comp_pos] );

    return unless @attrs;

    $self->error_msg
    (
        sprintf
        (
            "There should be %d columns, but found only %d on line:\n%s", 
            scalar @attrs, 
            scalar @values,
            $string,
        )
    )
        and return unless @attrs == @values;
    
    my %agp;
    @agp{ @attrs } = @values;

    return \%agp;
}

sub agp_to_string
{
    my ($self, $agp) = @_;

    $self->_enforce_instance;

    ( ref $agp eq 'HASH' )
    ? return $self->_hashref_to_string($agp)
    : return $self->_obj_to_string($agp);
}

sub _obj_to_string : PRIVATE
{
    my ($self, $agp) = @_;

    Finfo::Validate->validate
    (
        attr => 'agp object',
        value => $agp,
        isa => 'object',
        msg => 'fatal',
    );

    return join
    (
        "\t", 
        map { $agp->$_ } $self->agp_info->attributes_for_component_type( $agp->component_type ), 
    );
}

sub _hashref_to_string : PRIVATE
{
    my ($self, $agp) = @_;

    Finfo::Validate->validate
    (
        attr => 'agp hashref',
        value => $agp,
        ds => 'hashref',
        msg => 'fatal',
    );

    return join
    (
        "\t", 
        map { $agp->{$_} } $self->agp_info->attributes_for_component_type
        (
            $agp->{component_type} 
        ), 
    );
}

sub fasta_index_to_string
{
    my ($self, $fasta_index) = @_;

    $self->_enforce_instance;
    
    Finfo::Validate->validate
    (
        attr => 'fasta index',
        value => $fasta_index,
        isa => 'object Finishing::Assembly::AGP::FastaIndex',
        msg => 'fatal',
    );

    return join
    (
        "\t", 
        map { $fasta_index->$_ } $self->agp_info->attributes_for_component_type
        (
            $fasta_index->component_type 
        ), 
    );
}

1;

=pod

=head1 Name

Finishing::Assembly::AGP::Utils

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

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Finishing/Assembly/AGP/Utils.pm $
#$Id: Utils.pm 30518 2007-11-30 22:45:36Z ebelter $
