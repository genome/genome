package Finishing::Assembly::Commands::ExportChromosomeAgp;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::Commands::Base';

use Data::Dumper;
use Finishing::Assembly::AGP::Writer;

my %assembly_name :name(assembly_name:r)
    :desc('Name of assembly to create ');
my %assem :name(_assembly:p);
my %chr_name :name(chromosome_name:r)
    :isa(string)
    :desc('Chromosome name to export');
my %chr :name(_chromosome:p)
    :isa('object');
my %_chr_name :name(_chromosome_name:p) # the name to print in the file
    :isa(string);
my %ordered :name(ordered:o)
    :isa(boolean)
    :default(0)
    :desc('Export the ordered contigs on the chromosome.  Default is to export the *unordered* contigs');
my %writer :name(writer:r)
    :isa('object Finishing::Assembly::AGP::Writer');
my %ob :name(_object_begin:p)
    :isa('int')
    :default(1);
my %pn :name(_part_num:p)
    :isa('int')
    :default(0);

sub START
{
    my $self = shift;

    my $organism = $self->_assembly->organism;
    my $chromosome = $organism->get_chromosome($self->chromosome_name);
    $self->fatal_msg
    (
        sprintf
        (
            'Chromosome (%s) not found for organism (%s)',
            $self->chromosome_name,
            $organism->name,
        )
    ) unless $chromosome;
    $self->_chromosome($chromosome);

    return 1;
}

sub execute
{
    my $self = shift;

    if ( $self->ordered )
    {
        $self->_chromosome_name(  'chr' . $self->_chromosome->name );
        return $self->_export_ordered_chromosome;
    }
    else
    {
        $self->_chromosome_name(  'chr' . $self->_chromosome->name . '_random' );
        return $self->_export_unordered_chromosome;
    }
}

#- UNORDERED -#
sub _export_unordered_chromosome
{
    my $self = shift;

    my $contig_iterator = $self->_assembly->contigs
    (
        chromosome_id => $self->_chromosome->id,
        ordered_on_chromosome => 0, # TODO add NULL
    );

    $self->_write_contig( $contig_iterator->next );

    while ( my $contig = $contig_iterator->next )
    {
        if ( my @gaps = $contig->left_gaps->all )
        {
            $self->_write_gaps(@gaps);
        }
        else
        {
            $self->_write_fake_gap;
        }

        $self->_write_contig($contig);
    }

    return 1;
}

sub _write_fake_gap
{
    my $self = shift;

    my $gap_length = 25;
    return $self->writer->write_one
    ( 
        {
            object => $self->_chromosome_name,
            object_begin => $self->_object_begin,
            object_end => $self->_add_to_object_begin( $gap_length ) - 1,
            part_number => $self->_increment_part_num,
            component_type => 'N',
            gap_length => $gap_length,
            gap_type => 'contig',
            linkage => 'no',
        }
    );
}

#- ORDERED -#
sub _export_ordered_chromosome
{
    my $self = shift;

    my $gap = $self->_write_first_left_gaps_and_contig;
    while ( my $contig = $gap->right_contig )
    {
        $self->_write_contig($contig);
        my @gaps = $contig->right_gaps->all;
        $self->_write_gaps(@gaps);
        $gap = $gaps[0];
    }

    return 1;
}

sub _write_first_left_gaps_and_contig : PRIVATE
{
    my $self = shift;

    my $chromosome = $self->_chromosome;
    my $first_scaffold = $chromosome->first_scaffold_for_assembly( $self->_assembly );
    $self->fatal_msg
    (
        "No first scaffold for chromosome " . $self->_chromosome_name
    ) unless $first_scaffold;

    my $contig = $first_scaffold->scaffold->first_contig;
    $self->_write_gaps( $contig->left_gaps->all );
    $self->_write_contig($contig);
    
    return $self->_write_gaps( $contig->right_gaps->all);
}

#- SHARED -#
sub _write_contig : PRIVATE
{
    my ($self, $contig) = @_;

    $self->fatal_msg("No contig to write") unless $contig;
    
    my $contig_length = $contig->length;
    #my $contig_length = $contig->length('unpadded');
    my $scaffold = $contig->scaffold;

    $self->writer->write_one
    ( 
        {
            object => $self->_chromosome_name,
            object_begin => $self->_object_begin,
            object_end => $self->_add_to_object_begin( $contig_length ) - 1,
            part_number => $self->_increment_part_num,
            component_type => 'W',
            component_id => sprintf
            (
                'Contig%d.%d', 
                $scaffold->scaffold_num, 
                $contig->contig_num
            ),
            component_begin => 1,
            component_end => $contig_length,
            orientation => $scaffold->orientation,
        }
    );

    return 1;
}

sub _write_gaps : PRIVATE
{
    my ($self, @gaps) = @_;

    return unless @gaps;
    
    foreach my $gap ( @gaps )
    {
        my $gap_length = $gap->length;
        $self->writer->write_one
        ( 
            {
                object => $self->_chromosome_name,
                object_begin => $self->_object_begin,
                object_end => $self->_add_to_object_begin( $gap_length ) - 1,
                part_number => $self->_increment_part_num,
                component_type => 'N',
                gap_length => $gap_length,
                gap_type => $gap->type,
                linkage => $gap->linkage,
            }
        );
    }

    return $gaps[0];
}

sub _increment_part_num : PRIVATE
{
    my $self = shift;

    $self->_part_num( $self->_part_num + 1 );
}

sub _add_to_object_begin : PRIVATE
{
    my ($self, $object_begin) = @_;

    $self->_object_begin( $self->_object_begin + $object_begin );
}

1;

=pod

=head1 Name

Finishing::Assembly::Commands::CreateOrganism

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

