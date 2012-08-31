package Finishing::Assembly::Commands::ImportChromosomeAgp;

use strict;
use warnings;
no warnings 'reserved';


use base('Finishing::Assembly::Commands::Base');

use Bio::SeqIO;
use Data::Dumper;
use Finishing::Assembly::AGP::Reader;

my %assem_name :name(assembly_name:r)
    :isa(string)
    :desc('Assembly name');
my %organism_name :name(organism_name:r)
    :isa(string)
    :desc('Organism name');
my %assem :name(_assembly:p)
    :isa('object');
my %org :name(_organism:p)
    :isa('object');
my %reader :name(reader:r)
    :isa('object Finishing::Assembly::AGP::Reader');
my %create_assembly :name(create_assembly:o)
    :isa(boolean)
    :default(0)
    :desc('Creates the scaffolds and contigs from the objects in the file.  Default is to update the existing assembly\'s scaffolds and contigs.');
my %ordered :name(ordered:o)
    :isa(boolean)
    :default(0)
    :desc('The AGP file consists of ordered contigs.  First scaffold and gaps between scaffolds will not be stored.');

my %on_first_scaffold :name(_on_first_scaffold:p)
    :isa(boolean)
    :default(1);

# on error, die is called several times, so only execute the fatal msg once
my $count = 0; 
local $SIG{__DIE__} = sub
{
    my ($msg) = @_;
    chomp $msg; 
    __PACKAGE__->fatal_msg($msg, { 'caller' => [ caller ] }) unless $count++; 
};

sub execute
{
    my $self = shift;
    
    $self->_organism( $self->_factory->get_organism($self->organism_name));
    $self->_assembly($self->_organism->get_assembly($self->assembly_name));

    my $method = ( $self->ordered )
    ? '_import_ordered_agps'
    : '_import_unordered_agps';

    return $self->$method;
}

sub _import_unordered_agps
{
    my $self = shift;

    $self->_factory->schema->txn_do
    (
        sub
        {
            my ($contig, $gap);
            while ( my $agp = $self->reader->next )
            {
                if ( exists $agp->{component_id} )
                {
                    $contig = $self->_factory->schema->txn_do
                    (
                        sub{ return $self->_create_or_update_contig_from_agp($agp); } 
                    );
                    my $var = 1;
                    $self->_update_contigs_left_gaps($contig, $gap) if $gap;
                }
                else
                {
                    if ( $agp->{gap_type} ne 'contig' )
                    {
                        ($gap) = $self->_factory->schema->txn_do
                        (
                            sub
                            {
                                return $self->_create_gaps_for_contig($contig, 'right', $agp);
                            } 
                        );
                    }
                    else
                    {
                        $gap = undef;
                    }
                }
            }
        }
    );

    return 1;
}

sub _import_ordered_agps
{
    my $self = shift;

    $self->_factory->schema->txn_do
    (
        sub
        {
            my ($contig, @left_gap_agps);
            while ( my $agp = $self->reader->next )
            {
                if ( exists $agp->{component_id} )
                {
                    $contig = $self->_factory->schema->txn_do
                    (
                        sub{ return $self->_create_or_update_contig_from_agp($agp); } 
                    );

                    $self->_factory->schema->txn_do
                    (
                        sub
                        {
                            return $self->_create_gaps_for_contig
                            (
                                $contig,
                                'left',
                                @left_gap_agps
                            ); 
                        } 
                    );

                    last;
                }
                push @left_gap_agps, $agp;
            }

            my @gaps;
            while ( my $agp = $self->reader->next )
            {
                if ( exists $agp->{component_id} )
                {
                    $contig = $self->_factory->schema->txn_do
                    (
                        sub{ return $self->_create_or_update_contig_from_agp($agp); } 
                    );
                    $self->_update_contigs_left_gaps($contig, @gaps) if @gaps;
                    @gaps = ();
                }
                else
                {
                    push @gaps, $self->_factory->schema->txn_do
                    (
                        sub
                        {
                            return $self->_create_gaps_for_contig($contig, 'right', $agp);
                        } 
                    );
                }
            }

            return 1;
        }
    );

    return 1;
}

sub _determine_method_for_object  #TODO, this no longer produces correct results.  The get_<obj> method is now a source method, with a non-hash argument, while the get_or_create is a related method, and still needs hash as args
{
    my ($self, $object) = @_;

    return sprintf
    (
        '%s_%s',
        ( $self->create_assembly ) ? 'get_or_create' : 'get',
        $object
    );
}

sub _create_gaps_for_contig : PRIVATE
{
    my ($self, $contig, $gap_side, @agps) = @_;

    my $create_gap_method = sprintf('create_%s_gap', $gap_side);

    my @gaps;
    foreach my $agp ( @agps )
    {
        push @gaps, $contig->$create_gap_method
        (
            'length' => $agp->{gap_length},
            type => $agp->{gap_type},
            linkage => $agp->{linkage}
        );
    }

    return @gaps;
}

sub _update_contigs_left_gaps : PRIVATE
{
    my ($self, $contig, @gaps) = @_;

    foreach my $gap ( @gaps )
    {
        $gap->set_right_contig($contig);
    }

    return 1;
}

sub _create_or_update_contig_from_agp : PRIVATE
{
    my ($self, $agp) = @_;

    my $organism = $self->_organism;
    my $assembly = $self->_assembly;

    my $comp_id =  $agp->{component_id};
    $comp_id =~ s/Contig//;
    my ($scaffold_num, $contig_num) = split(/\./, $comp_id);

    my $chromosome_name = $agp->{object};
    $chromosome_name =~ s/chr//i;
    $chromosome_name =~ s/_random//i;
    my $chromosome = $organism->get_or_create_chromosome(name => $chromosome_name);
    #TODO, need to fix the determine method to allow assembly creation
    #my $scaffold_method = $self->_determine_method_for_object('scaffold');
    #my $scaffold = $assembly->$scaffold_method(scaffold_num => $scaffold_num);
    my $scaffold = $assembly->get_scaffold($scaffold_num);
    #####################
    $scaffold->orientation( $agp->{orientation} ); 
    $scaffold->update;

    # only record the first scaffold
    if ( $self->_on_first_scaffold and $self->ordered )
    {
        $chromosome->set_first_scaffold($scaffold);
        $self->_on_first_scaffold(0);
    }

    #TODO needs to be fixed for create_assembly
    #my $contig_method = $self->_determine_method_for_object('contig');
    #my $contig = $scaffold->$contig_method
    #(
    #    assembly_id => $assembly->id,
    #    contig_num => $contig_num
    #);
    my $contig = $scaffold->get_contig( $contig_num);
    ######################
    $contig->set_chromosome($chromosome, $self->ordered);
    if ( $self->create_assembly )
    {
        $contig->update_attribute('length', $agp->{component_end});
    }
    $contig->update;

    return $contig;
}

1;

=pod

=head1 Name

Finishing::Assembly::AGP::DBIxImporter

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

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Finishing/Assembly/AGP/Importer.pm $
#$Id: Importer.pm 31152 2007-12-18 23:18:16Z ebelter $
