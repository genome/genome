package Finishing::Assembly::Consed::Schema;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use Finishing::Assembly::Ace::Schema;
use Finishing::Assembly::Phd::Ball;
use Finishing::Assembly::Phd::Directory;
use Finishing::Assembly::Phd::FastaAndQualDB;
use Finishing::Assembly::Consed::AssembledRead;

my %ace_schema :name(ace_schema:r) :isa(object);
my %phd_schemas :name(phd_schemas:r) :ds(aryref) :isa(object);

sub connect
{
    my ($class, $sources) = @_;

    my $acefile = delete $sources->{acefile}
        or $class->fatal_msg("Need acefile");
    
    my @types;
    if ( my $priority = delete $sources->{priority} )
    {
        @types = @$priority;
    }
    else
    {
        @types = keys %$sources;
    }

    my %types_and_schemas = 
    (
        directory => 'Directory',
        ball => 'Ball',
        fasta_and_qual_db => 'FastaAndQualDB',
    );
    
    my @phd_schemas;
    foreach my $type ( @types )
    {
        my $phd_schema_class = 'Finishing::Assembly::Phd::' . $types_and_schemas{$type};
        push @phd_schemas, $phd_schema_class->connect($sources->{$type});
    }

    return $class->new
    (
        ace_schema => Finishing::Assembly::Ace::Schema->connect($acefile),
        phd_schemas => \@phd_schemas,
    );
}

sub disconnect
{
    return 1;
}

sub get_assembly
{
    my $self = shift;

    my $assembly = $self->ace_schema->get_assembly;

    my $orig_assembled_reads = $assembly->{_assembled_reads};
    $assembly->{_assembled_reads} = sub
    {
        my $reads = $orig_assembled_reads->();
        my @reads;
        while ( my $read = $reads->next )
        {
            push @reads, Finishing::Assembly::Consed::AssembledRead->new
            (
                ace_source => $read,
                phd_source => sub{ return $self->_get_phd( $read->phd_file ); },
            );
        }

        return Finfo::Iterator->new(objects => \@reads);
    };

    my $orig_get_assembled_read = $assembly->{_get_assembled_read};
    $assembly->{_get_assembled_read} = sub
    {
        my $read = $orig_get_assembled_read->($_[0]);
        return Finishing::Assembly::Consed::AssembledRead->new
        (
            ace_source => $read,
            phd_source => sub{ return $self->_get_phd( $read->phd_file ); },
        );
    };

    return $assembly;
}

sub _get_phd : PRIVATE
{
    my ($self, $file_name) = @_;

    foreach my $schema ( @{ $self->phd_schemas } )
    {
        my $phd = $schema->get_phd($file_name);
        return $phd if $phd;
    }

    return;
}

1;

=pod

=head1 Name

Finishing::Assembly::Consed::Schema

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

