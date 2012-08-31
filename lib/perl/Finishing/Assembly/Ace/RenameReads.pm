package Finishing::Assembly::Ace::RenameReads;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use Finishing::Assembly::Ace::Exporter;
use Finishing::Assembly::Ace::Schema;
use NCBI::TraceArchive::Trace;

require File::Copy;

my %file :name(file:r) :isa(file_r);
#my %n2g :name(convert_ncbi_to_gsc:o) :isa(boolean) :default(1);
#my %rm_scf :name(remove_scf:o) :isa(boolean) :default(1);

sub execute
{
    my $self = shift;

    my $tmp_ace = $self->file . '.tmp';
    unlink $tmp_ace if -e $tmp_ace;
    my $xporter = Finishing::Assembly::Ace::Exporter->new(file => $tmp_ace);

    my $ace_schema = Finishing::Assembly::Ace::Schema->connect($self->file);
    my $assembly = $ace_schema->get_assembly;
    my $contigs = $assembly->contigs;

    while ( my $contig = $contigs->next )
    {
        my $base_segments = $contig->base_segments;
        foreach my $bs ( @$base_segments )
        {
            my $trace = NCBI::TraceArchive::Trace->new(name => $bs->{read_name});
            $bs->{read_name} = $trace->gsc_name;
        }
        $contig->base_segments($base_segments);

        my $reads = $contig->assembled_reads;
        while ( my $read = $reads->next )
        {
            my $name = $read->name;
            my $trace = NCBI::TraceArchive::Trace->new(name => $name);
            my $gsc_name = $trace->gsc_name;
            next if $gsc_name eq $name;
            my $new_read = $read->rename
            (
                new_name => $gsc_name,
                update_phd => 1,
                update_chromat => 1,
            );
        }
        $xporter->export_contig(contig => $contig);
    }
    
    $xporter->close;
    
    $ace_schema->disconnect;
    
    #unlink $self->file;
    #File::Copy::copy($tmp_ace, $self->file);
    #unlink $tmp_ace;

    return 1;
}

1;

=pod

=head1 Name

Finishing::Assembly::Ace::RenameReads

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

