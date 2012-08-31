package Finishing::Assembly::Commands::ExportReads;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use Finishing::Assembly::Commands::ExportGscReads;
use Finishing::Assembly::Phd::Exporter;
use Finishing::Assembly::Phd::FastaAndQualDB;
use NCBI::TraceArchive;
use NCBI::TraceArchive::ProcessTraces;

my %read_names :name(read_names:r) :ds(aryref);
my %chromat_dir :name(chromat_dir:r) :isa(dir_rw);
my %phd_dir :name(phd_dir:r) :isa(dir_rw);
my %pri_srcs :name(read_sources:o)
    :ds(aryref)
    :isa([ 'in_list', __PACKAGE__->valid_read_sources ])
    :default([ __PACKAGE__->valid_read_sources ]);
my %missed :name(_missed_names:p) :ds(aryref) empty_ok(1);

sub valid_read_sources
{
    return (qw/ chimp gsc trace_archive /);
}

sub execute
{
    my $self = shift;

    my $read_names = $self->read_names;

    foreach my $src ( @{ $self->read_sources } )
    {
        my $method = '_' . $src;
        $read_names = $self->$method($read_names);
        last unless @$read_names;
    }

    $self->_missed_names($read_names);
    
    return 1;
}

sub missed_names
{
    my $self = shift;

    return $self->_missed_names
}

sub _gsc
{
    my ($self, $read_names) = @_;

    my $read_xporter = Finishing::Assembly::Commands::ExportGscReads->new
    (
        read_names => $read_names,
        chromat_dir => $self->chromat_dir,
        phd_dir => $self->phd_dir,
    );

    $read_xporter->execute;

    return $read_xporter->missed_names;
}

sub _trace_archive
{
    my ($self, $read_names) = @_;

    my $trace_archive = NCBI::TraceArchive->new
    (
        dir => $self->chromat_dir,
        src_names => $read_names,
    );

    $trace_archive->execute;

    my $retrieved_traces = $trace_archive->retrieved_traces;
    if ( @$retrieved_traces )
    {
        my $processor =  NCBI::TraceArchive::ProcessTraces->new
        (
            traces => $retrieved_traces,
            trace_location => $self->chromat_dir,
            chromat_dir => $self->chromat_dir,
            phd_dir => $self->phd_dir,
        );
        $processor->execute;
    }

    return $trace_archive->missed_gsc_names;
}

sub _chimp
{
    my ($self, $read_names) = @_;

    my (@gotten_names, @missed_names);
    
    my $chromat_dir = $self->chromat_dir;
    my $phd_dir = $self->phd_dir;

    my $factory = Finishing::Assembly::Factory->connect
    (
        'phd_fnqdb',
        '/gscmnt/temp113/finishing/scratch/chimp/2.1_051011/clipped_fasta_and_qual/chimp_read_fa_qual_index.sqlite'
    );

    foreach my $name ( @$read_names )
    {
        my $phd = $factory->get_phd($name);

        push @missed_names, $name
            and next unless $phd;

        my $xporter = Finishing::Assembly::Phd::Exporter->new
        (
            read => $phd,
            file => sprintf('%s/%s', $phd_dir, $phd->phd_file),
        );
        $xporter->execute;

        push @gotten_names, $name;
    }
    
    if ( @gotten_names )
    {
        my $trace_archive = NCBI::TraceArchive->new
        (
            dir => $chromat_dir,
            src_names => \@gotten_names,
            leave_scf_ext => 1,
        );

        $trace_archive->execute;

        foreach my $name ( @{ $trace_archive->retrieved_ncbi_names } )
        {
            my $trace = sprintf('%s/%s.scf', $chromat_dir, $name);
            system "gzip -f $trace";       
        }
    }

    return \@missed_names;
}

1;

=pod

=head1 Name

Finishing::Assembly::Commands::GetScfsAndPhds

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

