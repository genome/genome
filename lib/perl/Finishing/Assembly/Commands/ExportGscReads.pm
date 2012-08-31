package Finishing::Assembly::Commands::ExportGscReads;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use IO::File;
use NCBI::TraceArchive::Trace;

require Compress::Zlib;
require GSCApp; 

my %read_names :name(read_names:r) :ds(aryref);
my %ncbi_name  :name(ncbi_name:o) :isa(boolean) :default(1);
my %chromat_dir :name(chromat_dir:o) :isa(dir_rw);
my %phd_dir :name(phd_dir:o) :isa(dir_rw);
my %missed :name(_missed_names:p) :ds(aryref) empty_ok(1);

sub START
{
    my $self = shift;

    $self->fatal_msg() unless $self->phd_dir or $self->chromat_dir;

    App->init unless App::Init->initialized;
    
    return 1;
}

sub execute
{
    my $self = shift;

    my @missed_names;
    my %read_names;

    if ($self->ncbi_name) {
        for my $name ( @{ $self->read_names } ) {
            my $trace = NCBI::TraceArchive::Trace->new(name => $name);
            $read_names{ $trace->ncbi_name } = 1;
        }
    }
    else {
        map{$read_names{$_} = 1}@{$self->read_names};
    }

    my @reads = GSC::Sequence::Read->get
    (
        trace_name => [ keys %read_names ],
    );

    my $chromat_dir = $self->chromat_dir;
    my $phd_dir = $self->phd_dir;

    foreach my $read ( @reads )
    {
        my $name = $read->trace_name;
        my $scf_file;
        if ( $chromat_dir )
        {
            my $scf_name = $read->default_file_name('scf');
            $scf_name =~ s/^WIBR\-//g;
            my $scf_file = sprintf('%s/%s.gz', $chromat_dir, $scf_name);

            unlink $scf_file if -e $scf_file;

            if ( my $scf_content = $read->scf_content )
            {
                my $scf_fh = IO::File->new("> $scf_file");
                $self->fatal_msg("Can't open scf ($scf_file)\n$!")
                    and return unless $scf_fh;
                $scf_fh->print( Compress::Zlib::memGzip($scf_content) );
                $scf_fh->close;
            }
        }

        if ( $phd_dir )
        {
            my @read_edits = ( $read );
            if ( $read->isa('GSC::Sequence::ReadEdit') )
            {
                push @read_edits, $read->get_read, $read->get_previous_edits;
            }

            foreach my $read_edit ( @read_edits )
            {
                my $phd_name = $read_edit->default_file_name('phd');
                $phd_name =~ s/^WIBR\-//g;
                my $phd_file = sprintf('%s/%s', $phd_dir, $phd_name);

                unlink $phd_file if -e $phd_file;

                if ( my $phd_content = $read_edit->phd_content )
                {
                    my $phd_fh = IO::File->new("> $phd_file");
                    $self->fatal_msg("Can't open phd ($phd_file):\n$!") unless $phd_fh;
                    $phd_fh->print($phd_content);
                    $phd_fh->close;
                }
                elsif ( $scf_file and -s $scf_file )
                {
                    system "phred $scf_file -pd $phd_dir -process_nomatch -nocall";
                }
            }
        }

        delete $read_names{$name};
    }

    $self->_missed_names([ keys %read_names ]);

    return 1;
}

sub missed_names
{
    my $self = shift;

    return $self->_missed_names;
}

1;

=pod

=head1 Name

Finishing::Assembly::GSC::ReadExporter

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

