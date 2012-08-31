package Finishing::Assembly::Commands::SyncAcePhds;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use Finishing::Assembly::Ace::Exporter;
use Finishing::Assembly::Ace::Schema;
use Finishing::Assembly::Phd::Directory;

require File::Copy;

my %file :name(file:r) :isa(file_r);
my %phd_dir :name(phd_dir:r) :isa(dir_r);
my %missed :name(_missed_list:p) :ds(aryref) empty_ok(1);


sub execute
{
    my $self = shift;
    
    my @list;
    
    my $phd_schema = Finishing::Assembly::Phd::Directory->connect($self->phd_dir);
    
    my $tmp_ace = $self->file . '.tmp';
    unlink $tmp_ace if -e $tmp_ace;
    my $xporter = Finishing::Assembly::Ace::Exporter->new(file => $tmp_ace);

    my $ace_schema = Finishing::Assembly::Ace::Schema->connect($self->file);
    my $assembly = $ace_schema->get_assembly;
    my $contigs = $assembly->contigs;

    while ( my $contig = $contigs->next )
    {
        my $reads = $contig->assembled_reads;
        while ( my $read = $reads->next )
        {
            my $phd = $phd_schema->get_phd( $read->phd_file );
            unless ($phd) {
                $self->warn_msg("no phd: ".$read->phd_file); #error?
                push @list, $read->phd_file; 
                next;
            }
            #next if $read->time eq $phd->time;
            $read->time( $phd->time );
        }
        $xporter->export_contig(contig => $contig);
    }
    
    $xporter->close;
    
    $ace_schema->disconnect;    
    $self->_missed_list(\@list);
    
    rename $tmp_ace, $self->file;
    #unlink $self->file;
    #File::Copy::copy($tmp_ace, $self->file);
    #unlink $tmp_ace;

    return 1;
}


sub missed_names {
    return shift->_missed_list;
}

1;

=pod

=head1 Name

Finishing::Assembly::Ace::SyncPhdTimes

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

