package Finishing::Assembly::Commands::MakeChromosomeBlastdb;

use strict;
use warnings;
no warnings 'reserved';


use base 'Finishing::Assembly::Commands::Base';

use Bio::SeqIO;

my %asm_name :name(asm_name:r) :isa(string);
my %org_name :name(org_name:r) :isa(string);
my %chr_name :name(chr_name:r) :isa(string);
my %out_name :name(out_name:o) :isa(string);
my %out_dir  :name(out_dir:o)  :isa(dir_rw);


sub execute {
    my $self = shift;
    
    my $org = $self->_factory->get_organism($self->org_name);
    my $asm = $org->get_assembly($self->asm_name);
    
    my $chr = $org->get_chromosome($self->chr_name);
    $self->fatal_msg('can not get chr.'.$self->chr_name.' from DB') unless $chr;
    
    my $dir = $self->out_dir || '.';
    chop $dir if $dir =~ /\/$/;

    my $out_name  = $self->out_name || 'chr'.$self->chr_name.'_contigs';
    my $out_file = "$dir/$out_name";
    
    my $out_io = Bio::SeqIO->new(-format => 'fasta', -file => ">$out_file");
    my $ctg_itr = $asm->contigs(chromosome_id => $chr->id);

    my $ct = 0;
    while (my $ctg = $ctg_itr->next) {
        $ct++;
        print $ct.'...' if $ct % 1000 == 0;
        $out_io->write_seq(Bio::Seq->new(-seq=>$ctg->unpadded_base_string, -id => $ctg->name));
    }
    print "$ct\n";

    $self->fatal_msg("$out_file not existing") unless -s $out_file;

    my $ec = system "xdformat -n -I $out_file";
    $self->fatal_msg("xdformat failed on $out_file") if $ec;

    return 1;
}

1;

=pod

=head1 Name

Finishing::Assembly::Commands::MakeChromosomeBlastdb

=head1 Description

This module is to make chromosome contig blast DB out of database  

=head1 Synopsis

$mkdb = Finishing::Assembly::Commands::MakeChromosomeBlastdb->new(
    db       => 'dwrac_asm',
    org_name => 'pan_troglodytes',
    asm_name => '3.0',
    chr_name => '21',
);

$mkdb->execute;

=head1 Disclaimer

Copyright (C) 2008 Washington University Genome Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Feiyu Du> I<fdu@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$



