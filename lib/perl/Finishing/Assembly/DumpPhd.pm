package Finishing::Assembly::DumpPhd;

use strict;
use warnings;
no warnings 'reserved';
use Data::Dumper;

use Finfo::Std;

use File::Slurp;

use IO::File;
use Bio::SeqIO;

use GSC::IO::Assembly::PhdDB;
use GSC::IO::Assembly::Phd::Writer;
use GSC::IO::Assembly::Ace;

use GSCApp;

#attributes

my %base_dir :name(base_dir:o) :isa(dir_rw) :default('/gscmnt/temp113/finishing/scratch/macaque_merge/test');
my %phd_dir :name(phd_dir:o) :isa(dir_w);
my %grab_dir :name(grab_dir:o) :isa(dir_r);

sub START{
    my $self = shift;
    $self->phd_dir($self->base_dir.'/phd_dir') unless $self->phd_dir;
}

sub get_phd{
    my ($self, $phd_name)=@_;
    my $phd_path = $self->phd_dir."/$phd_name";
    $phd_path = $self->grab_dir."/$phd_name" unless -e $phd_path;
    die "Can't find $phd_name in ".$self->phd_dir." or ".$self->grab_dir unless -e $phd_path;
    my $fh = IO::File->new("< $phd_path");
    my $reader= Finishing::Assembly::Phd::Reader->new();
    return $reader->read($fh);
}

sub dump_phd{
    my ($self, %p) = @_;
    my $read_list = delete $p{'reads'};
    my $contig = delete $p{'contig'};# or $self->error_msg("Can't dump without contig");
	my $ace_contig = delete $p{'ace_contig'};
    my $ao;
    if (defined $contig){
        my $acefile = $contig->acefile;
        $self->info_msg("Acefile: $acefile");
        $ao=GSC::IO::Assembly::Ace->new(input_file => $acefile);
        $ace_contig = $ao->get_contig($contig->id);
    }
    my %reads = %{$ace_contig->reads};
    my @dump_list;
    if ($read_list){
        @dump_list = @$read_list;
    }else{
        @dump_list = keys %reads;
    }
    foreach my $read_name (@dump_list) {
        my $read = $reads{$read_name};
        my $seq = $read->sequence->unpadded_base_string;
        my $length = length $seq;
        my $phd_name = $read->ace_read->{'description'}->{'PHD_FILE'};
        my $time = $read->ace_read->{'description'}->{'TIME'};

        #my $phd_name  = $fasta->display_id;
        #my $time      = $fasta->desc;
        #my $seq       = $fasta->seq;
        #my $length    = $fasta->length;

        my $phd_file  = $self->phd_dir."/$phd_name";
        #print "\n$phd_name";

        if (-s $phd_file) {
            #print "$phd_file already exist\n";
            print ".";
            next;
        }

        my $pdo = Finishing::Assembly::PhdDB->new;

        my $po;
        eval {$po = $pdo->get_phd($phd_name)};

        if ($po) {
            my $bc;
            eval {$bc = $po->base_count};

            unless ($bc) {
                #print "$phd_name has no content in DB, now try makeup\n";
                $self->makeup($phd_name, $time, $seq, $length);
                next;
            }

            unless ($length == $bc) {
                #print "$phd_name from ace and from DB have different seq length\n";
                next;
            }
            my $phd_time = $po->comments->{TIME};
            unless ($time eq $phd_time) {
                $po->comments->{TIME} = $time;
                #print "$phd_name uses the time stamp in acefile\n";
            }
            unless ($seq eq $po->sequence->unpadded_base_string) {
                $po->sequence->unpadded_base_string($seq);
                #print "$phd_name uses the sequence in acefile\n";
            }
            my $fh  = IO::File->new("> $phd_file") 
                or die "can't write to $phd_file\n";
            my $phd = Finishing::Assembly::Phd::Writer->new();
            $phd->write($fh, $po);
            $fh->close;
        }
        else {
            #print "$phd_name has no content in DB\n";
            $self->makeup($phd_name, $time, $seq, $length);
        }
    }
}

sub makeup {
    my ($self, $phd_name, $time, $seq, $length) = @_;
    my ($si, $si_len);
    my $phd_file = $self->phd_dir."/$phd_name";

    $time = 'TIME: '.$time;

    my $name = $phd_name;
    $name =~ s/\.phd\./-/;

    eval {$si = GSC::Sequence::Item->get(sequence_item_name => $name)};

    unless ($si) {
        print "\ncan't makeup from db for $phd_name : no sequence item for $name\nTrying to generate from ace";
        $self->makeup_from_ace($phd_name, $time, $seq, $length);
        return;
    }

    eval {$si_len = $si->seq_length};
    unless ($si_len) {
        print "\ncan't makeup from db for $phd_name : no sequence_item length\nTrying to generate from ace";
        $self->makeup_from_ace($phd_name, $time, $seq, $length);
        return;
    }

    unless ($length == $si_len) {
        print "makeup $phd_name has diff seq length from the one in acefile\nTrying to generate from ace";
        $self->makeup_from_ace($phd_name, $time, $seq, $length);
        return;
    }

    require Bio::Seq::Quality;

    $name =~ s/\-\d+$//;

    my $sq = Bio::Seq::Quality->new(
        -seq  => $seq, 
        -qual => $si->sequence_quality_string,
        -id   => $name,			    
    );

    my $tmp = '/tmp/'.$phd_name;

    my $phd_io = Bio::SeqIO->new(
        -format => 'phd', 
        -file   => ">$tmp",
    );
    $phd_io->write_seq(-Quality => $sq);
    #is ok phd file contain upper case base seq ?

    my $text = read_file($tmp);
    $text =~ s/CHROMAT_FILE:.*?\n/CHROMAT_FILE: $name\n/;
    $text =~ s/TIME:.*?\n/$time\n/;
    write_file($phd_file, $text);

    unlink $tmp;

    #print "$phd_name is successfully made\n";
    return;
}

sub makeup_from_ace{
    my ($self, $phd_name, $time, $seq, $length) = @_;
    
    my $phd_file = $self->phd_dir."/$phd_name";

    $time = 'TIME: '.$time;

    my $name = $phd_name;
    $name =~ s/\.phd\./-/;

    my @qual;
    for (1..$length){
        push @qual, 15;
    }
        
    require Bio::Seq::Quality;

    $name =~ s/\-\d+$//;

    my $sq = Bio::Seq::Quality->new(
        -seq  => $seq, 
        -qual => \@qual,
        -id   => $name,			    
    );

    my $tmp = '/tmp/'.$phd_name;

    my $phd_io = Bio::SeqIO->new(
        -format => 'phd', 
        -file   => ">$tmp",
    );
    $phd_io->write_seq(-Quality => $sq);
    #is ok phd file contain upper case base seq ?

    my $text = read_file($tmp);
    $text =~ s/CHROMAT_FILE:.*?\n/CHROMAT_FILE: $name\n/;
    $text =~ s/TIME:.*?\n/$time\n/;
    write_file($phd_file, $text);

    unlink $tmp;

    return;
}
1;
