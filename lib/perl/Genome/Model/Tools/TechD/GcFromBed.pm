package Genome::Model::Tools::TechD::GcFromBed;

use strict;
use warnings;

use Genome;

#can't use this until on Perl 5.10
#use Bio::DB::Sam;

class Genome::Model::Tools::TechD::GcFromBed {
    is => 'Command',
    has => [
        bed_file => {
            is => 'String',
            doc => 'A BED format file of regions to interegate G+C content.',
        },
        fasta_file => {
            is => 'String',
            doc => 'The FASTA file that has been indexed with faidx and cooresponds to the bed file',
            is_optional => 1,
            default_value => '/gscmnt/gc4096/info/model_data/2741951221/build101947881/all_sequences.fa',
        },
        output_file => {
            is => 'String',
            doc => 'The output file which will contain the full line from the bed_file with G+C and DNA(optional) appended(tab delimited).',
        },
        sequence =>  {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
        },
    ],
};

sub execute {
    my $self = shift;
    unless ($] > 5.010) {
        die "Bio::DB::Sam requires perl 5.10 or greater!";
    }
    require Bio::DB::Sam;
    unless (-f $self->fasta_file .'.fai') {
        die('FASTA file '. $self->fasta_file .' has not beed indexed by faidx.');
    }
    my $fai = Bio::DB::Sam::Fai->load($self->fasta_file);
    my $output_fh = IO::File->new($self->output_file,'w');
    my $bed_fh = IO::File->new($self->bed_file,'r');
    while (my $line = $bed_fh->getline) {
        chomp($line);
        my ($chr,$start,$end) = split("\t",$line);
        my $seq = $fai->fetch($chr .':'. ($start+1) .'-'. $end);
        my $gc = get_gc_from_string($seq);
        print $output_fh $line ."\t". $gc;
        if ($self->sequence) {
            print $output_fh "\t". $seq;
        }
        print $output_fh "\n";
    }
    $output_fh->close;
    return 1;
}

sub get_gc_from_string {
    my $dna_string = shift;
    map {$_ = 0} my ($count, $A, $C, $G, $T);
    foreach my $nt (split ('', $dna_string)) {
        $A++ if ($nt =~ /A/i);
        $C++ if ($nt =~ /C/i);
        $G++ if ($nt =~ /G/i);
        $T++ if ($nt =~ /T/i);
        $count++;
    }
    my $gc = sprintf("%.02f",(($G + $C)  / $count));
    return $gc;
};


1;
