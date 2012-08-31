package Genome::Model::Tools::Sv::AluFilter;

use warnings;
use strict;

use Genome;
use IO::File;
use Carp;

class Genome::Model::Tools::Sv::AluFilter {
    is => 'Command',
    has_input => [
    sv_calls => {
        is => 'String',
        doc => 'path to merged SV call file you would like to filter',
    },
    output_file => {
        is => 'String',
        doc => 'a possibly reduced set of SV calls that are NOT in ALU regions',
    },
    ],
    doc => 'Remove SV DEL events in ALU regions from an input file.',
};

sub help_detail {
    return <<EOS
    This tool looks through files downloaded here: /gscmnt/sata194/info/sralign/UCSC/data/ to decide which type of region an SV DEL event is in. If the DEL is in an ALU sequence, it is not printed into the output file. If it is NOT in an ALU region, it will be printed into the output file, and thus preserved.
EOS
}

sub execute {

    #parse input arguments and open filehandles
    my $self = shift;
    my $sv_calls = $self->sv_calls;
    unless (-s $sv_calls) {
        $self->error_message("Input file in MERGED BreakDancer format required.");
        return;
    }
    my $input_fh = new IO::File $sv_calls,"r";
    my $output_file = $self->output_file;
    my $out_fh = new IO::File $output_file,"w" or die "Could not open output file $output_file for writing.";
    my @acceptable_sv_types = (qw/CTX DEL INS INV ITX/);

    #read through input file, writing to output file as you go
    #Format: #ID     CHR1    OUTER_START     INNER_START     CHR2    INNER_END       OUTER_END       TYPE    ORIENTATION     MINSIZE MAXSIZE SOURCE  SCORES  Copy_Number
    while (my $line = $input_fh->getline) {

        #retain headers
        if ($line =~ /^#/) { 
            print $out_fh $line; 
            next; 
        }

        #split line
        my ($id,$chr1,$bpA,undef,$chr2,$bpB,undef,$type) = split /\t/,$line;

        #check that split worked by checking for acceptable sv types
        unless (grep { /^$type$/ } @acceptable_sv_types) {
            $self->error_message("Unrecognized SV type $type.");
            return;
        }

        #only checking DEL events smaller than 10kbp
        unless ($type eq "DEL") {
            print $out_fh "$line";
            next;
        }
        if ($bpB-$bpA > 10000) {
            print $out_fh "$line";
            next;
        }

        #check appropriately sized DEL event for ALU region  
        my $is_alu = $self->is_alu_polymorphism($chr1, $bpA, $bpB);
        print $out_fh "$line" unless $is_alu;

    }

    return 1;
}

sub is_alu_polymorphism {
    my $self = shift;
    my ($chr, $bpA, $bpB) = @_;
    my $is_alu = 0;

    #open alu file
    my $dir = "/gscmnt/sata194/info/sralign/UCSC/data";
    my $file = "$dir/chr$chr"."_rmsk.txt";
    my $alu_fh = new IO::File $file,"r" or die "Could not open ALU file $file\n";

    while (my $line = $alu_fh->getline) {
	chomp $line;
	my (undef, undef, undef, undef, undef,  undef, $start, $stop, undef, undef, undef, $element) = split /\s+/, $line;
	next if ($element ne "LINE" && $element ne "SINE" && $element ne "LTR" && $element ne "DNA");
	last if ($stop > $bpB + 1000);
	if ( abs($start-$bpA) <= 50 && abs($stop-$bpB) <= 50 ) { 
            $is_alu = 1;
            last;
        }
    }
    return $is_alu;
}

1;
=cut
#Legacy JW Code
    if ( $bpB-$bpA > 10000 ) { next; }
    my ($element, $start, $stop) = aluPolymorphism($chr, $bpA, $bpB);
    if ( defined $element ) { 
	print "$line\t$element\t", abs($start-$bpA), "\t", abs($stop-$bpB), "\n";
