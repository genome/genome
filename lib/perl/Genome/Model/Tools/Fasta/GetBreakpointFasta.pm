package Genome::Model::Tools::Fasta::GetBreakpointFasta;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Fasta::GetBreakpointFasta {
    is => 'Command',                    
    has => [ # specify the command's properties (parameters) <--- 
	     

	     bp_ids => {
		 type  =>  'String',
		 doc   =>  "provide a breakpoint id or several ids each seperated by a comma",
	     },
	     name => {
		 type  =>  'String',
		 doc   =>  "Give a base name for your output file. Defualt will use first breakpoint id",
		 is_optional  => 1,
	     },
	     flank => {
		 type  =>  'Number',
		 doc   =>  "give the number of bases you would like on either side of your breakpoints, default set at 500",
		 is_optional  => 1,
	     },

	     bases_per_line => {
		 type  =>  'Number',
		 doc   =>  "give the number of bases you wolud like on each line of the fasta, default set at 50",
		 is_optional  => 1,
	     },

	     organism => {
		 type  =>  'String',
		 doc   =>  "provide the organism either mouse or human; default is human",
		 is_optional  => 1,
		 default => 'human',
	     },
	],
	
    
};

sub help_brief {
    return <<EOS
  This tool was design to retrieve sequence from NCBI Human Build 36, sequence from an individual or list of breakpoint ids.
EOS
}

sub help_synopsis {
    return <<EOS

running with optional minimum input...

gmt fasta get-breakpoint-fasta --bp-ids 
 
...will provide you with a fasta file of sequence from your breakpoint_ids

EOS
}

sub help_detail {
    return <<EOS 

running...
gmt fasta get-breakpoint-fasta --bp-ids chr3:103057567-103057567

    will return the result as a file with the default name
       chr3:103057567-103057567.1.500.fasta

running...
gmt fasta get-breakpoint-fasta --bp-ids chr3:103057567-103057567,chr3:103058567-103058567
    will return the result as a file with the default name
       chr3:103057567-103057567.2.500.fasta

if you use the --name option .2.500.fasta will be appended to it were the 2 is the number of breakpoints you are retrieving sequence for and the 500 is the length of flank retrieved.  --flank is optional with the default set at 500.


EOS
}


sub execute {

    my $self = shift;
    
    my $bp_ids = $self->bp_ids;
    my $flank = $self->flank;
    unless ($flank) {$flank = 500;}
    my $name = $self->name;
    my $bases_per_line = $self->bases_per_line;
    unless($bases_per_line) {$bases_per_line = 50;}

    my @bps = split(/\,/,$bp_ids);
    my $bp_n = @bps;

    my $organism = $self->organism;
    
    my $out;
    if ($name) {
	$out = "$name.$bp_n.$flank.fasta";
	open(OUT,">$out");
    }
    
    for my $bp_id (@bps) {
	
	unless ($out) {
	    $out="$bp_id.$bp_n.$flank.fasta";
	    open(OUT,">$out");
	}
	
	my ($chromosome,$breakpoint1,$breakpoint2);
	if ($bp_id =~ /chr([\S]+)\:(\d+)\-(\d+)/) { 
	    $chromosome = $1;
	    $breakpoint1 = $2;
	    $breakpoint2 = $3;
	} else { die "please check the format of your breakpoint id\n"; }
	
	my $start = $breakpoint1 - $flank;
	my $stop = $breakpoint2 + $flank;
	my ($mid) = &get_ref_base($chromosome,$breakpoint1,$breakpoint2,$organism);
	my ($lstop) = $breakpoint1 - 1;
	my ($rstart) = $breakpoint2 + 1;
	
	my ($lseq) = &get_ref_base($chromosome,$start,$lstop,$organism);
	my ($rseq) = &get_ref_base($chromosome,$rstart,$stop,$organism);
	
	my $fullseq = $lseq . $mid . $rseq;

	my $header;
	if ($organism eq "human") {
	    $header = "\>$bp_id.fasta Flank $flank NCBI Human Build 36, Chr:$chromosome, Coords $start-$stop, Ori (+)";
	} else {
	    $header = "\>$bp_id.fasta Flank $flank NCBI Mouse Build 37, Chr:$chromosome, Coords $start-$stop, Ori (+)";
	}

	print OUT qq($header\n);
	
	my @seq = split(//,$fullseq);
	my $n = 0;
	for my $base (@seq) {
	    $n++;
	    print OUT qq($base);
	    if ($n == $bases_per_line) {
		print OUT qq(\n);
		$n = 0;
	    }
	}
	print OUT qq(\n);
	
	
	#print OUT qq($lseq\[$mid\]$rseq\n);
	#print qq(see the result in fasta file $out\n);
	
	
    }
    
    print qq(see the result in fasta file $out\n);
}


sub get_ref_base {
    
    my ($chr_name,$chr_start,$chr_stop,$organism) = @_;

#used to generate the refseqs;
    use Bio::DB::Fasta;

    my $RefDir;
    if ($organism eq "human"){
	$RefDir = "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/";
    } else {
	$RefDir = "/gscmnt/sata147/info/medseq/rmeyer/resources/MouseB37/";
    }
    my $refdb = Bio::DB::Fasta->new($RefDir);
    
    my $seq = $refdb->seq($chr_name, $chr_start => $chr_stop);
    $seq =~ s/([\S]+)/\U$1/;
    
    if ($seq =~ /N/) {warn "your sequence has N in it\n";}
    
    return $seq;
    
}



