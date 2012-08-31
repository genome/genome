package Genome::Model::Tools::Fasta::HugoToFasta;

use strict;
use warnings;

use Genome;

use Bio::Seq;
use Bio::DB::Fasta;

class Genome::Model::Tools::Fasta::HugoToFasta {
    is => 'Command',                       
    has => [ 
	organism => {
	    type  =>  'String',
	    doc   =>  "provide the organism either mouse or human; default is human",
	    is_optional  => 1,
	    valid_values => ['human','mouse'],
	    default => 'human',
	},
	version => {
	    type  =>  'String',
	    doc   =>  "provide the imported annotation version; default for human is 54_36p_v2 and for mouse is 54_37g_v2",
	    is_optional  => 1,
	    default => '54_36p_v2',
	},
	gene => {
	    type  =>  'String',
	    doc   =>  "input a hugo gene name",
	},
	provide_utr => {
	    is =>  'Boolean',
	    doc   =>  "this will write the utr-sequence to the fasta file",
	    is_optional  => 1,
	},
	write_all_transcripts => {
	    is =>  'Boolean',
	    doc   =>  "this will write a fasta for all transcripts",
	    is_optional  => 1,
	},
	include_protein_fasta => {
	    is =>  'Boolean',
	    doc   =>  "this will write a protein fasta",
	    is_optional  => 1,
	},
     ], 
};

sub help_synopsis {
    return <<EOS

gmt fasta hugo-to-fasta -h

EOS
}

sub help_detail {
    return <<EOS 

	This tool will return the transcript fasta with the longest coding sequence for your hugo gene name. Currently, gene transcript relations are only seen from version 54_36p_v2 from human. 

EOS
}

sub execute {

    my $self = shift;
    
    my $gene = $self->gene;
    my $organism = $self->organism;
    my $version = $self->version;

    my $transcripts = &get_transcripts($gene);
    my @transcripts = split(/\,/,$transcripts);

    my $fastas;

    my @results;
    for my $transcript (@transcripts) {
	
	my $output = "$gene.$transcript";

	my ($info) = Genome::Model::Tools::Annotate::TranscriptSequence->execute(transcript => $transcript, no_stdout => '1', fasta => '1', organism => $organism, version => $version);
	unless ($info) {$self->error_message("Couldn't get info for $transcript"); next;}
	my $transcript_info = $info->{result};
	unless ($transcript_info) {$self->error_message("Couldn't get transcript info for $transcript"); next;}

	my $strand = $transcript_info->{$transcript}->{-1}->{strand};
	my $sequence = $transcript_info->{$transcript}->{-1}->{sequence};
	unless ($sequence) {$self->error_message("Couldn't find a the sequence for $transcript"); next;}
	my $protien_seq = $transcript_info->{$transcript}->{-1}->{protien_seq};
	my $length = length($sequence);
	my $fasta = $transcript_info->{$transcript}->{-1}->{FASTA};
	my $chromosome = $transcript_info->{$transcript}->{-1}->{chromosome};
	unless ($protien_seq) {$self->error_message("Couldn't find a protein sequence for $transcript"); next;}

	unless ($protien_seq =~ /^M/ && $protien_seq =~ /\*$/) {
	    $length = -1/$length;
	}
	
	my $protein_fasta;
	if ($self->include_protein_fasta) {
	    $protein_fasta = "\>$output.protein.fasta\n$protien_seq\n";
	}
	
	if ($self->write_all_transcripts) {

	    if ($self->provide_utr) {
		$fasta = &get_utr_fasta($self,$transcript_info,$transcript);
	    }

	    open(OUT,">$output.fasta") || $self->error_message("Could not open and write to $output.fasta") && next;
	    print OUT qq($fasta);
	    close OUT;
	    push(@results,"$output.fasta");
	    if ($self->include_protein_fasta) {
		open(OUT,">$output.protein.fasta") || $self->error_message("Could not open and write to $output.fasta") && next;
		print OUT qq($protein_fasta);
		close OUT;
		push(@results,"$output.protein.fasta");
	    }

	} else {
	    $fastas->{$length}->{$transcript}->{fasta} = $fasta;
	    if ($self->provide_utr) {
		$fastas->{$length}->{$transcript}->{transcript_info} = $transcript_info;
	    }
	    if ($self->include_protein_fasta) {
		$fastas->{$length}->{$transcript}->{protien_seq} = $protein_fasta;
	    }
	}
    }

    unless ($self->write_all_transcripts) {
	my $n = 1;
	foreach my $length ( sort {$b<=>$a} keys %{$fastas}) {
	    next if $n > 1;
	    foreach my $transcript (sort {$b gt $a} keys %{$fastas->{$length}}) {
		my $output = "$gene.$transcript";
		my $fasta = $fastas->{$length}->{$transcript}->{fasta};

		next if $n > 1;

		if ($self->provide_utr) {
		    my $transcript_info = $fastas->{$length}->{$transcript}->{transcript_info};
		    $fasta = &get_utr_fasta($self,$transcript_info,$transcript);
		}
		
		open(OUT,">$output.fasta") || $self->error_message("Could not open and write to $output.fasta") && next;
		print OUT qq($fasta);
		close OUT;
		push(@results,"$output.fasta");

		$n++;
		
		if ($self->include_protein_fasta) {
		    my $protein_fasta = $fastas->{$length}->{$transcript}->{protien_seq};
		    open(OUT,">$output.protein.fasta") || $self->error_message("Could not open and write to $output.fasta") && next;
		    print OUT qq($protein_fasta);
		    close OUT;
		    push(@results,"$output.protein.fasta");

		}
	    }
	}
    }
    unless (@results == 1) {
	my $last = pop(@results);
	push(@results,"and");
	push(@results,$last);
    }

    my $result = join ', ' , @results;
    $result =~ s/, and,/ and/;
    print qq(Your results have been writen to $result.\n);
    return @results;
}


sub get_transcripts {

    my ($gene) = @_;

    # bdericks This might need to be changed to point at new 54_36p_v2 version
    my $GTDir = "/gscmnt/200/medseq/biodb/shared/misc/annotation/54_36p";

    my $gtdb = Bio::DB::Fasta->new($GTDir);
    my $transcripts = $gtdb->seq($gene, 1 => 1000);

    return unless $transcripts;
    return $transcripts;

}

sub get_utr_fasta {

    my ($self,$transcript_info,$transcript) = @_;
    my $organism = $self->organism;
    my $gene = $self->gene;

    my $strand = $transcript_info->{$transcript}->{-1}->{strand};
    my $chromosome = $transcript_info->{$transcript}->{-1}->{chromosome};
    
    my (@positions);
    
    if ($strand eq "+1") {
	foreach my $pos (sort {$a<=>$b} keys %{$transcript_info->{$transcript}}) {
	    unless ($pos == -1) {
		push(@positions,$pos);
	    }
	}
    } elsif ($strand eq "-1") {
	foreach my $pos (sort {$b<=>$a} keys %{$transcript_info->{$transcript}}) {
	    unless ($pos == -1) {
		push(@positions,$pos);
	    }
	}
    } else {
	die "strand not define for transcript $transcript\n";
    }
    
    my @seq;
    for my $pos (@positions) {
	
	my $frame = $transcript_info->{$transcript}->{$pos}->{frame};
	my $aa_n = $transcript_info->{$transcript}->{$pos}->{aa_n};
	my $base = $transcript_info->{$transcript}->{$pos}->{base};
	my $masked_base = $transcript_info->{$transcript}->{$pos}->{masked_base};
	my $trans_pos = $transcript_info->{$transcript}->{$pos}->{trans_pos};
	
	my $range = $transcript_info->{$transcript}->{$pos}->{range};
	my ($r_start,$r_stop) = split(/\-/,$range);
	
	if ($base =~ /(\d+)\:\S/) {
	    my $coord = $1; 
	    $base = &get_utr_seq($coord,$strand,$organism,$chromosome);
	}
	push(@seq,$base);
    }
    my $sequence = join '' , @seq;
    my $fasta = "\>$transcript\|$strand\|Chromosome:$chromosome\|$gene\n$sequence\n";
   return $fasta; 
}

sub get_utr_seq {
    
    my ($pos,$strand,$organism,$chromosome) = @_;
    my $base;
    if ($strand eq "-1") {
	my $seq = &get_ref_base($pos,$pos,$chromosome,$organism);
	my $rev = &reverse_complement_allele($seq);
	$base = $rev;
    } else {
	my $seq = &get_ref_base($pos,$pos,$chromosome,$organism);
	$base = $seq;
    }
    return($base);
}

sub reverse_complement_allele {
    my ($allele_in) = @_;

    unless ($allele_in =~ /[ACGT]/) { return ($allele_in); }
    my $seq_1 = new Bio::Seq(-seq => $allele_in);
    my $revseq_1 = $seq_1->revcom();
    my $rev1 = $revseq_1->seq;
    
    my $r_base = $rev1;
    
    return $r_base;
}

sub get_ref_base {

    my ($chr_start,$chr_stop,$chr_name,$organism) = @_;
    my $RefDir;
    if ($organism eq "human"){
	$RefDir = "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/";
    } else {
	$RefDir = "/gscmnt/sata147/info/medseq/rmeyer/resources/MouseB37/";
    }
    my $refdb = Bio::DB::Fasta->new($RefDir);
    my $seq = $refdb->seq($chr_name, $chr_start => $chr_stop);
    $seq =~ s/([\S]+)/\U$1/;

    return $seq;
}

1;
