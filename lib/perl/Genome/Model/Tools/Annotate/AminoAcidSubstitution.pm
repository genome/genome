package Genome::Model::Tools::Annotate::AminoAcidSubstitution;

use strict;
use warnings;
use Genome;

use Bio::Seq;
use Bio::Tools::CodonTable;
use Bio::DB::Fasta;

class Genome::Model::Tools::Annotate::AminoAcidSubstitution {
    is  => 'Genome::Model::Tools::Annotate',
    has => [
        transcript => {
            type     => 'String',
 	    is_optional  => 1,
           doc      => "Provide the transcript name.",
        },
	gene => {
            type     => 'String',
	    is_optional  => 1,
            doc      => "Provide the hugo gene name.",
        },
	
	amino_acid_substitution => {
	    type      => 'String',
            is_input => 1,
	    doc       => "Provide the amino acid substitution or substitution's for the transcript represented in this format => P2249A or Q165H,S166C.",
	},
	 organism => {
	    type  =>  'String',
	    doc   =>  "Provide the organism either mouse or human; default is human.",
	    is_optional  => 1,
	    default => 'human',
	},
	  version => {
	    type  =>  'String',
	    doc   =>  "Provide the imported annotation version; default for human is 54_36p_v2 and for mouse is 54_37g_v2.",
	    is_optional  => 1,
	    default => '54_36p_v2',
	},
	   output => {
	    type  =>  'String',
	    doc   =>  "Provide a file name to write your results to. \".txt\" will be appended to it. Default is to print to stdout.",
	    is_optional  => 1,
	},
	
	],
};

sub help_synopsis { 
"
gmt annotate amino-acid-substitution -transcript ENST00000269305 -amino-acid-substitution S166C
        or
gmt annotate amino-acid-substitution -gene TP53 -amino-acid-substitution S166C
"
}

sub help_detail {
    return <<EOS

This tool was designed to identify all posible base changes in a codon that could produce a given amino acid substitution. It will also provide the frame and genomic coordinates of the bases that change.
	
EOS
}

sub execute { 
    
    my $self = shift;
    
    my $transcript = $self->transcript;
    my $gene = $self->gene;

    unless ($transcript || $gene) {
	$self->error_message("\nyou need to provide either a hugo gene name or a transcript name\n\n"); return;
    }

    my $amino_acid = $self->amino_acid_substitution;
    my $organism = $self->organism;
    my $version = $self->version;
    if ($organism eq "mouse" && $version eq "54_36p_v2") { $version = "54_37g_v2"; }

    my @transcripts;
    if ($gene) {
	my $transcripts = &get_transcripts($gene,$version);
	@transcripts = split(/\,/,$transcripts);
    } else {
	push(@transcripts,$transcript);
    }

    my $output = $self->output;
    if ($output) {
	open(OUT,">$output.txt") || $self->error_message("couldn't open the output file $output.txt") && return;
    }
    my @results;
    for my $transcript (@transcripts) {

	my $TranscriptSequence = Genome::Model::Tools::Annotate::TranscriptSequence->create(transcript => $transcript, organism => $organism, version => $version, no_stdout => "1");
	unless ($TranscriptSequence) { $self->error_message("couldn't create a transcript sequence object for $transcript"); next;}
	
	my ($transcript_info) = $TranscriptSequence->execute();
	unless ($transcript_info) { $self->error_message("couldn't execute the transcript sequence object for $transcript"); next;}
	
	my @positions = &get_positions ($transcript_info,$transcript);
	unless (@positions) { $self->error_message("couldn't extract positions from the transcript sequence object for $transcript"); next;}
	
	my @amino_acid_subs = split(/\,/,$amino_acid);
	for my $nsprotein (@amino_acid_subs) { #nsprotein nonsynonymous protein
	    
	    my ($taa,$protein_number,$daa) = $nsprotein =~ /^(\D)([\d]+)(\D)$/;
	    unless ($taa && $protein_number && $daa) {
		$self->error_message("\n$nsprotein is an invalid format. The amino acid change should be represented in this format => P2249A. $nsprotein will be skipped.\n\n");
		if ($output) {
		    print OUT qq(\n$nsprotein is an invalid format. The amino acid change should be represented in this format => P2249A. $nsprotein will be skipped.\n\n);
		} 
		next;
	    }
	    $taa =~ s/(\D)/\U$1/;
	    $daa =~ s/(\D)/\U$1/;
	    unless ($taa =~ /[C,H,I,M,S,V,A,G,L,P,T,R,F,Y,W,D,N,E,Q,K]/ && $daa =~ /[C,H,I,M,S,V,A,G,L,P,T,R,F,Y,W,D,N,E,Q,K]/) {
		$self->error_message("\n$nsprotein is an invalid format. The amino acids most be one of twenty found in a protein chain. $nsprotein will be skipped.\n\n");
		if ($output) {
		    print OUT qq(\n$nsprotein is an invalid format. The amino acids most be one of twenty found in a protein chain. $nsprotein will be skipped.\n\n);
		} 
		next;
	    }
	    
	    my ($p1,$p2,$p3,$b1,$b2,$b3) = &get_codon ($transcript_info,$transcript,$protein_number,@positions);
	    unless ($p1 && $p2 && $p3 && $b1 && $b2 && $b3) {  
		$self->error_message("\nCouldn't identify the target codon $protein_number in $transcript.  $nsprotein will be skipped for $transcript.\n\n");
		if ($output) {
		    print OUT qq(\nCouldn't identify the target codon $protein_number in $transcript.  $nsprotein will be skipped for $transcript.\n\n);
		} 
		next;
	    }
	    
	    my ($result) = &get_result($p1,$p2,$p3,$b1,$b2,$b3,$transcript,$taa,$protein_number,$daa);
	    unless ($result) {
		$self->error_message("\nNo result was found for $nsprotein in $transcript. $nsprotein will be skipped for $transcript.\n\n");
		if ($output) {
		    print OUT qq(\nNo result was found for $nsprotein in $transcript.  $nsprotein will be skipped for $transcript.\n\n);
		} 
		next;
	    }
	    
	    if ($output) {
		print OUT qq(\n$result\n\n);
	    } else {
		print qq(\n$result\n\n);
	    }
	    push(@results,$result);
	}
    }
    if ($output) {
	print qq(Your results have been printed in $output.txt\n);
	close OUT;
    }
    return unless @results;
    return 1;
    
}

sub get_result {

    my @result;

    my ($p1,$p2,$p3,$b1,$b2,$b3,$transcript,$taa,$protein_number,$daa) = @_;

    my $codon = "$b1$b2$b3";
    my $display_id = $transcript;
    my $newcodon = Bio::Seq->new( -display_id => $display_id, -seq => $codon );
    my $aa = $newcodon->translate->seq;
    
    unless ($taa eq $aa) {my $result = "The amino acid $taa input for position doesn\'t match the expected amino acid identified as $aa in $transcript"; return $result; }

    
    my $myCodonTable = Bio::Tools::CodonTable->new();
    my @pcodons = $myCodonTable->revtranslate($daa);
    
    my $all_combo_codons = join ' or ' , @pcodons;
    $all_combo_codons =~ s/([\S\s]+)/\U$1/;
    

    my $line = "A mutation in amino acid number $protein_number of $transcript causing a nonsynonymous change from $aa to $daa could occur by changing the codon $codon to $all_combo_codons.\n";
    push (@result,$line);
    

    for my $c (@pcodons) {
	$c =~ s/(\S+)/\U$1/;

	my @line;	

	my ($nb1,$nb2,$nb3) = split(//,$c);
	
	$line = "To get the new amino acid with the codon $c,";
	push (@line,$line);

	my ($d1,$d2,$d3);
	
	unless ($b1 eq $nb1) {
	    $line = "the first frame at $p1 would change from $b1 to $nb1";
	    push (@line,$line);
	    $d1 = 1;
	}
	unless ($b2 eq $nb2) {
	    if ($d1) {
		$line = "and";
		push (@line,$line);
	    }
	    $d2 = 1;
	    $line = "the second frame at $p2 would change from $b2 to $nb2";
	    push (@line,$line);
	}
	unless ($b3 eq $nb3) {
	    if ($d1 || $d2) {
		$line = "and";
		push (@line,$line);
	    }
	    $d3 = 1;
	    $line = "the thrid frame at $p3 would change from $b3 to $nb3";
	    push (@line,$line);
	}
	$line = join ' ' , @line;

	$line = "$line.";
	undef@line;
	push (@result,$line);
    }
    my $result = join "\n" , @result;
    return($result);
}

sub get_positions {

    my ($transcript_info,$transcript) = @_;
    my $strand = $transcript_info->{$transcript}->{-1}->{strand};
    my @positions;
    if ($strand eq "+1") {
	foreach my $pos (sort {$a<=>$b} keys %{$transcript_info->{$transcript}}) {
	    unless ($pos == -1) {push(@positions,$pos);}
	}
    } else {
	foreach my $pos (sort {$b<=>$a} keys %{$transcript_info->{$transcript}}) {
	    unless ($pos == -1) {push(@positions,$pos);}
	}
    }

    return @positions;   
    
}

sub get_transcripts {

    my ($gene,$version) = @_;

    my $GTDir = "/gscmnt/200/medseq/biodb/shared/misc/annotation/$version";

    my $gtdb = Bio::DB::Fasta->new($GTDir);
    my $transcripts = $gtdb->seq($gene, 1 => 1000);

    return unless $transcripts;
    return $transcripts;

}

sub get_codon {

    my ($transcript_info,$transcript,$protein_number,@positions) = @_;
    
    my ($p1,$p2,$p3,$b1,$b2,$b3);
    for my $pos (@positions) {
	
	my ($exon,$region) = split(/\,/,$transcript_info->{$transcript}->{$pos}->{exon});
	my $frame = $transcript_info->{$transcript}->{$pos}->{frame};
	my $aa_n = $transcript_info->{$transcript}->{$pos}->{aa_n};
	my $base = $transcript_info->{$transcript}->{$pos}->{base};
	my ($trans_pos) = $transcript_info->{$transcript}->{$pos}->{trans_pos};
	
	my $codon;

	if ($region eq "cds_exon" && $aa_n eq $protein_number) {
	    
	    if ($frame == 1) {
		$p1 = $pos;
		$b1 = $base;
	    } elsif ($frame == 2) {
		$p2 = $pos;
		$b2 = $base;
	    } elsif ($frame == 3) {
		$p3 = $pos;
		$b3 = $base;
		$codon = "$b1$b2$b3";
	    }

	    if ($codon) {

		return ($p1,$p2,$p3,$b1,$b2,$b3);

	    }
	}
    }
    return;
}

1;

