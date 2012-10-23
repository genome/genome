package Genome::Model::Tools::WuBlast::OligoHit;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::WuBlast::OligoHit {
    is => 'Command',                    
    has => [                                # specify the command's properties (parameters) <--- 
        blastdb     => { type => 'String',      doc => "give the full path to your blast database or use the default human build 36", is_optional => 1 },
        fasta     => { type => 'String',      doc => "give your fasta file to be blasted" },
	blast_options     => { type => 'String'      ,doc => "give all or none blast options as a quoted string", is_optional => 1  },
	blast_results     => { type => 'String'      ,doc => "result of the blasted fasta not intended as an input", is_optional => 1  },
    ], 
};

sub help_brief {                            # keep this to just a few words <---
    "provide a fasta file and get human B36 blast result back"                 
}

sub help_synopsis {                         # replace the text below with real examples <---
    return <<EOS
gmt wu-blast oligo-hit --fasta subject_fasta
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

	please provide the subject_fasta file in standard fasta format 
	for example gmt wu-blast oligo-hit --fasta /gscmnt/sata147/info/medseq/rmeyer/scratch/21_46356247_46356978.c1.refseq.fasta
	should result in
chr21   21_46356247_46356978.c1.refseq.fasta    0.0     100.0   100.0   1464.4  0.      46944323        732     732     46356247-46356978:1-732(100);
	where the line breaks down like this
subject_id  query_id  subject_cov  query_cov  percent_identity  bit_score  p-value  subject_length  query_length  alignment_bases  HSPs(subject:query)

EOS
}

sub execute {                               # replace with real execution logic.
    my $self = shift;

    my $blastdb = $self->blastdb;
    unless ($blastdb) { $blastdb = "/gscmnt/200/medseq/analysis/software/resources/B36/HS36.fa"; }

    my $fasta = $self->fasta;
    unless ($fasta) { sub help_detail; }

    my $blast_options = $self->blast_options;
    unless ($blast_options) {
	my $base_count=0;
	my $ref_bases;
	open(FASTA,"$fasta");
	while (<FASTA>) {
	    chomp;
	    my $seq = $_;
	    unless ($seq =~ /^>/) {
		my @sequence = split(//, $seq);
		foreach my $base (@sequence){
		    $base =~ s/($base)/\U$1/; 
		    $base_count++;
		    $ref_bases->{$base_count}->{base}=$base;
		}
	    }
	}
	
	my $n = $base_count + 1;
	my $s = $base_count;
	my $s2 = $base_count;
	
	$blast_options = "-nogap -M=1 -N=-$n -S2=$s2 -S=$s topcomboN=1";
	
    }

    print "Running the blast command:\n";
	

    my $result = `blastn /gscmnt/200/medseq/analysis/software/resources/B36/HS36.fa $fasta $blast_options | blast2gll -s`;
    
    $self->blast_results($result);
    print qq($result\n);
    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

1;

