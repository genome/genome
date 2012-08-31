package Genome::Model::Tools::RefSeq::Fasta;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::RefSeq::Fasta {
    is => 'Command',                       
    has => [ 
	refseq_fasta => {
            type  =>  'String',
            doc  => "refseq fasta file",
	},
	no_stdout => {
	    is => 'Boolean',
	    doc   =>  "Use this option if you do not want the info to print to stdout. Default is to print to stdout.",
	    is_optional  => 1,
	},

	], 
};

sub help_brief {                            

"A tool to parse your refseq fasta file header"

}

sub help_synopsis { 

    return <<EOS

	gmt consed ace-reference -h

EOS
}
sub help_detail {
    return 'This tool was designed to parse the header of a refseq.fasta file header that was use to create the fake trace in a consed ace file. The header should have a format similar to this ">chr3:103057567:103057936.refseq.fasta Chr:3, Coords 103057567-103057936, Ori (+), comment"
it will return the fasta name, chromosome, genomic start and stop coordinates, genomic_coord used to convert to reference coordinate, the length, and the orientation ';
}


sub execute {

    my $self = shift;

    my $refseq_fasta = $self->refseq_fasta;
    unless (-f $refseq_fasta) {$self->error_message("could see the refseq fasta file");return;}
    my $refseq_header = &parse_header($self,$refseq_fasta);

    my $orientation = $refseq_header->{orientation};
    my $genomic_coord = $refseq_header->{genomic_coord};
    my $chromosome = $refseq_header->{chromosome};
    my $length = $refseq_header->{length};
    my $start = $refseq_header->{start};
    my $stop = $refseq_header->{stop};
    my $name = $refseq_header->{name};
    unless ($self->no_stdout) {
	print qq(name => $name, chromosome => $chromosome, orientation => $orientation, genomic_coord => $genomic_coord, start => $start, stop => $stop, length => $length\n);
    }
    return unless $refseq_header;
    return $refseq_header;

}

sub parse_header {
    
    my $refseq_header;

    my ($self,$refseq_fasta) = @_;

    open(REF,$refseq_fasta) || $self->error_message("couldn't open the refseq fasta file") && return;
    while (<REF>) {
	chomp;
	my $line=$_;
	next unless ($line =~ /\>/);
	
	my ($name) = $line =~ /\>([\S]+)/;
	my ($orientation,$genomic_coord,$ref_seq_length);
	
	my ($fisrt_coord,$second_coord) = $line =~ /Amplicon\_Coords\:[\s]+(\d+)\S(\d+)/; ## Throw back from irregularly formated headers 
	unless ($fisrt_coord && $second_coord) {
	    ($fisrt_coord,$second_coord) = $line =~ /Coords[\s]+(\d+)\S(\d+)/;
	}
	my ($chromosome) = $line =~ /Chr\:([\S]+)/;
	$chromosome =~ s/\,//;

	my ($ncbi_build_no) = $line =~ /NCBI Build[\s]+(\d+)/;
	my ($gene) = $line =~ /GeneName:(\S+),/;
	my ($gene_id) = $line =~ /GeneID:(\S+),/;

	unless ($gene) { $gene = "-"; }
	unless ($gene_id) { ($gene_id) = $line =~ /Target (\S+)\,/; }
	unless ($gene_id) { $gene_id = "-"; }

	if ($line=~ /Ori[\s]+\(\+\)/) {
	    $orientation="plus";
	    $genomic_coord = $fisrt_coord - 1;
	    #} elsif ($ori eq "-") {
	} elsif ($line=~ /Ori[\s]+\(\-\)/) {
	    $orientation="minus";
	    $genomic_coord = $second_coord + 1;
	}
	if ($fisrt_coord > $second_coord) { #WARNING:
	    $ref_seq_length = ($fisrt_coord - $second_coord + 1);
	} else {
	    $ref_seq_length = ($second_coord - $fisrt_coord + 1);
	}
	
	unless ($orientation && $genomic_coord && $chromosome && $ref_seq_length && $fisrt_coord && $second_coord && $name) {
	    close REF;
	    $self->error_message("couldn't parse the information from the refseq fasta file header");
	    return;
	}

	unless ($ncbi_build_no) {$ncbi_build_no = "-";}
	$refseq_header->{orientation}=$orientation;
	$refseq_header->{genomic_coord}=$genomic_coord;
	$refseq_header->{chromosome}=$chromosome;
	$refseq_header->{length}=$ref_seq_length;
	$refseq_header->{start}=$fisrt_coord;
	$refseq_header->{stop}=$second_coord;
	$refseq_header->{name}=$name;
	$refseq_header->{ncbi_build_no}=$ncbi_build_no;
	$refseq_header->{gene}=$gene;
	$refseq_header->{gene_id}=$gene_id;
	$refseq_header->{header}=$line;

	if ($refseq_header) {
	    close REF;
	    return $refseq_header;
	}
    }
    close REF;
    
    return unless $refseq_header;
    return $refseq_header;
    
}


1;
