package Genome::Model::Tools::Consed::AceReference;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Consed::AceReference {
    is => 'Command',                       
    has => [ 
	ace_file => {
            type  =>  'String',
            doc  => "manditory; ace file name",
	},
	name_and_number => {
            type  =>  'Boolean',
            doc  => "gets the name and number of the contig with your reference sequence in it; reference defined by a read ending in .c1",
  	    is_optional  => 1,
	},
	no_stdout => {
	    is => 'Boolean',
	    doc   =>  "Use this option if you do not want the info to print to stdout. Default is to print to stdout.",
	    is_optional  => 1,
	},
	population_file => {
            type  =>  'String',
	    doc   =>  "provide a name for your pop.txt file",
	    is_optional  => 1,
	},
	 pretty_source_1 => {
            type  =>  'String',
            doc  => 'pretty-source-1 to set the sample start pos',
	    is_optional  => 1,
	    default => 1,
        },
        pretty_source_2 => {
            type  =>  'String',
            doc  => 'pretty-source-2 to set the sample stop pos default is 10 this value defines the end of the sample name in the population file',
	    is_optional  => 1,
	    default => 20,
        },

	], 
};

sub help_brief {                            

"A tool to get info about your reference contig"

}

sub help_synopsis { 

    return <<EOS

	gmt consed ace-reference -h

EOS
}
sub help_detail {
    return 'This tool was designed to get information about ace files that were assembled by aligning reads under a reference sequence "a fake trace ending .c1 derived from a slice out of a reference geneome"';
}


sub execute {

    my $self = shift;

    my $ace_reference;

    if ($self->name_and_number || $self->population_file) {
	($ace_reference) = &name_and_number($self,$ace_reference);
	my $Contig_number = $ace_reference->{Contig_number};
	my $reseqid = $ace_reference->{reseqid};
	unless ($self->no_stdout) {
	    print qq(Contig_number => $Contig_number, Refseq_id => $reseqid\n);
	}
    }

    return unless $ace_reference;
    return $ace_reference;

}


sub name_and_number {

    use Genome::Model::Tools::Pcap::Ace;

    my ($self,$ace_reference) = @_;
    
    my $ace_file = $self->ace_file;
    unless (-f $ace_file) {$self->error_message("could see the ace file");return;}
    my $ao = new Genome::Model::Tools::Pcap::Ace(input_file => $ace_file);
    
    my $ace_ref;
    my @number;
    my @name;


    my ($pretty_source_1,$pretty_source_2) = &source($self);
    foreach my $contig_number (@{ $ao->get_contig_names }) {
	my $contig = $ao->get_contig($contig_number);
	
	if (grep { /\.c1$/ } keys %{ $contig->reads }) {
	    push(@number,$contig_number);
	}
	
	foreach my $read_name (keys %{ $contig->reads }) {
	    if ($read_name =~ /(\S+\.c1)$/) {
		push(@name,$read_name);
	    } elsif ($self->population_file) {
		my $sample = substr($read_name,$pretty_source_1,$pretty_source_2);
		$ace_reference->{sample}->{$sample}=1;
	    }
	}
    }

    unless (@name && @number) {$self->error_message("couldn't find the contig name and number");return;}

    if (@name > 1) {$self->error_message("there is more than one reference sequence in your assembly");}
    if (@number > 1) {$self->error_message("there is more than one contig with a reference sequence in your assembly");}
    
    my $name = join 'name' , @name;
    my $number = join 'number' , @number;
    
    $ace_reference->{Contig_number}=$number;
    $ace_reference->{reseqid}=$name;

    if ($self->population_file) {
	my $pop = &write_population_file($self,$ace_reference);
	unless ($pop) {$self->error_message("failed to write the population file\n");}
    }
    
    return ($ace_reference);
    
}

sub write_population_file {

    my ($self,$ace_reference) = @_;
    my $population_file = $self->population_file;
    open (POP,">$population_file") || $self->error_message("couldn't open $population_file for writing") && return;
    my $reseqid = $ace_reference->{reseqid};
    print POP qq($reseqid\n);
    foreach my $sample (sort keys %{$ace_reference->{sample}}) {
	print POP qq($sample\n);
    }
    close POP;
    return unless (-f $population_file);
    return $population_file;
}

sub source {

    my ($self) = @_;
    
    my $pretty_source_1 = $self->pretty_source_1;
    my $pretty_source_2 = $self->pretty_source_2;

    if ($pretty_source_1) {
	my $d = $pretty_source_1;
	$pretty_source_1 = ($d - 1);
    } else {
	$pretty_source_1 = 0;
    }
    if ($pretty_source_2) {
	if ($pretty_source_1 != 0) {
	    my $d = $pretty_source_2;
	    $pretty_source_2 = ( $d - $pretty_source_1 );
	}
    } else {
	$pretty_source_2 = 10;
	if ($pretty_source_1 != 1) {
	    my $d = $pretty_source_2;
	    $pretty_source_2 = $d - $pretty_source_1;
	}
    }
    return ($pretty_source_1,$pretty_source_2);   
}

1;
