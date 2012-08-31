package Genome::Model::Tools::Annotate::GeneRegions;

use strict;
use warnings;
use Genome;

use Bio::DB::Fasta;

class Genome::Model::Tools::Annotate::GeneRegions {
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
	   output => {
	    type  =>  'String',
	    doc   =>  "provide a file name to write you transcript information to. Default is to print to stdout.",
	    is_optional  => 1,
	},
	list => {
	    type  =>  'String',
	    doc   =>  "input list tab/space delimited with chromosome ie {1,2,...,22,X,Y,M}, Start coordinate, Stop coordinate. Or a hugo gene name",
	},
	mode => {
	    type  =>  'String',
	    doc   =>  "enter the source for your search",
	    valid_values => ['coordinates','gene'],
	},
    ], 
};


sub help_synopsis {
    return <<EOS

gmt annotate gene-region -h

EOS
}

sub help_detail {
    return <<EOS 

	This tool will find all genes and transcripts along with region types that coincides with the coordinates or gene names in your list.

	gmt annotate gene-region -mode coordinates -list some_file
	
	gmt annotate gene-region -mode gene -list some_file
	gmt annotate gene-region -mode gene -list EGFR
	gmt annotate gene-region -mode gene -list EGFR,KRAS

EOS
}

sub execute {

    my $self = shift;
    my $output = $self->output;
    if ($output) {
	open(OUT,">$output") || $self->error_message("\n\nCouldn't open the file $output for writting\n\n.") && return;
    }

    my $gene_list = &get_gene_list($self);
    unless ($gene_list) {$self->error_message( "\n\nCouldn't generate a gene list\n\n.");return;}

    my $organism = $self->organism;
    my $version = $self->version;

    if ($organism eq "mouse") {
	if ($version eq "54_36p_v2") {$version = "54_37g_v2";}
    }
    my $results = 0;

    foreach my $hugo_name (sort keys %{$gene_list}) {
	foreach my $line (sort keys %{$gene_list->{$hugo_name}}) {
	    unless ($hugo_name eq $line) {
		if ($output) {
		    print OUT qq($line\n);
		} else {
		    print qq($line\n);
		}
	    }
	}

	my $transcripts = &get_transcripts($hugo_name);
	unless ($transcripts) {print qq(no transcripts found for $hugo_name\n);next;}
	my @transcripts = split(/\,/,$transcripts);
	
	for my $transcript (@transcripts) {
	    
	    my ($info) = Genome::Model::Tools::Annotate::TranscriptSequence->execute(transcript => $transcript, no_stdout => '1',version => $version, organism => $organism);
	    unless ($info) {print qq(no info found for $hugo_name $transcript\n); next;}
	    my $transcript_info = $info->{result};
	    unless ($transcript_info) {print qq(no transcript info found for $hugo_name $transcript\n); next;}
	    
	    my $source_line = $transcript_info->{$transcript}->{'-1'}->{source_line};
	    my $strand = $transcript_info->{$transcript}->{'-1'}->{strand};
	    my $transcript_start = $transcript_info->{$transcript}->{'-1'}->{first_base};
	    my $transcript_stop = $transcript_info->{$transcript}->{'-1'}->{last_base};
	    my $g_id = $transcript_info->{$transcript}->{'-1'}->{gene_id};
	    my $chromosome = $transcript_info->{$transcript}->{'-1'}->{chromosome};
	    my $hugo_gene_name = $transcript_info->{$transcript}->{'-1'}->{hugo_gene_name};
	    
	    my $cds;
	    foreach my $pos (sort {$a<=>$b} keys %{$transcript_info->{$transcript}}) {
		unless ($pos == -1) {
		    my $codon = $transcript_info->{$transcript}->{$pos}->{aa_n};
		    my $base = $transcript_info->{$transcript}->{$pos}->{base};
		    my $exon = $transcript_info->{$transcript}->{$pos}->{exon};
		    my $frame = $transcript_info->{$transcript}->{$pos}->{frame};
		    my $range = $transcript_info->{$transcript}->{$pos}->{range};
		    my ($exon_number,$exon_type) = split(/\,/,$exon);
		    $cds->{$exon_number}->{$exon_type}=$range;
		}
	    }
	    
	    my $threeprimeendutr; #5'to3'
	    foreach my $exon_number (sort {$a<=>$b} keys %{$cds}) {
		my $utr_range = $cds->{$exon_number}->{utr_exon};
		my $cds_range = $cds->{$exon_number}->{cds_exon};
		
		my $fiveprimeendutr;
		if ($utr_range) {
		    unless ($threeprimeendutr) {
			$utr_range =~ s/\-/\t/;
			$fiveprimeendutr = 1;
			$results =1;
			if ($output) {
			    print OUT qq($hugo_name\t$strand\t$chromosome\t$transcript\t$exon_number\tutr_exon\t$utr_range\n);
			} else {
			    print qq($hugo_name\t$strand\t$chromosome\t$transcript\t$exon_number\tutr_exon\t$utr_range\n);
			} 
		    }
		}
		if ($cds_range) {
		    $threeprimeendutr = 1;
		    $cds_range =~ s/\-/\t/;
		    $results =1;
		    if ($output) {
			print OUT qq($hugo_name\t$strand\t$chromosome\t$transcript\t$exon_number\tcds_exon\t$cds_range\n);
		    } else {
			print qq($hugo_name\t$strand\t$chromosome\t$transcript\t$exon_number\tcds_exon\t$cds_range\n);
		    } 
		}
		if ($utr_range && $threeprimeendutr) {
		    unless ($fiveprimeendutr) {
			$utr_range =~ s/\-/\t/;
			$results =1;
			if ($output) {
			    print OUT qq($hugo_name\t$strand\t$chromosome\t$transcript\t$exon_number\tutr_exon\t$utr_range\n);
			} else {
			    print qq($hugo_name\t$strand\t$chromosome\t$transcript\t$exon_number\tutr_exon\t$utr_range\n);
			}
		    }
		}
	    }
	}
    }
    return $results;
}


sub get_transcripts {

    my ($gene) = @_;

    # bdericks: This might need to be changed to new 54_36p_v2 version
    my $GTDir = "/gscmnt/200/medseq/biodb/shared/misc/annotation/54_36p";

    my $gtdb = Bio::DB::Fasta->new($GTDir);
    my $transcripts = $gtdb->seq($gene, 1 => 1000);

    return unless $transcripts;
    return $transcripts;

}

sub get_gene_list {

    my ($self) = @_;
    my $mode = $self->mode;

    my $gene_list;
    if ($mode eq "coordinates") {
	$gene_list = &parse_coordinates($self);
    } elsif ($mode eq "gene") {
	$gene_list = &parse_genes($self);
    }
    return unless $gene_list;
    return ($gene_list);
}

sub parse_genes {
    
    my $gene_list;
    my ($self) = @_;
    my $list = $self->list;
    
    if (-f $list) {
	open (LIST,"$list") || $self->error_message( "\n\nCouldn't open your list $list\n\n.") && return;
	while (<LIST>) {
	    chomp;
	    my $gene = $_;
	    $gene_list->{$gene}->{$gene}=1;
	} 
	close LIST;
	
    } else {
	my @genes = split(/\,/,$list);
	for my $gene (@genes) {
	    $gene_list->{$gene}->{$gene}=1;
	}
    }
    
    return unless ($gene_list);
    return ($gene_list);

}

sub parse_coordinates {
    
    my ($self) = @_;
    my $organism = $self->organism;
    my $version = $self->version;

    if ($organism eq "mouse") {
	if ($version eq "54_36p_v2") {$version = "54_37g_v2";}
    }
    my $list = $self->list;
    my $anno_db = "NCBI-$organism.combined-annotation";

    unless (-f $list) {
	system qq(gmt annotate gene-regions --help);
	$self->error_message("\nYour file $list was not found.\t\n\tPlease check that your list is in place and try again.\n\n\n");
	return;
    }
    
    my @genes;
    open(LIST,"$list") || $self->error_message( "\n\nCouldn't open your list $list\n\n.") && return;
    my $gene_list;

    my $coords;
    while (<LIST>) {
	chomp;
	my $line = $_;
	my ($chr,$start,$stop) = (split(/[\s]+/,$line))[0,1,2];

	unless ($chr =~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|M)$/ && $start =~ /^[\d]+$/ && $stop =~ /^[\d]+$/) {
	    $self->error_message("$line does not represent valid coordinates\n");
	    next;
	}
	$coords->{$chr}->{$start}->{$stop}=$line;
    }
    close LIST;

    foreach my $chr (sort keys %{$coords}) {
	my $ti = Genome::Model->get(name => $anno_db)->build_by_version($version)->transcript_iterator(chrom_name => $chr);
	unless ($ti) {$self->error_message("\ncouldn't get a transcript iterator for $chr\n\n");next;}
	my $transcript_window =  Genome::Utility::Window::Transcript->create(iterator => $ti);
	unless ($transcript_window) {$self->error_message("\ncouldn't get a transcript window for $chr\n\n");next;}
	
	foreach my $start (sort {$a<=>$b} keys %{$coords->{$chr}}) {
	    foreach my $stop (sort {$a<=>$b} keys %{$coords->{$chr}->{$start}}) {
		my $line = $coords->{$chr}->{$start}->{$stop};
		for my $t ($transcript_window->scroll($start,$stop)) {
		    my $gene = $t->gene;
		    my $gene_name;
		    if ($gene) {
			my @g_id = $t->gene_id;
			my ($gene_id) = $g_id[0];
			
			$gene_name = $gene->name || "unknown\-$gene_id";
		    }
		    
		    $gene_list->{$gene_name}->{$line}=1;
		}
	    }
	}
	
	undef($ti);
	undef($transcript_window);
	
    }
    
    return unless ($gene_list);
    return ($gene_list);
    
}

1;
