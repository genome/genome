package Genome::Model::Tools::Snp::ScreenNonsynonymousSnps;

use strict;
use warnings;
use Genome;
use Bio::Seq;
use Bio::SeqIO;


class Genome::Model::Tools::Snp::ScreenNonsynonymousSnps {
    is => 'Command',                    
    has => [ 
	
	list => {
	    is =>  'string',
	    doc   =>  "a file listing the nonsynonymous snps you want to run sift and polyphen on",
	},
	
	transcript => {
	    is  =>  'String',
	    doc   =>  "provide the transcript id",
	},
	organism => {
	    type  =>  'String',
	    doc   =>  "provide the organism either mouse or human; default is human",
	    is_optional  => 1,
	    default => 'human',
	},
	version => {
	    type  =>  'String',
	    doc   =>  "provide the imported annotation version; default for human is 54_36p_v2 and for mouse is 54_37g_v2",
	    is_optional  => 1,
	    default => '54_36p_v2',
	},
    nr_db => {
        type => 'String',
        doc => "use a different non-redundant protein db",
        is_optional => 1,
        default => "/gscmnt/sata847/info/genome-apipe/nr/nr",
    },
	]
	    
};

sub help_brief {
    return <<EOS
	This tool was design to retrieve derive the protien sequence for a transcript using the sequence from the datawarhouse and subsequently run sift and then polyphen on the list of nonsynonymous snps.
EOS
}

sub help_synopsis {
    return <<EOS


gmt snp screen-nonsynonymous-snps --transcript --list 


EOS
}


sub help_detail {
    return <<EOS 

gmt snp screen-nonsynonymous-snps --transcript --list 

EOS
}


sub execute {
    
    my $self = shift;
    
    my $transcript = $self->transcript;
    my $list = $self->list;

    my $organism = $self->organism;
    my $version = $self->version;
    if ($organism eq "mouse" && $version eq "54_36p_v2") { $version = "54_37g_v2"; }

    
########################### prep NS input list for sift and polyphen
    my $sortedlist;
    if (-f $list) {
	open(LIST,$list);
	while (<LIST>) {
	    chomp;
	    my $ns_snp = $_;
	    my ($sub1)=$ns_snp=~/^([A-Z])/;
	    my ($sub2)=$ns_snp=~/([A-Z])$/;
	    my ($sub3)=$ns_snp=~/(\d+)/;
	    $sortedlist->{$sub3}->{$ns_snp}->{sub1}=$sub1;
	    $sortedlist->{$sub3}->{$ns_snp}->{sub2}=$sub2;
	} close (LIST);
	open(SL,">SIFT_LIST");
	open(PL,">POLYPHEN_LIST");
	foreach my $sub3 (sort {$a<=>$b} keys %{$sortedlist}) {
	    foreach my $ns_snp (sort keys %{$sortedlist->{$sub3}}) {
		my $sub1 = $sortedlist->{$sub3}->{$ns_snp}->{sub1};
		my $sub2 = $sortedlist->{$sub3}->{$ns_snp}->{sub2};
		my $aa = "$sub1$sub3$sub2";
		print SL qq($aa\n);
		print PL "0\t0\t$transcript\t$sub3\t$sub1\t$sub2\n";
	    }
	}
	close SL;
	close PL;
    } else {
	print qq(Could not find the nonsynonymous snps list\n);
	return 0;
    }
########################### get the protein sequence
    
 
    my $eianame = "NCBI-" . $organism . ".ensembl";
    my $gianame = "NCBI-" . $organism . ".genbank";

    my $ensembl_build = Genome::Model::ImportedAnnotation->get(name => $eianame)->build_by_version($version);
    my ($ensembl_data_directory) = $ensembl_build->determine_data_directory;

    my $genbank_build = Genome::Model::ImportedAnnotation->get(name => $gianame)->build_by_version($version);
    my ($genbank_data_directory) = $genbank_build->determine_data_directory;

    my $t;
    if ($transcript =~/^ENS/){ #ENST for Human ENSMUST
	($t) = Genome::Transcript->get( transcript_name =>$transcript, data_directory => $ensembl_data_directory, reference_build_id => $ensembl_build->reference_sequence_id);
    }else{
	($t) = Genome::Transcript->get( transcript_name =>$transcript, data_directory => $genbank_data_directory, reference_build_id => $genbank_build->reference_sequence_id)
    }
    unless ($t) {print qq(\nCould not find a transcript object for $transcript from the $organism data warehouse\nWill exit the program now\n\n);exit(1);}
    

    my $tseq = $t->cds_full_nucleotide_sequence;
    
    unless ($tseq) {
	print qq(Could not find the nucleotide sequence for $transcript in the dw\n);
	return 0;
    }
    
    my $fasta = "$transcript.fasta";
    open(FA,">$fasta");
    print FA qq(>$transcript.fasta\n$tseq\n);
    close (FA);
    
    my $display_id = $transcript;
    
    my $newseq = Bio::Seq->new( -display_id => $display_id, -seq => $tseq );
    my $protein = $newseq->translate->seq;
    
    if ( $protein =~ /^[A-Z]+\*$/ ) {
	
	my $protein_file = "$transcript.protein.fasta";
	open(PFA,">$protein_file");
	print PFA qq(>$transcript $transcript\n$protein\n);
	close (PFA);
	
	my $db_location = $self->nr_db;

	if (-f "$transcript.protein.time.error") {system qq(rm $transcript.protein.time.error);}

	my $siftcmd = "SIFT.csh $protein_file $db_location SIFT_LIST";
	print "$siftcmd \nRunning...\n\n";
	
	`$siftcmd > ./SIFT_$protein_file.out`;

	if (-f "$transcript.protein.SIFTprediction") {
	    print qq(Sift results can be viewed in $transcript.protein.SIFTprediction\n);
	} else {
	    print qq(SIFT failed to produce a result\n);
	}

	my $prediction_file = $protein_file;
	$prediction_file =~ s/.fasta$/.PPHprediction/;
	my $pphdire=qq(/gscmnt/200/medseq/analysis/software/scripts/pph.sh);
	my $polycmd = "$pphdire -s $protein_file POLYPHEN_LIST";
	
	print "$polycmd \nRunning...\n\n";
	`$polycmd > ./$prediction_file`;
	if (-f "$transcript.protein.PPHprediction") {
	    print qq(polyphen results can be viewed in $transcript.protein.PPHprediction\n);
	} else {
	    print qq(polyphen failed to produce a result\n);
	}

	#clean up
	system qq(rm SIFT_LIST POLYPHEN_LIST);

    } else {
	
	print qq(protein sequence does not conform to requirements and sift will not run on it.\n);
	
    }
}


1;
