package Genome::Model::Tools::Annotate::Translocation;

use strict;
use warnings;
use Genome;
use IPC::Run;

use Bio::Seq;
use Bio::Tools::CodonTable;
use Bio::DB::Fasta;

class Genome::Model::Tools::Annotate::Translocation {
    is => 'Command',                       
    has => [ 
	breakdancer => {
	    type  =>  'String',
	    doc   =>  "provide breakdancer output",
	    is_optional  => 1,
	},
	fasta => {
	    type  =>  'String',
	    doc   =>  "provide a fasta contig for the translocation",
	    is_optional  => 1,
	},
	organism => {
	    is => 'String',
	    doc   =>  "use this option if you're working with mouse",
	    is_optional  => 1,
	    default => 'human',
	},
	version => {
	    type  =>  'String',
	    doc   =>  "provide the imported annotation version; default for human is 54_36p_v2 and for mouse is 54_37g_v2",
	    is_optional  => 1,
	},

    ], 
};
     

sub help_brief {   
#return <<EOS
	"gmt annotation translocation"
#EOS  
}

sub help_synopsis {
    return <<EOS

	gmt annotation translocation -breakdancer

	or

	gmt annotation translocation -fasta
EOS
}

sub help_detail { 
    return <<EOS 

	you can use the -fasta option if there is anassembly fasta for your translocation
	otherwise you can use the breakdancer output with the -breakdancer option
	organsim mouse or human
EOS
}


sub execute {

    my $self = shift;
    my $fasta = $self->fasta;
    my $breakdancer = $self->breakdancer;
    my $organism = $self->organism;

    unless ($fasta && -e "$fasta" || $breakdancer) { 
	print "Could not find an input file.\n";
	exit(1);
    }
    
    
    #my ($np1,$np2,$np1_orientation,$np2_orientation,$name);
    my ($bp1_chr,$bp1_pos,$bp2_chr,$bp2_pos,$np1_orientation,$np2_orientation,$name,$np1,$np2);
    if ($breakdancer) {
	($bp1_chr,$bp1_pos,$bp2_chr,$bp2_pos,$np1_orientation,$np2_orientation,$name) = &parse_breakdancer($breakdancer);
	$np1 = $bp1_pos;
	$np2 = $bp2_pos;
    }
    
    if ($fasta) {
	unless ($bp1_chr && $bp1_pos && $bp2_chr && $bp2_pos) {
	    ($bp1_chr,$bp1_pos,$bp2_chr,$bp2_pos) = $fasta =~ /(chr[\S]+)\_(\d+)\-(chr[\S]+)\_(\d+)\-CTX/;
	}
	
	#$name = "$bp1_chr\_$bp1_pos\-$bp2_chr\_$bp2_pos";
	($name) = $fasta =~ /(\S+)\.fasta/;
	my $blast_out = "$name.parsed.blast.out";
	my $blast_db;
	if ($organism eq "human") {
	    $blast_db = "/gscmnt/200/medseq/analysis/software/resources/B36/HS36.fa";
	} else {
	    $blast_db = "/gscmnt/200/medseq/analysis/software/resources/MouseB37/MS37.fa";
	}

	#system qq(blastn /gscmnt/200/medseq/analysis/software/resources/B36/HS36.fa $fasta | blast2gll -s > $blast_out);
	system qq(blastn $blast_db $fasta | blast2gll -s > $blast_out);

	my $alignment_locations = &parse_blast_out($blast_out,$bp1_chr,$bp1_pos,$bp2_chr,$bp2_pos);
	
	my ($np1_alignment,$np2_alignment);
	($np1,$np1_alignment) = split(/\,/,$alignment_locations->{$bp1_chr}->{$bp1_pos});
	unless($np1) {print qq(Failed to identify a location on $bp2_chr for $bp2_pos\nProcessing will continue using $bp1_pos\n); $np1 = $bp1_pos;}
	($np2,$np2_alignment) = split(/\,/,$alignment_locations->{$bp2_chr}->{$bp2_pos});
	unless($np2) {print qq(Failed to identify a location on $bp2_chr for $bp2_pos\nProcessing will continue using $bp2_pos\n); $np2 = $bp2_pos;}
	
	my ($np1gs,$np1ge,$np1rs,$np1re) = $np1_alignment =~ /(\d+)\-(\d+)\:(\d+)\-(\d+)/;
	my ($np2gs,$np2ge,$np2rs,$np2re) = $np2_alignment =~ /(\d+)\-(\d+)\:(\d+)\-(\d+)/;
	
	my ($np1go,$np1ro,$np2go,$np2ro);
	if ($np1gs && $np1ge) {
	    $np1go = &get_alignment_orientation($np1gs,$np1ge);
	    $np1ro = &get_alignment_orientation($np1rs,$np1re);
	    if ($np1go) {
		$np1_orientation = $np1go;
	    }
	}
	if ($np2gs && $np2ge) {
	    $np2go = &get_alignment_orientation($np2gs,$np2ge);
	    $np2ro = &get_alignment_orientation($np2rs,$np2re);
	    if ($np2go) {
		$np2_orientation = $np2go;
	    } 
	} 
    }
    
    my $np1_rb = &get_ref_base($np1,$np1,$bp1_chr,$organism);
    my ($rev_np1_rb) = &reverse_complement_allele ($np1_rb);
    
    my $np2_rb = &get_ref_base($np2,$np2,$bp2_chr,$organism);
    my ($rev_np2_rb) = &reverse_complement_allele ($np2_rb);
    
    my $nbp1_chr = $bp1_chr;
    my $nbp2_chr = $bp2_chr;
    $nbp1_chr =~ s/chr//;
    $nbp2_chr =~ s/chr//;
    $nbp1_chr =~ s/x/X/;
    $nbp2_chr =~ s/y/Y/;
    $bp1_chr =~ s/x/X/;

    $bp2_chr =~ s/y/Y/;
    my $variant_file = "$name.annotation.list";
    open(ANNOLIST,">$variant_file");
    print ANNOLIST qq($nbp1_chr\t$np1\t$np1\t$np1_rb\t$rev_np1_rb\n);
    print ANNOLIST qq($nbp2_chr\t$np2\t$np2\t$np2_rb\t$rev_np2_rb\n);
    close(ANNOLIST);


    
    my $version = $self->version;
    unless ($version) { if ($organism eq "mouse") { $version = "54_37g_v2" ; } else { $version = "54_36p_v2" ; } }
    my $references_transcripts = "NCBI-$organism.combined-annotation/$version";

    my $annotation_file = "$name.annotated.list.csv";
    
    my $AnnotateTranscriptVariants = Genome::Model::Tools::Annotate::TranscriptVariants->create(variant_file=>$variant_file,output_file=>$annotation_file,flank_range=>"0",extra_details => "1", reference_transcripts =>$references_transcripts);

    $AnnotateTranscriptVariants->execute() || die 'did not get annotation for the annotation_file $annotation_file';

    my ($annotation_info) = &parse_annotation($annotation_file);
    
    my $bp1_transcript = $annotation_info->{$bp1_chr}->{$np1}->{transcript};
    my $bp2_transcript = $annotation_info->{$bp2_chr}->{$np2}->{transcript};
    
    unless ($bp1_transcript && $bp2_transcript) {
	if ($bp1_transcript) {
	    print qq(\nNo transcript was indentified by the annotator on $bp2_chr at $np2.\nWill exit the program now\n\n);
	    exit 1;
	} elsif ($bp2_transcript) {
	    print qq(\nNo transcript was indentified by the annotator on $bp1_chr at $np1.\nWill exit the program now\n\n);
	    exit 1;
	} else {
	    print qq(\nNo transcript was indentified by the annotator on $bp1_chr at $np1 nor on $bp2_chr at $np2.\nWill exit the program now\n\n);
	    exit 1;
	}
    }
    

    my ($AnnotateTranscriptSequence) = Genome::Model::Tools::Annotate::TranscriptSequence->create(transcript => "$bp1_transcript,$bp2_transcript", trans_pos => "$np1,$np2", organism => "$organism", version => "$version");
    my ($transcript_info) = $AnnotateTranscriptSequence->execute() || die 'some error message goes here';

    my $bp1_region = $transcript_info->{$bp1_transcript}->{-1}->{region};
    my $bp2_region = $transcript_info->{$bp2_transcript}->{-1}->{region};

    
    my $bp1_gene = $annotation_info->{$bp1_chr}->{$np1}->{gene};
    my $bp1_strand = $annotation_info->{$bp1_chr}->{$np1}->{strand};
    my ($bp1_flip) = &compare_orientations($bp1_strand,$np1_orientation);
    if ($bp1_flip eq "YES") {
	#print qq($bp1_transcript used for $bp1_chr $np1 will be flipped as a result of this transloction and will  not start with a Methionine and will likely hit a premature stop codon\n);
    }
    
    my $bp2_gene = $annotation_info->{$bp2_chr}->{$np2}->{gene};
    my $bp2_strand = $annotation_info->{$bp2_chr}->{$np2}->{strand};
    my ($bp2_flip) = &compare_orientations($bp1_strand,$np2_orientation);
    
        
#####  Does the orientation make sense. 
    
#1 # if bp1 eq 5' to 3' and bp2 eq 5' to 3' this has one product expected fairly straight forward  ==> 5' bp1 to 3' bp2 ($bp1_flip eq NO && $bp2_flip NO)
#2 # if bp1 eq 3' to 5' and bp2 eq 5' to 3' no promoter no product  ==> 3' bp1 to 3' bp2 ($bp1_flip eq YES && $bp2_flip NO)
#3 # if bp1 eq 3' to 5' and bp2 eq 3' to 5' this might be considered a double flip were we could reverse the transcript order  ==> 5' bp2 to 3' bp1 ($bp1_flip eq YES && $bp2_flip YES)
#4 # if bp1 eq 5' to 3' and bp2 eq 3' to 5' this may result in two possible products
    #a ==> 5' bp1 to 3' bp2
    #b ==> 5' bp2 to 3' bp1
    
### is the translocation position is in an intron or in an exon
###    if both coordinates fall in an intron then all applicable exons get used
###    if one is in an intron and the other is in an exon then all the sequence from the exon that the other is in will be omitted from the protien
###    if both coordinates fall in an exon then they should become contiguous
###            5' UTR can be used if that is the exon the second coordinate falls in when both coordinates fall in an exon
###            3' UTR can be used if a stop codon hasn't been identifed in the second transcript before reaching it
    
#####

    #my $translocation;
    #my $translocation_boundary;
    
    open(OUT,">$name.translocation.info.txt");  ### Results written to this file
    my ($translocation,$count);
    if ($bp1_flip eq "NO" && $bp2_flip eq "NO") {
	
	my $base_count=0;
	my $side = "BP1";
	($translocation,$count,$bp1_region) = &get_translocation($translocation,$transcript_info,$bp1_gene,$bp1_strand,$bp1_transcript,$base_count,$side,$bp1_flip);
	$side = "BP2";
	#$translocation_boundary=$count;
	($translocation,$count,$bp2_region) = &get_translocation($translocation,$transcript_info,$bp2_gene,$bp2_strand,$bp2_transcript,$count,$side,$bp2_flip);
	my $regionside = "BP2";
	&get_translocation_info($transcript_info,$translocation,$bp1_region,$bp2_region,$regionside,$organism,$bp1_gene,$bp1_strand,$bp1_transcript,$bp1_flip,$bp1_region,$bp2_gene,$bp2_strand,$bp2_transcript,$bp2_flip,$bp2_region,$np1,$np2,$np1_orientation,$np2_orientation,$name,$bp1_chr,$bp2_chr);
	
    }
    
    if ($bp1_flip eq "YES" && $bp2_flip eq "NO") {
	print qq(There is no promoter for this translocation, therefore no protein will be predicted\nWill exit the program now.\n);
	exit 1;
    }
    
    if ($bp1_flip eq "YES" && $bp2_flip eq "YES") {  ##This needs to be tested
	my $base_count=0;
	my $side = "BP1";
	($translocation,$count,$bp2_region) = &get_translocation($translocation,$transcript_info,$bp2_gene,$bp2_strand,$bp2_transcript,$base_count,$side,$bp2_flip);
	$side = "BP2";
	#$translocation_boundary=$count;
	($translocation,$count,$bp1_region) = &get_translocation($translocation,$transcript_info,$bp1_gene,$bp1_strand,$bp1_transcript,$count,$side,$bp1_flip);
	my $regionside = "BP1";
	&get_translocation_info($transcript_info,$translocation,$bp2_region,$bp1_region,$regionside,$organism,$bp1_gene,$bp1_strand,$bp1_transcript,$bp1_flip,$bp1_region,$bp2_gene,$bp2_strand,$bp2_transcript,$bp2_flip,$bp2_region,$np1,$np2,$np1_orientation,$np2_orientation,$name,$bp1_chr,$bp2_chr);
    }
    
    if ($bp1_flip eq "NO" && $bp2_flip eq "YES") {  ##This needs to be tested
	my $base_count=0;
	my $side = "BP1";
	($translocation,$count,$bp1_region) = &get_translocation($translocation,$transcript_info,$bp1_gene,$bp1_strand,$bp1_transcript,$base_count,$side,$bp1_flip);
	$side = "BP2";
	#$translocation_boundary=$count;
	($translocation,$count,$bp2_region) = &get_translocation($translocation,$transcript_info,$bp2_gene,$bp2_strand,$bp2_transcript,$count,$side,$bp2_flip);
	my $regionside = "BP2";
	&get_translocation_info($transcript_info,$translocation,$bp1_region,$bp2_region,$regionside,$organism,$bp1_gene,$bp1_strand,$bp1_transcript,$bp1_flip,$bp1_region,$bp2_gene,$bp2_strand,$bp2_transcript,$bp2_flip,$bp2_region,$np1,$np2,$np1_orientation,$np2_orientation,$name,$bp1_chr,$bp2_chr);
	
###AND
	
	undef($translocation);
	#undef($translocation_boundary);
	
	$base_count=0;
	$count=0;
	$side = "BP1";
	($translocation,$count,$bp2_region) = &get_translocation($translocation,$transcript_info,$bp2_gene,$bp2_strand,$bp2_transcript,$base_count,$side,$bp2_flip);
	$side = "BP2";
	#$translocation_boundary=$count;
	($translocation,$count,$bp1_region) = &get_translocation($translocation,$transcript_info,$bp1_gene,$bp1_strand,$bp1_transcript,$count,$side,$bp1_flip);
	$regionside = "BP1";
	&get_translocation_info($transcript_info,$translocation,$bp2_region,$bp1_region,$regionside,$organism,$bp1_gene,$bp1_strand,$bp1_transcript,$bp1_flip,$bp1_region,$bp2_gene,$bp2_strand,$bp2_transcript,$bp2_flip,$bp2_region,$np1,$np2,$np1_orientation,$np2_orientation,$name,$bp1_chr,$bp2_chr);
		
    }
    print qq(see the result printed in $name.translocation.info.txt\n);
    close OUT;
}
##########################################

sub get_translocation_info {
    my $translocation_boundary;
    my ($transcript_info,$translocation,$region1,$region2,$regionside,$organism,$bp1_gene,$bp1_strand,$bp1_transcript,$bp1_flip,$bp1_region,$bp2_gene,$bp2_strand,$bp2_transcript,$bp2_flip,$bp2_region,$np1,$np2,$np1_orientation,$np2_orientation,$name,$bp1_chr,$bp2_chr) = @_;
    my $bases_per_side=0;
    my $x = 0;
    my @translocation_dnaseq;
    
    my $frame = 0;
    my @codon;
    my $myCodonTable = Bio::Tools::CodonTable->new();
    my $stop_codon;
    my $stopped_short = 0;
    
    my $bp1exons;
    my $bp2exons;
    
    my ($tps1_exon,$tps2_exon);
    foreach my $n (sort {$a<=>$b} keys %{$translocation}) {
	
	if ($stop_codon) {
	    $stopped_short++;
	} else {
	    
	    #trv_type {3_prime_flanking_region,3_prime_untranslated_region,5_prime_flanking_region,5_prime_untranslated_region,intronic,missense,nonsense,silent,splice_site}
	    ###if ($translocation_boundary == $n) {print qq( boundary );}
	    
	    my $base = $translocation->{$n}->{base};
	    my $transcript = $translocation->{$n}->{transcript};
	    my ($exon,$region,$side) = split(/\,/,$translocation->{$n}->{exon});

	    #unless ($regionside eq $side) {if ($side eq "BP2") {$side = "BP1";} elsif ($side eq "BP1") {$side = "BP2";}}
		
	    if ($base =~ /\d/) {	
		if ($side eq "BP2") {
		    if ($region1 =~ /cds/ && $region2 =~ /5_prime_untranslated/) {
			my ($gpos,$strand) = split(/\:/,$base);
			my $refbase = &get_ref_base($gpos,$gpos,$bp2_chr,$organism);
			if ($strand eq "-1") { my $revrefbase = &reverse_complement_allele($refbase);$refbase=$revrefbase;}
			$base = $refbase;
		    }

		    if ($region =~ /3_prime_untranslated/) {
			my ($gpos,$strand) = split(/\:/,$base);
			my $refbase = &get_ref_base($gpos,$gpos,$bp2_chr,$organism);
			if ($strand eq "-1") { my $revrefbase = &reverse_complement_allele($refbase);$refbase=$revrefbase;}
			$base = $refbase;
		    }
		}
	    }
	    
	    #my $trans_pos = $translocation->{$base_count}->{trans_pos};
	    
	    my $trans_pos = $translocation->{$n}->{trans_pos};
	    if ($trans_pos) {
		print OUT qq(bases to the break point $bases_per_side\n);
		print qq(bases to the break point $bases_per_side\n);
		$translocation_boundary = $bases_per_side;
		$bases_per_side = 0;

		my ($tpe,$tps) = split(/\,/,$trans_pos);
		if ($tps eq "BP1" && $side eq "BP1" && $region1 =~ /cds/ && $region2 =~ /intron/ && $tpe == $exon) {
		    $tps1_exon = $exon;
		}
		if ($tps eq "BP2" && $side eq "BP2" && $region1 =~ /intron/ && $region2 =~ /cds/ && $tpe == $exon) {
		    $tps2_exon = $exon;
		}
	    }

	    if ($side eq "BP1" && $tps1_exon) {
		if ($tps1_exon eq $exon) {
		    $base = 1;
		}
	    }
	    if ($side eq "BP2" && $tps2_exon) {
		if ($tps2_exon eq $exon) {
		    $base = 1;
		}
	    }
	    
	    unless ($base =~ /\d/) {
		$bases_per_side++;

		if ($transcript eq $bp1_transcript) {
		    $bp1exons->{$exon}=1;
		    #if ($bp1_flip eq "YES") {
		#	my $revbase = &reverse_complement_allele($base);    ###Flip doesn't mean reverse complement
		#	$base=$revbase;
		    #}
		}
		if ($transcript eq $bp2_transcript) {
		    $bp2exons->{$exon}=1;
		    #if ($bp2_flip eq "YES") {
		#	my $revbase = &reverse_complement_allele($base);    ###
		#	$base=$revbase;
		    #}
		}
		
		#print qq($base);
		#unless ($region =~ /utr/){
		push(@translocation_dnaseq,$base);
		
		$x++;
		
		$frame++;
		if ($frame == 3) {
		    push(@codon,$base);
		    my $codon = join '' , @codon;
		    my $aa = $myCodonTable->translate($codon);
		    if ($aa eq "*") {
			$stop_codon = $aa;
		    }
		    $frame = 0;
		    undef(@codon);
		} else {
		    push(@codon,$base);
		}
	    }
	}
    }

    print OUT qq(bases from the break point $bases_per_side\n);
    print qq(bases from the break point $bases_per_side\n);
    print qq(\n);

    my $translocation_seq = join '',@translocation_dnaseq;
    
    my $protien_seq = $myCodonTable->translate($translocation_seq);
#my @pcodons = $myCodonTable->revtranslate($daa);
    
#$transcript_info->{$bp1_transcript}->{$nbp1}->{trans_pos}
    
    print OUT qq(Breakpoint 1; Gene=>$bp1_gene Transcript=>$bp1_transcript Chromosome=>$bp1_chr Coordinate=>$np1 Region=>$bp1_region Incorporated exons =>{);
    print qq(Breakpoint 1; Gene=>$bp1_gene Transcript=>$bp1_transcript Chromosome=>$bp1_chr Coordinate=>$np1 Region=>$bp1_region Incorporated exons =>{);
    foreach my $exon (sort {$a<=>$b} keys %{$bp1exons}) {
	print OUT qq( $exon);
	print qq( $exon);
    }
    my $bp1exon_total = $transcript_info->{$bp1_transcript}->{-1}->{exon_total};
    print OUT qq( } Total exons in the transcript =>$bp1exon_total\n);
    print qq( } Total exons in the transcript =>$bp1exon_total\n);

    print OUT qq(Breakpoint 2; Gene=>$bp2_gene Transcript=>$bp2_transcript Chromosome=>$bp2_chr Coordinate=>$np2 Region=>$bp2_region Incorporated exons =>{);
    print qq(Breakpoint 2; Gene=>$bp2_gene Transcript=>$bp2_transcript Chromosome=>$bp2_chr Coordinate=>$np2 Region=>$bp2_region Incorporated exons =>{);
    foreach my $exon (sort {$a<=>$b} keys %{$bp2exons}) {
	print OUT qq( $exon);
	print qq( $exon);
    }
    my $bp2exon_total = $transcript_info->{$bp2_transcript}->{-1}->{exon_total};
    print OUT qq( } Total exons in the transcript =>$bp2exon_total\n);
    print qq( } Total exons in the transcript =>$bp2exon_total\n);
    unless ($stopped_short == 0) {
	print OUT qq(In translation, the transcript hit a stop codon $stopped_short bases before the end of $bp2_transcript\n);
	print qq(In translation, the transcript hit a stop codon $stopped_short bases before the end of $bp2_transcript\n);
    }
    
    print OUT qq(\n\nHere is the new sequence that describes the translocation with the splice at base $translocation_boundary\n\>$name.nucleotide.seq\n$translocation_seq);
    print qq(\n\nHere is the new sequence that describes the translocation with the splice at base $translocation_boundary\n\>$name.nucleotide.seq\n$translocation_seq);

    
    my $looks;
    if (($protien_seq =~ /^M/) && ($protien_seq =~ /\*/)) {
	$looks = "The protien represented by $name starts with a Methionine and has a stop codon";
    } else {
	$looks = "The protien represented by $name either doesn't start with a Methionine or doesn't have a stop codon";
    }
    print OUT  qq(\n\n$looks\nHere is the protien sequence\n>$name.protien.seq\n$protien_seq\n\n);
    print qq(\n\n$looks\nHere is the protien sequence\n>$name.protien.seq\n$protien_seq\n\n);

    
    print OUT qq($bp1_gene $bp1_transcript flipped $bp1_flip gene strand $bp1_strand, alignment strand $np1_orientation. 
$bp2_gene $bp2_transcript flipped $bp2_flip gene strand $bp2_strand, alignment strand $np2_orientation.\n);
    
    print qq($bp1_gene $bp1_transcript flipped $bp1_flip gene strand $bp1_strand, alignment strand $np1_orientation. 
$bp2_gene $bp2_transcript flipped $bp2_flip gene strand $bp2_strand, alignment strand $np2_orientation.\n);
    
    my $y = $x/3;
    print OUT qq(\nthe new transcript has $x bases and $y amino acids\n\n);
    print qq(\nthe new transcript has $x bases and $y amino acids\n\n);


}


sub get_translocation {
    
    my ($pause,$resume,$bp1_region,$bp2_region);
    my ($translocation,$transcript_info,$gene,$strand,$transcript,$base_count,$side,$flip) = @_;

    #if (($strand eq "+1" && $flip eq "NO") || ($strand eq "-1" && $flip eq "YES")) {

    my (@positions);

    if ($strand eq "+1") {
	foreach my $pos (sort {$a<=>$b} keys %{$transcript_info->{$transcript}}) {
	    unless ($pos == -1) {push(@positions,$pos);}
	}
    } else {
	foreach my $pos (sort {$b<=>$a} keys %{$transcript_info->{$transcript}}) {
	    unless ($pos == -1) {push(@positions,$pos);}
	}
    }

    for my $pos (@positions) {
	
	my ($exon,$region) = split(/\,/,$transcript_info->{$transcript}->{$pos}->{exon});
	my $frame = $transcript_info->{$transcript}->{$pos}->{frame};
	my $aa_n = $transcript_info->{$transcript}->{$pos}->{aa_n};
	my $base = $transcript_info->{$transcript}->{$pos}->{base};
	my ($trans_pos) = $transcript_info->{$transcript}->{$pos}->{trans_pos};
	#my ($trans_pos,$region) = split(/\,/,$transcript_info->{$transcript}->{$pos}->{trans_pos});
	
	if ($trans_pos) {
	    unless ($trans_pos =~ /intron/) {
		if ($side eq "BP2") {
		    $resume = 1;
		    ($bp2_region) = $trans_pos =~ /\S+\,(\S+)/;
		    $translocation->{$base_count}->{trans_pos}="$exon,BP2";
		    
		    #$bp2_region = $region;
		}
	    }
	}
	if ($side eq "BP2") {
	    if ($resume) {
		
		#unless ($trans_pos =~ /utr/) {
		#print qq($gene $transcript $exon $frame $aa_n $base $pos $bp2_chr\n);
		$base_count++;
		$translocation->{$base_count}->{base}=$base;
		$translocation->{$base_count}->{transcript}=$transcript;
		$translocation->{$base_count}->{exon}="$exon,$region,BP2";
		#$translocation->{$base_count}=$base;
		#}
	    }
	}
	
	if ($side eq "BP1") {
	    unless ($pause) {
		#print qq($gene $transcript $exon $frame $aa_n $base $pos $bp1_chr\n);
		$base_count++;
		$translocation->{$base_count}->{base}=$base;
		$translocation->{$base_count}->{transcript}=$transcript;
		$translocation->{$base_count}->{exon}="$exon,$region,BP1";
		#$translocation->{$base_count}=$base;
	    }
	}
	if ($trans_pos) {
	    if ($side eq "BP1") {
		$pause=1;
		($bp1_region) = $trans_pos =~ /\S+\,(\S+)/;
		$translocation->{$base_count}->{trans_pos}="$exon,BP1";
		#$bp1_region = $region;
		
	    }
	    if ($trans_pos =~ /intron/) {
		if ($side eq "BP2") {
		    $resume = 1;
		    ($bp2_region) = $trans_pos =~ /\S+\,(\S+)/;
		    $translocation->{$base_count}->{trans_pos}="$exon,BP2";
		    #$bp2_region = $region;
		    
		}
	    }
	}
    }
    my $return_region;
    if ($side eq "BP2") {$return_region=$bp2_region;}else{$return_region=$bp1_region;}
    return($translocation,$base_count,$return_region);
}



sub parse_annotation {
    my $annotation_info;
    my ($annotation_file) = @_;
    open(ANO,"$annotation_file");
    while (<ANO>) {
	chomp;
	my $line = $_;
	my ($chromosome,$start,$stop,$ref,$var,$variant_type,$gene,$transcript,$transcript_species,$source,$tv,$strand,$Traans_stat,$trv_type,$c_pos,$aa,$cons_score,$domain) = split(/[\s]+/,$line); ##get extended annotation
#chromosome_name	start	stop	reference	variant	type	gene_name	transcript_name	transcript_species	transcript_source	transcript_version	strand	transcript_status	trv_type	c_position	amino_acid_change	ucsc_cons	domain	all_domains	flank_annotation_distance_to_transcript	intron_annotation_substructure_ordinal	intron_annotation_substructure_size	intron_annotation_substructure_position

	unless ($chromosome eq "chromosome_name") {
	    my $chr = "chr$chromosome";

	    if ($transcript =~ /\d+/) {
		#if ($trv_type =~ /untranslated/) {
		    #print qq($chr $start falls in the $trv_type of $gene $transcript and will effectively act as no transcript id\'ed for this site\n);
		#} else {

		    #if ($transcript eq "ENST00000370769" || $transcript eq "NM_002393") {
		    $annotation_info->{$chr}->{$start}->{gene}=$gene;
		    $annotation_info->{$chr}->{$start}->{transcript}=$transcript;
		    $annotation_info->{$chr}->{$start}->{strand}=$strand;
		    $annotation_info->{$chr}->{$start}->{trv_type}=$trv_type;
		    #}
		    #print qq($chr,$start,$stop,$ref,$gene,$transcript,$strand,$trv_type\n);
		#}
	    }
	}
    }
    return ($annotation_info);
}


sub parse_blast_out {
    my $locations;
    my $alignment_locations;
    my ($blast_out,$bp1_chr,$bp1_pos,$bp2_chr,$bp2_pos) = @_;
    open(IN,$blast_out);
    while (<IN>) {
	chomp;
	my $line = $_;
	my ($subject_id,$query_id,$subject_cov,$query_cov,$percent_identity,$bit_score,$p_value,$subject_length,$query_length,$alignment_bases,$HSPs) = split(/\t/,$line);

	$HSPs =~ s/\s//gi;
	if ($subject_id eq $bp1_chr) {
	    $locations=&get_locations($HSPs,$bp1_chr,$bp1_pos,$locations);
	}
	if ($subject_id eq $bp2_chr) {
	    $locations=&get_locations($HSPs,$bp2_chr,$bp2_pos,$locations);
	}
    }
    foreach my $chr (sort keys %{$locations}) {
	foreach my $apos (sort keys %{$locations->{$chr}}) {
	    my $anchore;
	    my $anchore_alignment;
	    foreach my $alignment (sort keys %{$locations->{$chr}->{$apos}}) {
		my $pos = $locations->{$chr}->{$apos}->{$alignment};
		unless ($anchore) {$anchore = $pos;$anchore_alignment=$alignment;}
		if (abs($apos - $pos) < abs($apos - $anchore)) {
		    $anchore = $pos;
		    $anchore_alignment = $alignment;
		}
	    }
	    #print qq($chr $apos => $anchore $anchore_alignment\n);
	    my $new_loc = "$anchore,$anchore_alignment";
	    $alignment_locations->{$chr}->{$apos}=$new_loc;
	}
    }
    return($alignment_locations);
}

sub get_locations {
    my ($HSPs,$chr,$pos,$locations) = @_;

    my @hsps = split(/\;/,$HSPs);

    for my $loc (@hsps) {

	my ($start,$stop) = $loc =~ /(\d+)\-(\d+)\:/;
	my $end;
	if ($pos eq $start) {
	    $end = $start;
	} elsif ($pos eq $stop) {  
	    $end = $stop;
	} elsif (abs($pos - $start) < abs($pos - $stop)) {
	    $end = $start;
	} else {
	    $end = $stop;
	}
		
	#print qq($chr $pos $end $loc\n);
	my $alignment = "$start-$stop";
	$locations->{$chr}->{$pos}->{$loc}=$end;

    }
    return($locations);
}


sub get_ref_base {

    my ($chr_start,$chr_stop,$chr_name,$organism) = @_;

    my $RefDir;
    if ($organism eq "human") {
	$RefDir = "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/";
    } else {
	$RefDir = "/gscmnt/sata147/info/medseq/rmeyer/resources/MouseB37/";
    }

    my $refdb = Bio::DB::Fasta->new($RefDir);

    my $seq = $refdb->seq($chr_name, $chr_start => $chr_stop);
    $seq =~ s/(\S)/\U$1/;

    return $seq;
}

sub ipc_run {

    my (@command) = @_;

    #my $comand_line = join (' ' , @command);
    #print qq(will now run $comand_line\n);

    my ($in, $out, $err);
    my ($obj) = IPC::Run::run(@command, \$in, \$out, \$err);
    if ($err) {
#	print qq($err\n);
    }
    if ($out) {
	return ($out);
	#print qq($out\n);
    }
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

sub get_alignment_orientation {

    my ($start,$stop) = @_;

    my $orientation;
    if ($start <= $stop) {
	$orientation = "plus";
    } else {
	$orientation = "minus";
    }
    return ($orientation);

}

sub parse_breakdancer {
    
    my ($bp1_chr,$bp1_pos,$bp2_chr,$bp2_pos,$np1_orientation,$np2_orientation,$name);
    my ($breakdanerfile) = @_;
    
    my ($bp1_c,$bp1_p,$bp1read_ori,$bp2_c,$bp2_p,$bp2read_ori);
    if (-f $breakdanerfile) {
	open(BDO,"$breakdanerfile") || die "could not open the $breakdanerfile file\nWill exit the program now\n";
	while (<BDO>) {
	    chomp;
	    my $line = $_;
#15      72113493        7+0-    17      35752361        0+2-  CTX     -267    36      2       tB|2    1.00  BreakDancer0.0.1r47    t1
	    ($bp1_c,$bp1_p,$bp1read_ori,$bp2_c,$bp2_p,$bp2read_ori) = (split(/[\s]+/,$line))[0,1,2,3,4,5];
	}
	close (BDO);
    } else {
	($bp1_c,$bp1_p,$bp1read_ori,$bp2_c,$bp2_p,$bp2read_ori) = (split(/[\s]+/,$breakdanerfile))[0,1,2,3,4,5];
    }
    unless ($bp1_c && $bp1_p && $bp1read_ori && $bp2_c && $bp2_p && $bp2read_ori) { die "could not identify the breakdancer info\nWill exit the program now\n"; }



    $bp1_chr = "chr$bp1_c";
    $bp2_chr = "chr$bp2_c";
    $bp1_pos = $bp1_p;
    $bp2_pos = $bp2_p;
    
    my ($support_bp1plus,$support_bp1minus) = $bp1read_ori =~ /([\d]+)\+([\d]+)\-/;
    if ($support_bp1plus >= $support_bp1minus) {
	$np1_orientation = "plus";
    } else {
	$np1_orientation = "minus";
    }
    my ($support_bp2plus,$support_bp2minus) = $bp2read_ori =~ /([\d]+)\+([\d]+)\-/;
    if ($support_bp2plus >= $support_bp2minus) {
	$np2_orientation = "minus";
	#$np2_orientation = "plus";
    } else {
	#$np2_orientation = "minus";
	$np2_orientation = "plus";
    }
    $name = "$bp1_chr\_$bp1_pos\-$bp2_chr\_$bp2_pos";
    
    return($bp1_chr,$bp1_pos,$bp2_chr,$bp2_pos,$np1_orientation,$np2_orientation,$name);
    
}


sub compare_orientations {

    my $flip;
    my ($transcript_strand,$alignment_orientation) = @_;

    if (($transcript_strand eq "+1" && $alignment_orientation eq "plus") || ($transcript_strand eq "-1" && $alignment_orientation eq "minus")) {
	$flip = "NO";
    } else {
	$flip = "YES";
    }

    #print qq($transcript_strand,$alignment_orientation,$flip\n);
    return($flip);
}


1;

### test case for a fusion gene => gmt annotate translocation-information -b "15       72113189       16+1-   17       35756596        14+18-  CTX" 

### same test case for a fusion gene but with a fasta => gmt annotate translocation -f /gscmnt/sata147/info/medseq/rmeyer/PROJECTS/TRANSLOCATIONS/chr15_72113189-chr17_35756596-CTX.fasta

 
=head1 TITLE

get_translocation_information

=head1 DESCRIPTION

This script can deal with a variety of inputs and dump out many fun and exciting files for your viewing pleasure!

=head1 Input Options:

a file containing one line of breakdancer output
and or
a contig fasta file representing the consensus of the reads forming the translocation

=head1 KNOWN BUGS

Please report bugs to <rmeyer@genome.wustl.edu>

=head1 AUTHOR

Rick Meyer <rmeyer@genome.wustl.edu>

=cut
