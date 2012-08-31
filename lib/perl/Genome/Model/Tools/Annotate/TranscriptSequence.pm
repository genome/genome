package Genome::Model::Tools::Annotate::TranscriptSequence;

use strict;
use warnings;
use Genome;

use Bio::Seq;
use Bio::Tools::CodonTable;
use Bio::DB::Fasta;

class Genome::Model::Tools::Annotate::TranscriptSequence {
    is => 'Command',                       
    has => [ 
	transcript => {
	    type  =>  'String',
	    doc   =>  "provide the transcript name",
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
	trans_pos => {
	    type  =>  'String',
	    doc   =>  "provide a coordinate of interest",
	    is_optional  => 1,
	},
	utr_seq => {
	    is => 'Boolean',
	    doc   =>  "use this flag if you would like to retriev the utr sequence for this transcript.",
	    is_optional  => 1,
	    default => 0,
	},
	masked => {
	    is => 'Boolean',
	    doc   =>  "use this option to mask_snps_and_repeats",
	    is_optional  => 1,
	},
	   output => {
	    type  =>  'String',
	    doc   =>  "provide a file name to write you transcript information to .txt will be appended to it. Default is to print to stdout.",
	    is_optional  => 1,
	},
	   fasta => {
	    type  =>  'Boolean',
	    doc   =>  "will write the nucleotide fasta to output.fasta. Default is to print to stdout.",
	    is_optional  => 1,
	},
	no_stdout => {
	    is => 'Boolean',
	    doc   =>  "Use this option if you do not want the info to print to stdout. Default is to print to stdout.",
	    is_optional  => 1,
	},


    ], 
};


sub help_brief {
    return <<EOS
	gmt annotate transcript-sequence will print information about and transcript in the lastest version of annotation data in our database
EOS
}

sub help_synopsis {
    return <<EOS
	gmt annotate transcript-sequence -transcript NM_001024809

	or for multiple transcripts, seperate each with a comma and if the trans-pos option is used the positions need to be in the order of the transcripts

	gmt annotate transcript-sequence -transcript NM_001024809,NM_033238 -trans-pos 35752360,72113517
EOS
}

sub help_detail { 
    return <<EOS 
	
	-trans_pos option will locate the given coordinate in the transcript and print a line of information at the bottom of your screen 
	-organism use this option if your transcript is from mouse otherwise, human will be assumed
	-version if you would prefer something other than the default for human is 54_36p_v2 and for mouse is 54_37g_v2
	-utr_seq will display the utr sequence for human build 36 or mouse build 37

	gmt annotate transcript-sequence -transcript ENSMUST00000102781 -organism mouse -utr-seq -trans-pos 100857095

EOS
}


sub execute {
    my ($transcript_info);

    my $self = shift;
    my $output = $self->output;

    if ($output) {
	open(OUT,">$output.txt") || $self->error_message( "\n\ncouldn't open the output file $output.txt\n\n") && return; 
	if ($self->fasta) {
	    open(FASTA,">$output.fasta") || $self->error_message( "\n\ncouldn't open the output file $output.fasta\n\n") && return;
	}
    }
    
    my $trans_poss = $self->trans_pos;
    #my $transcript = $self->transcript;
    my $transcripts = $self->transcript;

    my @trans = split(/\,/,$transcripts);
    my $n = 0;
    
    for my $transcript (@trans) {

	if ($n > 0) { unless ($self->no_stdout) { print qq(\n\n\n); }}

	my $trans_pos;
	if ($trans_poss) {
	    $trans_pos = (split(/\,/,$trans_poss))[$n];
	}
	$n++;
	my ($chromosome,$strand);
	($transcript_info,$chromosome,$strand) = &get_transcript_info($self,$transcript,$transcript_info,$trans_pos);
	my ($transcript_info,$trans_pos_number_line,$post_pos_bases) = &get_transcript($self,$transcript,$trans_pos,$transcript_info,$strand,$chromosome);
	$transcript_info->{$transcript}->{-1}->{strand}=$strand;

	if ($self->trans_pos) {
	    if ($trans_pos_number_line) {
		my $line = qq(\n$trans_pos_number_line. There are $post_pos_bases coding bases left in $transcript after $trans_pos\n);
		if ($self->output) {print OUT qq($line);}
		unless ($self->no_stdout) {print qq($line);}
		
		$transcript_info->{$transcript}->{-1}->{post_pos_bases}=$post_pos_bases;

	    } else {
		$transcript_info->{$transcript}->{-1}->{trans_posid}="not_ided";
		my $line = qq(\n$trans_pos was not located\n);
		if ($self->output) {print OUT qq($line);}
		unless ($self->no_stdout) {print qq($line);}

	    }
	}
    }
    if ($self->fasta) {
	unless ($self->no_stdout) {print qq(\n\n);}
	for my $transcript (@trans) {
	    my $fasta = $transcript_info->{$transcript}->{-1}->{FASTA};
	    if ($self->output) {
		print FASTA qq($fasta\n);
	    } else {
		unless ($self->no_stdout) {print qq($fasta\n);}
	    }
	}
    }
    
    close (OUT);
    return($transcript_info);
}

sub get_transcript {
	
    my $trans_pos_number_line;
    my ($self,$transcript,$trans_pos,$transcript_info,$strand,$chromosome) = @_;
    #my $transcript = $self->transcript;
    #my $trans_pos = $self->trans_pos;
    unless ($trans_pos) { $trans_pos = 0; }

    my $organism = $self->organism;

    my $version = $self->version;
    if ($organism eq "mouse" && $version eq "54_36p_v2") { $version = "54_37g_v2"; }

    my $utr_seq = $self->utr_seq;
    if ($utr_seq && $organism eq "human") { unless ($version =~ /\_36/) { print qq(can only get utr seq for human build 36);$utr_seq = 0; } }
    if ($utr_seq && $organism eq "mouse") { unless ($version =~ /\_37/) { print qq(can only get utr seq for mouse build 37);$utr_seq = 0; } }

    my $myCodonTable = Bio::Tools::CodonTable->new();
    my @seq;
    my @maskedseq;
    my @fullseq;
    my (@p3m,@p5m,@p3,@p5);

    my ($pexon,$pregion,$ppos);
    my $p5=0;
    my $p3=0;
    my ($coding_start,$coding_stop);

    #print qq($transcript $strand\n);
    
    my @positions;
    my $trans_pos_in=0;
    my $trans_pos_in_5utr=0;
    my $trans_pos_in_3utr=0;
    if ($strand eq "+1") {
	foreach my $pos (sort {$a<=>$b} keys %{$transcript_info->{$transcript}}) {
	    unless ($pos == -1) {push(@positions,$pos);}
	}
    } else {
	foreach my $pos (sort {$b<=>$a} keys %{$transcript_info->{$transcript}}) {
	    unless ($pos == -1) {push(@positions,$pos);}
	}
    }
    my $pre_pos_bases = 0;
    for my $pos (@positions) {
	my ($exon,$region) = split(/\,/,$transcript_info->{$transcript}->{$pos}->{exon});
	my $frame = $transcript_info->{$transcript}->{$pos}->{frame};
	my $aa_n = $transcript_info->{$transcript}->{$pos}->{aa_n};
	my $base = $transcript_info->{$transcript}->{$pos}->{base};
	my $masked_base;
	if ($self->masked) {$masked_base = $transcript_info->{$transcript}->{$pos}->{masked_base};}

	my ($trans_posid) = $transcript_info->{$transcript}->{$pos}->{trans_pos};
	
	if ($trans_posid) {
	    $pre_pos_bases = $trans_pos_in;
	    my ($trans_pos_n,$trans_pos_r)  =  split(/\,/,$trans_posid);
	    $trans_pos_number_line = qq(The position $trans_pos_n was in an $trans_pos_r and is in or after $region $exon, frame $frame, base $base, amino_acid_numuber $aa_n, and falls after $trans_pos_in_5utr bases of 5prime UTR, $trans_pos_in bases of coding sequence and $trans_pos_in_3utr bases of 3prime UTR);

	    $transcript_info->{$transcript}->{-1}->{trans_pos_n}=$trans_pos_n;
	    $transcript_info->{$transcript}->{-1}->{trans_pos_r}=$trans_pos_r;
	    $transcript_info->{$transcript}->{-1}->{region}=$region;
	    $transcript_info->{$transcript}->{-1}->{exon}=$exon;
	    $transcript_info->{$transcript}->{-1}->{frame}=$frame;
	    $transcript_info->{$transcript}->{-1}->{base}=$base;
	    $transcript_info->{$transcript}->{-1}->{aa_n}=$aa_n;
	    $transcript_info->{$transcript}->{-1}->{trans_pos_in_5utr}=$trans_pos_in_5utr;
	    $transcript_info->{$transcript}->{-1}->{trans_pos_in}=$trans_pos_in;
	    $transcript_info->{$transcript}->{-1}->{trans_pos_in_3utr}=$trans_pos_in_3utr;
	    $transcript_info->{$transcript}->{-1}->{trans_posid}=$trans_posid;

	}

	if ($region =~ /cds/) {
	    $trans_pos_in++;
	    unless ($coding_start) {$coding_start=$pos;}
	    $coding_stop=$pos;
	}
	
	#if ($base =~ /\d/) {

	if ($region =~ /utr/) {
	    if ($coding_start) {
		$p3++;
		$trans_pos_in_3utr++;
		push(@p3,$base);
		if ($self->masked) {push(@p3m,$masked_base);}
		
	    } else {
		$p5++;
		$trans_pos_in_5utr++;
		push(@p5,$base);
		if ($self->masked) {push(@p5m,$masked_base);}
	    }
	}
	#} else {
	unless ($region =~ /utr/) { 
	    push(@seq,$base);
	    if ($self->masked) {push(@maskedseq,$masked_base);}
	#}
	}
	
	my $range = $transcript_info->{$transcript}->{$pos}->{range};
	my ($r_start,$r_stop) = split(/\-/,$range);
	if ($pos == $r_stop) {
	    if ($region =~ /utr/) {
		if ($coding_start) {
		    my $line = qq($exon $region $range $p3\n);
		    if ($self->output) {print OUT qq($line);}
		    unless ($self->no_stdout) {print qq($line);}
		    $p3=0;
		    if ($utr_seq) { 
			my $us = join '' , @p3;
			unless ($self->no_stdout) {print qq($us\n\n);}
			if ($self->output) {print OUT qq($us\n\n);}
			undef(@p3);
			if ($self->masked) {
			    my $mus = join '' , @p3m;
			    unless ($self->no_stdout) {print qq($mus\n\n);}
			    if ($self->output) {print OUT qq($mus\n\n);}
			    undef(@p3m);
			}
			#&print_utr_seq($r_start,$r_stop,$organism,$strand,$chromosome,$self);
		    }
		} else {
		    my $line = qq($exon $region $range $p5\n);
		    if ($self->output) {print OUT qq($line);}
		    unless ($self->no_stdout) {print qq($line);}
		    $p5=0;
		    if ($utr_seq) { 
			my $us = join '' , @p5;
			unless ($self->no_stdout) {print qq($us\n\n);}
			if ($self->output) {print OUT qq($us\n\n);}
			undef(@p5);
			if ($self->masked) {
			    my $mus = join '' , @p5m;
			    unless ($self->no_stdout) {print qq($mus\n\n);}
			    if ($self->output) {print OUT qq($mus\n\n);}
			    undef(@p5m);
			}
			#&print_utr_seq($r_start,$r_stop,$organism,$strand,$chromosome,$self); }
		    }
		}
	    }
	    if ($region =~ /cds/) {
		my $cds = join '' , @seq;
		my $length = length($cds);

		my $line =  qq($exon $region $range $length\n$cds\n\n);
		if ($self->output) {print OUT qq($line);}
		unless ($self->no_stdout) {print qq($line);}

		my $maskedcds;
		if ($self->masked) {
		    $maskedcds = join '' , @maskedseq;
		    unless ($self->no_stdout) {print qq($maskedcds\n\n);}
		    if ($self->output) {print OUT qq($maskedcds\n\n);}
		}

		push(@fullseq,$cds);
		undef(@seq);
		undef(@maskedseq);
	    }
	}
    }
    my $post_pos_bases = $trans_pos_in - $pre_pos_bases - 1;

    my $sequence = join '' , @fullseq;
    my $protien_seq = $myCodonTable->translate($sequence);

    $transcript_info->{$transcript}->{-1}->{protien_seq} = $protien_seq;
    $transcript_info->{$transcript}->{-1}->{sequence} = $sequence;

    my $hugo_gene_name = $transcript_info->{$transcript}->{-1}->{hugo_gene_name};
    
    if ($self->fasta) {
	if ($sequence) {
	    $transcript_info->{$transcript}->{-1}->{FASTA} = "\>$transcript\|$strand\|Chromosome:$chromosome\|$hugo_gene_name\|\n$sequence";
	} else {
	    $transcript_info->{$transcript}->{-1}->{FASTA} = "\>$transcript\|$strand\|Chromosome:$chromosome\|$hugo_gene_name\|";
	}
    }

    my $line = qq(\n\>$transcript.dna.fasta\n$sequence\n\n\>$transcript.protien.fasta\n$protien_seq\n);
    if ($self->output) {print OUT qq($line);}
    unless ($self->no_stdout) {print qq($line);}

    return ($transcript_info,$trans_pos_number_line,$post_pos_bases);

}

sub print_utr_seq {
    
    my ($r_start,$r_stop,$organism,$strand,$chromosome,$self) = @_;
    if ($strand eq "-1") {
	my $seq = &get_ref_base($r_stop,$r_start,$chromosome,$organism);
	my $rev = &reverse_complement_allele($seq);
	my $line = qq($rev\n\n);
	if ($self->output) {print OUT qq($line);}
	unless ($self->no_stdout) {print qq($line);}
    } else {
	my $seq = &get_ref_base($r_start,$r_stop,$chromosome,$organism);
	my $line = qq($seq\n\n);
	if ($self->output) {print OUT qq($line);}
	unless ($self->no_stdout) {print qq($line);}
    }
}

sub get_transcript_info {
    my ($strand,$chromosome);
    my ($self,$transcript,$transcript_info,$trans_pos) = @_;
    #($transcript_info) = @_;
    #my $self = shift;
    #my $transcript = $self->transcript;
    #my $trans_pos = $self->trans_pos;
    unless ($trans_pos) { $trans_pos = 0; }

    my $organism = $self->organism;

    my $version = $self->version;
    if ($organism eq "mouse") { if ($version eq "54_36p_v2") { $version = "54_37g_v2";}}
    my $genome;
    if ($self->masked) {
	if ($organism eq "human") {
	    $genome = GSC::Sequence::Genome->get(sequence_item_name => 'NCBI-human-build36');
	} else {
	    $genome = GSC::Sequence::Genome->get(sequence_item_name => 'NCBI-mouse-buildC57BL6J');
	}
    }



    my ($ncbi_reference) = $version =~ /\_([\d]+)/;

    my $eianame = "NCBI-" . $organism . ".ensembl";
    my $gianame = "NCBI-" . $organism . ".genbank";
    my $build_source = "$organism build $ncbi_reference version $version";

    my $ensembl_build = Genome::Model::ImportedAnnotation->get(name => $eianame)->build_by_version($version);
    my ($ensembl_data_directory)= $ensembl_build->determine_data_directory;
    my $genbank_build = Genome::Model::ImportedAnnotation->get(name => $gianame)->build_by_version($version);
    my ($genbank_data_directory) = $genbank_build->determine_data_directory;

    my $t;
    if ($transcript =~/^ENS/){ #ENST for Human ENSMUST
	($t) = Genome::Transcript->get( transcript_name =>$transcript, data_directory => $ensembl_data_directory, reference_build_id => $ensembl_build->reference_sequence_id);
    }else{
	($t) = Genome::Transcript->get( transcript_name =>$transcript, data_directory => $genbank_data_directory, reference_build_id => $genbank_build->reference_sequence_id);
    }
    unless ($self->no_stdout) {print qq($ensembl_data_directory\n\n\n);}
    unless ($t) {print qq(\nCould not find a transcript object for $transcript from the $organism data warehouse\nWill exit the program now\n\n);;exit(1);}

    my $tseq = $t->cds_full_nucleotide_sequence;
    my @substructures = $t->ordered_sub_structures;
    
    my $total_substructures = @substructures;
    my $t_n = 0; #substructure counter
    
    $strand = $t->strand;
    $chromosome = $t->chrom_name;
    $transcript_info->{$transcript}->{-1}->{chromosome}=$chromosome;
    my $info;
    $info->{$transcript}->{strand}=$strand;

    my $data_directory = $t->data_directory;
    my $gene_id = $t->gene_id;
    my $source = $t->source;
    my $transcript_status = $t->transcript_status;

    my $gene = $t->gene;
    my $hugo_gene_name = $gene->hugo_gene_name;
    unless ($hugo_gene_name) {$hugo_gene_name = "unlisted";}

    my $line = qq(Hugo gene name $hugo_gene_name, Gene Id $gene_id, Transcript name $transcript, Chromosome $chromosome, Strand $strand, Transcript status $transcript_status, Transcript source $source $build_source\n\n\n);
    if ($self->output) {print OUT qq($line);}
    unless ($self->no_stdout) {print qq($line);}

    $transcript_info->{$transcript}->{-1}->{source_line} = qq(Hugo gene name $hugo_gene_name, Gene Id $gene_id, Transcript name $transcript, Chromosome $chromosome, Strand $strand, Transcript status $transcript_status, Transcript source $source $build_source);
    $transcript_info->{$transcript}->{-1}->{hugo_gene_name} = $hugo_gene_name;
    $transcript_info->{$transcript}->{-1}->{gene_id} = $gene_id;

    if (@substructures) {
	#print qq($transcript $total_substructures  $strand  $chr $trans_pos\n);

	while ($t_n < $total_substructures) {
	    my $t_region = $substructures[$t_n];
	    $t_n++;
	    
	    my $tr_start = $t_region->{structure_start};
	    my $tr_stop = $t_region->{structure_stop};
	    my $range = "$tr_start\-$tr_stop";
	    my $structure_type = $t_region->{structure_type};

	    #print qq($structure_type $range\n);

	    if ($t_region->{structure_type} =~ /exon/) {
		unless ($transcript_info->{$transcript}->{-1}->{first_base}) {
		    $transcript_info->{$transcript}->{-1}->{first_base}=$tr_start;
		}
		$transcript_info->{$transcript}->{-1}->{last_base}=$tr_stop;

		my $trv_type = $t_region->{structure_type};
		my @nucleotide_array = split(//,$t_region->nucleotide_seq);
		
		my $base_n;
		if ($strand eq "-1") { $base_n=@nucleotide_array; } else {$base_n=-1;}
		
		my @masked_array;
		if ($self->masked) {
		    #print qq(looking for mask_snps_and_repeats $chromosome $tr_start $tr_stop\n);
		    
		    my $chr = $genome->get_chromosome($chromosome);
		    my $masked_seq = $chr->mask_snps_and_repeats(begin_position       => $tr_start, 
								 end_position         => $tr_stop,
								 sequence_base_string => $t_region->nucleotide_seq);
		    @masked_array = split//,$masked_seq;
		}

		for my $n ($tr_start..$tr_stop) {

		    my ($refbase,$masked_refbase);
		    if ($t_region->{structure_type} =~ /cds/) {
			if ($strand eq "-1") {$base_n--;} else {$base_n++;}
			$refbase = $nucleotide_array[$base_n];
			if ($self->masked) {$masked_refbase = $masked_array[$base_n];}
		    } elsif ($self->utr_seq) {
			if ($strand eq "-1") {$base_n--;} else {$base_n++;}
			$refbase = $nucleotide_array[$base_n];
			if ($self->masked) {$masked_refbase = $masked_array[$base_n];}
		    } else {
			$refbase = "$n:$strand";
			if ($self->masked) {$masked_refbase = "$n:$strand";}
		    }
		    $info->{$transcript}->{ref_base}->{$n}="$trv_type,$refbase";
		    if ($self->masked) {$info->{$transcript}->{masked_ref_base}->{$n}="$trv_type,$masked_refbase";}
		    $info->{$transcript}->{range}->{$n}="$tr_start-$tr_stop";
		    if ($strand eq "-1") {
			$info->{$transcript}->{range}->{$n}="$tr_stop-$tr_start";
		    }
		}
	    }
	}
    } else {
	my $line = qq(\nCould not find substructures in the transcript object for $transcript from the $organism data warehouse\nWill exit the program now\n\n);
	if ($self->output) {print OUT qq($line);}
	unless ($self->no_stdout) {print qq($line);}
	exit 1;	
    }
    
    my $exon=0;
    my $previous_coord;
    my $frame=0;
    my $aa_n=1;

    my @positions;
    if ($info->{$transcript}->{strand} eq "-1") {
	foreach my $gcoord (sort {$b<=>$a} keys %{$info->{$transcript}->{ref_base}}) {
	    push(@positions,$gcoord);
	}
    } else {
	foreach my $gcoord (sort {$a<=>$b} keys %{$info->{$transcript}->{ref_base}}) {
	    push(@positions,$gcoord);
	}
    }

    for my $gcoord (@positions) {
	my ($region,$base) = split(/\,/,$info->{$transcript}->{ref_base}->{$gcoord});
	my ($range) = $info->{$transcript}->{range}->{$gcoord};
	my ($maskedregion,$masked_base);
	if ($self->masked) {($maskedregion,$masked_base) = split(/\,/,$info->{$transcript}->{masked_ref_base}->{$gcoord});}

	unless ($previous_coord) {$previous_coord = $gcoord;}

	if ($info->{$transcript}->{strand} eq "-1") {
	    unless ($gcoord + 1 == $previous_coord) {
		$exon++;
	    }
	} else {
	    unless ($gcoord - 1 == $previous_coord) {
		$exon++;
	    }
	}

	if ($region =~ /utr/) {
	    $frame = "-";
	} else {
	    $frame++;
	}
	
	if ($trans_pos == $gcoord) {
	    $transcript_info->{$transcript}->{$gcoord}->{trans_pos}="$trans_pos,$region";
	    #print qq(trans_pos $trans_pos = $gcoord\n);
	} else {
	    if ($info->{$transcript}->{strand} eq "-1") {
		if ($trans_pos < $previous_coord && $trans_pos > $gcoord) {
		    $transcript_info->{$transcript}->{$previous_coord}->{trans_pos}="$trans_pos,intron";
		}
	    } else {
		if ($trans_pos > $previous_coord && $trans_pos < $gcoord) {
		    $transcript_info->{$transcript}->{$previous_coord}->{trans_pos}="$trans_pos,intron";
		}
	    }
	}
	
	$previous_coord = $gcoord;
	$transcript_info->{$transcript}->{$gcoord}->{exon}="$exon,$region";
	$transcript_info->{$transcript}->{$gcoord}->{frame}=$frame;
	$transcript_info->{$transcript}->{$gcoord}->{aa_n}=$aa_n;
	$transcript_info->{$transcript}->{$gcoord}->{base}=$base;
	if ($self->masked) {$transcript_info->{$transcript}->{$gcoord}->{masked_base}=$masked_base;}
	$transcript_info->{$transcript}->{$gcoord}->{range}=$range;
	$transcript_info->{$transcript}->{-1}->{exon_total}=$exon;
	if ($frame eq "3") {$frame=0;$aa_n++;}
    }
    return($transcript_info,$chromosome,$strand);
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



=head1 TITLE

TranscriptInformation

=head1 DESCRIPTION

This script will get transcript information

=head1 Input Options:

transcript
trans_pos
organism

=head1 KNOWN BUGS

Please report bugs to <rmeyer@genome.wustl.edu>

=head1 AUTHOR

Rick Meyer <rmeyer@genome.wustl.edu>

=cut
