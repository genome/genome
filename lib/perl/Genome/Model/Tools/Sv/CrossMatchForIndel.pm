package Genome::Model::Tools::Sv::CrossMatchForIndel;

use strict;
use warnings;
use Genome;
use Bio::SeqIO;

=cut
my %opts = (m=>0.02,s=>0,h=>100,u=>30,b=>0,z=>0);
getopts('x:s:m:h:iu:g:r:c:db:z:',\%opts);
die("
Usage:   getCrossMatchIndel.pl <crossmatch alignment result>\n
Output:  breakpoint1 (chr/pos), breakpoint2 (chr/pos), contig breakpoint, Size, Type, Contig, AlnScore
Options:
         -x STRING  start position of the reference [chr_pos]
         -b INT     reference concatenation position [$opts{b}]
         -z INT     Add variant size by [$opts{z}] bp if breakpoints are on the opposite side of the concatenated reference
         -i         Only detect Indels
         -s INT     minimal size of the indels to report
         -m INT     percent substitution rate in the flanking alignment [$opts{m}]
         -h INT     Maximal microhomology size for Non-homologous end joining [$opts{h}] (bp)
         -u INT     minimal length of unique sequence in the alignment [$opts{u}] (bp)
         -r FILE    fasta file that contains the local reference sequence
         -c FILE    Assembly fasta file containing all the contigs
         -d         Output microhomology debugging information
") unless (@ARGV);
=cut

class Genome::Model::Tools::Sv::CrossMatchForIndel {
    is  => 'Genome::Model::Tools::Sv',
    has => [
        cross_match_file => {
            type => 'String',
            doc  => 'The input cross match out file to parse',
        },
        local_ref_seq_file => {
            type => 'String',
            doc  => 'fasta file that contains the local reference sequence',
        },
        assembly_contig_file => {
            type => 'String',
            doc  => 'Assembly fasta file containing all the contigs',
        },
    ],
    has_optional => [
        per_sub_rate  => {
            type => 'Number',
            doc  => 'percent substitution rate in the flanking alignment',
            default_value => 0.02,
        },
        min_indel_size => {
            type => 'Number',
            doc  => 'minimal size of the indels to report',
            default_value => 0,
        },
        min_unique_length => {
            type => 'Number',
            doc  => 'Minimal length(bp) of unique sequence in the alignment',
            default_value => 30,
        },
        ref_start_pos => {
            type => 'String',
            doc  => 'start position of the reference [chr_pos]',
        },
        ref_cat_pos  => {
            type => 'Number',
            doc  => 'reference concatenation position',
            default_value => 0,
        },
        variant_size => {
            type => 'Number',
            doc  => 'Add variant size by bp if breakpoints are on the opposite side of the concatenated reference',
            default_value => 0,
        },
        indel_only   => {
            type => 'Boolean',
            doc  => 'Only detect Indels',
            default_value => 0,
        },
        max_microhomology_size => {
            type => 'Number',
            doc  => 'Maximal microhomology size(bp) for Non-homologous end joining',
            default_value => 100,
        },
        microhomology_log_file => {
            type => 'String',
            doc  => 'Output microhomology debugging information',
        },
        _refseq => {
            type => 'SCALAR',
        },
    ],
};
            

#TODO get rid of out_fh(output_file) and directly return an output
#scalar string

sub execute {
    my $self = shift;
    my $cm = Genome::Model::Tools::Sv::ParseCrossMatch->create(
        input_file => $self->cross_match_file
    );
    
    my @DCposes = keys %{$cm->dcpos}; #discrepant position
    my $swap_chrom = 0;

    my ($refseq, @vars, @refbases);

    if ($self->local_ref_seq_file){
        my $local_ref_seq = $self->local_ref_seq_file;
        unless (-s $local_ref_seq) {
            $self->error_message('local ref seq fasta file '. $local_ref_seq .' is invalid');
            return;
        }
        
        my $refseq_stream = Bio::SeqIO->newFh(-file => $local_ref_seq, -format => 'Fasta');
        my $refseq_obj = <$refseq_stream>;
        $refseq   = uc($refseq_obj->seq);
        @refbases = split //, $refseq;
    }

    my $log_fh;
    if ($self->microhomology_log_file) {
        $log_fh = Genome::Sys->open_file_for_writing($self->microhomology_log_file) or return;
    }

    my $len_refseq = length($refseq);
    $self->_refseq($refseq);

    my ($chr1,$refpos1,$chr2,$refpos2,$pretype,$presize,$preori);
    if ($self->ref_start_pos) {
        ($chr1,$refpos1,$chr2,$refpos2,$pretype,$presize,$preori) = split /\_/, $self->ref_start_pos;
    }

#examining gapped indels
    for my $refpos (@DCposes){
        next if !defined $cm->dcpos->{$refpos};
        #my $refbase=substr $refseq, $refpos-1, 1;
        my @dcs = @{$cm->dcpos->{$refpos}};
        for my $dc (@dcs){
            my ($type, $size) = $dc->{type} =~ /([DI])\-*(\d*)/;
            if (defined $size && $size =~ /\S+/) {
                my ($read_len,$trimmed_readlen,$nbp_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$strand) = $self->_ReadStats($dc->{read}, $cm);
                my $var;
                my $altseq     = $self->_GetContig($dc->{read});
                my @altbases   = split //, $altseq;
                my $len_altseq = length($altseq);
                my $bkpos      = $dc->{rpos};

                if ($type eq 'D') {
	                $type = 'DEL';
	                my ($idx1, $idx2) = (0, 0);

            	    #Search for microhomology
	                if ($strand eq '+') {
                        #$bkpos++;
                        #while ($refbases[$refpos-$idx1-1-1] eq $refbases[$refpos+$size-$idx1-1-1] && $refpos-$idx1-1-1>=0){
                        while ($refbases[$refpos-$idx1-1-1] eq $refbases[$refpos+$size-$idx1-1-1] &&
		                    $refbases[$refpos-$idx1-1-1] eq $altbases[$bkpos-$idx1-1] && 
                            $bkpos-1>=0 && $refpos-$idx1-1-1>=0
                        ){
	                        $idx1++;
	                    }
	                    my $P5_homology = $idx1;
	                    if ($refbases[$refpos+$idx2-1] eq $altbases[$bkpos+$idx2+1-1]) {
	                        while ($refbases[$refpos+$idx2-1] eq $refbases[$refpos+$idx2+$size-1] && 		  
                                $refbases[$refpos+$idx2-1] eq $altbases[$bkpos+$idx2+1-1] &&
		                        $bkpos+$idx2-1<$len_altseq && $refpos+$idx2+$size-1<$len_refseq
                            ){
	                            $idx2++;
	                        }
                        }
                        my $P3_homology = $idx2;
                        $var->{refpos1}       = $refpos-$P5_homology;
	                    $var->{refpos2}       = $refpos+$size+$P3_homology-1;
	                    $var->{bkpos1}        = $bkpos-$P5_homology+1;
	                    $var->{bkpos2}        = ($P5_homology+$P3_homology==0) ? '-' : $bkpos+$P3_homology;
                        $var->{microhomology} = $P5_homology+$P3_homology;
                    }
	                else {
                        #while ($refbases[$refpos-$idx1-1] eq $refbases[$refpos-$size-$idx1-1] && $refpos-$size-$idx1-1>=0){
	                    my $refseq2 = $refseq; 
                        $refseq2 =~ tr/ACGT/TGCA/;
	                    my @refbases2 = split //,$refseq2;
	                    while ($refbases2[$refpos-$idx1-1] eq $refbases2[$refpos-$size-$idx1-1] &&
		                    $refbases2[$refpos-$idx1-1] eq $altbases[$bkpos+$idx1+1-1] &&
		                    $bkpos+$idx1-1<$len_altseq && $refpos-$size-$idx1-1>=0
                        ){
	                        $idx1++;
                        }
	                    my $P3_homology = $idx1;
	                    if ($refbases2[$refpos] eq $altbases[$bkpos-1]) {
	                        while ($refbases2[$refpos+$idx2] eq $refbases2[$refpos-$size+$idx2] && 
                                $refbases2[$refpos+$idx2] eq $altbases[$bkpos-$idx2-1] &&
		                        $bkpos-$idx2+1-1>=0 && $refpos+$idx2+1<$len_refseq
                            ){
	                            $idx2++;
                            }
                        }
	                    my $P5_homology = $idx2;
	                    $var->{refpos1}       = $refpos-$size-$P3_homology+1;
	                    $var->{refpos2}       = $refpos+$P5_homology;
	                    $var->{bkpos1}        = $bkpos-$P5_homology+1;
	                    $var->{bkpos2}        = ($P5_homology+$P3_homology==0) ? '-' : $bkpos+$P3_homology;
                        $var->{microhomology} = $P5_homology+$P3_homology;
                        # change according to Ken's requirement:
#                         print STDERR "original: $var->{bkpos2}\n";
#                        $var->{refpos2} -= 1;
#                        $var->{bkpos2} -= 1;
#                        print STDERR "Now: $var->{bkpos2}\n";
                    }
                }
                else {
	                $type = 'INS';
	                my ($idx1, $idx2) = (0, 0);
	                #Search for microhomology
	                if ($strand eq '+') {
	                    while ($altbases[$bkpos-$idx1-1-1] eq $altbases[$bkpos+$size-$idx1-1-1] &&
		                    $refbases[$refpos-$idx1-1] eq $altbases[$bkpos-$idx1-1-1] &&
		                    $bkpos-$idx1-1-1>=0 && $refpos-$idx1-1>=0
                        ){
	                        $idx1++;
                        }
	                    my $P5_homology = $idx1;
	                    if ($refbases[$refpos+$idx2] eq $altbases[$bkpos+$idx2-1]) {
	                        while ($altbases[$bkpos+$idx2-1] eq $altbases[$bkpos+$idx2+$size-1] &&
		                        $refbases[$refpos+$idx2] eq $altbases[$bkpos+$idx2-1] &&
		                        $bkpos+$idx2+$size-1<$len_altseq && $refpos+$idx2<$len_refseq
                            ){
	                            $idx2++;
                            }
                        }
	                    my $P3_homology = $idx2;
	                    $var->{refpos1}       = $refpos-$P5_homology+1;
	                    $var->{refpos2}       = $refpos+$P3_homology+1;
	                    $var->{microhomology} = $P5_homology+$P3_homology;
	                    $var->{bkpos1}        = $bkpos-$var->{microhomology};
	                    $var->{bkpos2}        = $bkpos+$size-1;
                    }
	                else {
	                    my $refseq2 = $refseq; 
                        $refseq2 =~ tr/ACGT/TGCA/;
	                    my @refbases2 = split //,$refseq2;
	                    while ($altbases[$bkpos+$idx1-1] eq $altbases[$bkpos+$idx1+$size-1] &&
		                    $refbases2[$refpos-$idx1-1-1] eq $altbases[$bkpos+$idx1-1] &&
		                    $bkpos+$idx1-1<$len_altseq && $refpos-$idx1-1-1>=0
                        ){
	                        $idx1++;
                        }
	                    my $P3_homology = $idx1;
	                    if ($refbases2[$refpos-1] eq $altbases[$bkpos-1-1]) {
	                        while ($altbases[$bkpos-$idx2-1-1] eq $altbases[$bkpos-$idx2+$size-1-1] &&
		                        $refbases2[$refpos+$idx2-1] eq $altbases[$bkpos-$idx2-1-1] &&
		                        $bkpos-$idx2-1-1>=0 && $refpos+$idx2-1<$len_refseq
                            ){
	                            $idx2++;
                            }
                        }
	                    my $P5_homology = $idx2;
	                    $var->{refpos1}       = $refpos-$P3_homology;
	                    $var->{refpos2}       = $refpos+$P5_homology;
	                    $var->{bkpos1}        = $bkpos;
	                    $var->{microhomology} = $P5_homology+$P3_homology;
	                    $var->{bkpos2}        = $bkpos+$size+$var->{microhomology}-1;
                    }
	            }
                ($var->{type},$var->{size},$var->{read},$var->{score}) = ($type,$size,$dc->{read},$dc->{score});
                $var->{orientation} = '+-';
                $var->{scar}=0;

                $nbp_indel-=$size;
                $n_indel--;
                my $fraction_aligned = ($nbp_aligned||0)/($trimmed_readlen||1);
                ($var->{read_len},$var->{fraction_aligned},$var->{n_seg},$var->{n_sub},$var->{n_indel},$var->{nbp_indel},$var->{strand}) = ($trimmed_readlen,$fraction_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$strand);
                
                my @alnstrs;
                for my $aln (@{$cm->_align->{$dc->{read}}->{aln}}){
	                push @alnstrs, $aln->{line};
                }
                $var->{alnstrs} = join ',', @alnstrs;
                push @vars, $var;
            }
        }
    }

#examining split-reads alignment for intra-chromosomal rearrangement
    for my $read (keys %{$cm->_align}) {
        my @alns = @{$cm->_align->{$read}->{aln}};
        my @alnstrs;
        my ($read_len,$trimmed_readlen,$nbp_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$strand) = $self->_ReadStats($read, $cm);
        
        for (my $i=0;$i<$#alns;$i++) {
            my $aln1 = $alns[$i];
            
            for (my $j=$i+1;$j<=$#alns;$j++) {
                my $aln2 = $alns[$j];
                my ($perc_mismatch1,$perc_mismatch2) = (0,0);
                for my $mistype ('perc_sub','perc_del','perc_ins') {
	                $perc_mismatch1 += $aln1->{$mistype};
	                $perc_mismatch2 += $aln2->{$mistype};
                }
                next unless $perc_mismatch1 < $self->per_sub_rate * 100 and $perc_mismatch2 < $self->per_sub_rate * 100;
                
                my ($refpos1, $refpos2, $type, $size, $orientation);
                my $score = $aln1->{score} + $aln2->{score};
                
                my $uniq_query1 = $aln2->{r_start} - $aln1->{r_start};
                my $uniq_query2 = $aln2->{r_end}   - $aln1->{r_end};
                
                next unless $uniq_query1 >= $self->min_unique_length and $uniq_query2 >= $self->min_unique_length;
                
                my $sbj_overlap  = $aln1->{r_end} - $aln2->{r_start} + 1;
                my $sbj_homology = $sbj_overlap>0 ? $sbj_overlap : 0;
                my $sbj_gain     = $sbj_overlap<0 ? -$sbj_overlap : 0;
                my ($bkpos1,$bkpos2) = ($aln2->{r_start}, $aln1->{r_end});

                my @tmpalns = ($aln1,$aln2);
                ($nbp_aligned,$n_sub,$n_indel,$nbp_indel,$strand) = &AlignStats(@tmpalns);
                my $fraction_aligned = ($nbp_aligned||0)/($trimmed_readlen||1);
                
                if ($aln1->{refseq} eq $aln2->{refseq}){  # intra-chromosomal
	                my $readspan = $aln2->{r_start} - $aln1->{r_start};
	                my $refspan  = $aln2->{ref_start} - $aln1->{ref_start};
	                $size = abs($refspan) - $readspan;
	                if ($aln1->{orientation} eq $aln2->{orientation}) {   # indel, tandem dup
	                    if ($sbj_overlap > 0) {
	                        ($aln1, $aln2) = $self->FixMicrohomology($read,$sbj_overlap,$aln1,$aln2);
	                        $sbj_overlap   = $aln1->{r_end} - $aln2->{r_start} + 1;
                            $sbj_homology  = ($sbj_overlap>0) ? $sbj_overlap : 0;
	                        $sbj_gain      = ($sbj_overlap<0) ? -$sbj_overlap : 0;

                        }
	                    $readspan = $aln2->{r_start}-$aln1->{r_start};
	                    $refspan  = $aln2->{ref_start}-$aln1->{ref_start};
                        my $ref_overlap = ($aln1->{orientation} eq 'U')?$aln1->{ref_end}-$aln2->{ref_start}+1:$aln2->{ref_start}-$aln1->{ref_end}+1;

	                    if ($sbj_gain>0) {
	                        $bkpos1 = $aln1->{r_end}+1;
	                        $bkpos2 = $aln2->{r_start}-1;
                        }
	                    if ($sbj_overlap==0) {  # no microhomology
	                        $bkpos1 = $aln1->{r_end}+1;
	                        $bkpos2 = '-';
                        }
	                    $size = abs($refspan) - $readspan;
	                    $orientation = '+-';
                        
	                    if ($size>0 || $ref_overlap<0) {  # Deletion, loss of bases in the reference
	                        next if $aln1->{orientation} eq 'U' && $refspan < $self->min_unique_length ||   # ignore cross-mapping
		                        $aln1->{orientation} eq 'C' && $refspan > -($self->min_unique_length);
	                        $type = 'DEL';
	                        if ($sbj_gain == 0) {
	                            if ($aln1->{orientation} eq 'U'){
		                            $refpos1 = $aln1->{ref_end}   - $sbj_overlap + 1;
		                            $refpos2 = $aln2->{ref_start} + $sbj_overlap - 1;
                                }
	                            else {
		                            $refpos1 = $aln2->{ref_start} - $sbj_overlap + 1;
		                            $refpos2 = $aln1->{ref_end}   + $sbj_overlap - 1;
                                }
                            }
	                        else {  #non-template insertion
	                            if ($aln1->{orientation} eq 'U') {
		                            $refpos1 = $aln1->{ref_end}   + 1;
		                            $refpos2 = $aln2->{ref_start} - 1;
                                }
	                            else {
		                            $refpos1 = $aln2->{ref_start} + 1;
		                            $refpos2 = $aln1->{ref_end}   - 1;
                                }
                            }
                        }
	                    else {  # Insertion/Duplication
	                        if($aln1->{orientation} eq 'U') {
		                        $refpos1 = $aln2->{ref_start};
		                        $refpos2 = $aln1->{ref_end};
                            }
	                        else {
		                        $refpos1 = $aln1->{ref_end};
		                        $refpos2 = $aln2->{ref_start};
                            }
                            if ($ref_overlap >= -$size-$sbj_gain) {  # tandem duplication
	                            $type = 'ITX';
                            }
                            else {
                                $type = 'INS';   # insertion
                            }
                        }
                    }
	                else {  #Inversion
	                    $type = 'INV';
	                    if ($aln1->{orientation} eq 'U') {
	                        $orientation = '++';
	                        $strand = '+';
	                        if ($aln1->{ref_end} < $aln2->{ref_start}){   # 1->|    <-2|
	                            $refpos1 =  $aln1->{ref_end} - $sbj_overlap + 1;
	                            $refpos2 =  $aln2->{ref_start};
	                            $refpos2 += $aln2->{r_rest} if $aln2->{ref_rest} <= 0;
	                            $size    =  $refpos2 - $refpos1 - $sbj_overlap + 1;
                            }
	                        else{   # |<-2    1->|
	                            $refpos1 =  $aln2->{ref_end};
	                            $refpos1 -= $aln2->{r_rest} if $aln2->{ref_rest} <= 0;
	                            $refpos2 =  $aln1->{ref_end};
	                            $refpos2 += $aln1->{r_rest} if $aln1->{ref_rest} <= 0;
	                            $size    =  $refpos2 - $refpos1 + 1;
                            }
                        }
	                    else{
	                        $orientation = '--';
	                        $strand      = '-';
	                        if ($aln1->{ref_end} > $aln2->{ref_start}) {   # |2->    |<-1
	                            $refpos1 =  $aln2->{ref_start};
	                            $refpos1 -= $aln2->{r_rest} if $aln2->{ref_rest} <= 0;
	                            $refpos2 =  $aln1->{ref_end} + $sbj_overlap - 1;
	                            $size    =  $refpos2 - $refpos1 - $sbj_overlap + 1;
                            }
	                        else {   # |<-1    2->|
	                            $refpos1 =  $aln1->{ref_end};
	                            $refpos1 -= $aln1->{r_rest} if $aln1->{ref_rest} <= 0;
	                            $refpos2 =  $aln2->{ref_end};
	                            $refpos2 += $aln2->{r_rest} if $aln2->{ref_rest} <= 0;
	                            $size    =  $refpos2 - $refpos1 + 1;
                            }
                        }

                    }
                }
                else {  #Inter-chromosomal
	                $size = 1;
	                $type = 'CTX';
                        #my ($chrom1) = $aln1->{refseq} =~ /\.([\w\d]+)\.fa*/;
                        #my ($chrom1)=($aln1->{refseq}=~/[\.\/]([\w\d]+)\.fa*/); # to handle mouse reference as well
                        my ($chrom1) = ($aln1->{refseq} =~ /(\S+):\d+-\d+/);
                    #if (&GLess($aln1->{refseq}, $aln2->{refseq})) { #keep the repeat in the lower chromosome
                    if (!defined $chr1 || $chrom1 eq $chr1) {    
	                    if ($aln1->{orientation} eq 'U'){
	                        if ($aln2->{orientation} eq 'U') {
	                            $orientation = '+-';
                                $refpos1 = $aln1->{ref_end}-$sbj_overlap+1;
	                            $refpos2 = $aln2->{ref_start}+$sbj_overlap-1;
                            }
	                        else{
	                            $orientation = '++';
                                $refpos1 = $aln1->{ref_end}-$sbj_overlap+1;
	                            $refpos2 = $aln2->{ref_start}-$sbj_overlap+1;
                            }
                        }
	                    else{
	                        if ($aln2->{orientation} eq 'U') {
	                            $orientation='--';
                                $refpos1 = $aln1->{ref_end}+$sbj_overlap-1;
	                            $refpos2 = $aln2->{ref_start}+$sbj_overlap-1;
                            }
	                        else{
	                            $orientation='-+';
  	                            $refpos1 = $aln1->{ref_end}+$sbj_overlap-1;
	                            $refpos2 = $aln2->{ref_start}-$sbj_overlap+1;
                            }
                        }
                    }
	                else {
	                    if ($aln1->{orientation} eq 'U') {
	                        if ($aln2->{orientation} eq 'U') {
	                            $orientation = '-+';
                                $refpos1 = $aln1->{ref_end}-$sbj_overlap+1;
	                            $refpos2 = $aln2->{ref_start}+$sbj_overlap-1;
                            }
	                        else {
	                            $orientation = '++';
                                $refpos1 = $aln1->{ref_end}-$sbj_overlap+1;
	                            $refpos2 = $aln2->{ref_start}-$sbj_overlap+1;
                            }
                        }
	                    else {
	                        if ($aln2->{orientation} eq 'U') {
	                            $orientation = '--';
                                $refpos1 = $aln1->{ref_end}+$sbj_overlap-1;
	                            $refpos2 = $aln2->{ref_start}+$sbj_overlap-1;
                            }
	                        else {
	                            $orientation = '+-';
                                $refpos1 = $aln1->{ref_end}+$sbj_overlap-1;
	                            $refpos2 = $aln2->{ref_start}-$sbj_overlap+1;
                            }
                        }
	                    my $tmp  = $refpos1;
                        $refpos1 = $refpos2;
                        $refpos2 = $tmp;
	                    $tmp     = $aln1;
                        $aln1    = $aln2;
                        $aln2    = $tmp;
                    }
                }
                if ($bkpos2 ne '-' && $bkpos1 > $bkpos2) {
	                my $tmp = $bkpos1;
                    $bkpos1 = $bkpos2;
                    $bkpos2 = $tmp;
                }
                my $var;
                ($var->{type},$var->{size},$var->{chr1},$var->{refpos1},$var->{chr2},$var->{refpos2},$var->{orientation},$var->{bkpos1},$var->{bkpos2},$var->{read},$var->{score},$var->{scar},$var->{read_len},$var->{fraction_aligned},$var->{n_seg},$var->{n_sub},$var->{n_indel},$var->{nbp_indel},$var->{strand},$var->{microhomology}) = ($type,abs($size),$aln1->{refseq},$refpos1,$aln2->{refseq},$refpos2,$orientation,$bkpos1,$bkpos2,$read,$score,$sbj_gain,$trimmed_readlen,$fraction_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$strand,$sbj_homology);
                $var->{alnstrs} = join ',', $aln1->{line}, $aln2->{line};
                push @vars, $var if $score > 0;
            }
        }
    }

    #Find best answer (highest score)
    my $maxscore = 0;
    my $mindist  = 1e10;
    my ($seq, $bestvar, $out_str);
    
    for my $var (@vars) {
        next if $var->{type} eq 'INV' && $self->indel_only;
        my $dist;
        my $sametype = 0;
        if (defined $pretype && ($var->{type} eq $pretype ||
		    $var->{type} eq 'INS' && $pretype eq 'ITX' ||
		    $var->{type} eq 'ITX' && $pretype eq 'INS')) {
            $sametype = 1;
        }
        if (defined $pretype && defined $presize) {
            $dist = $sametype ? abs($presize - $var->{size}) : 1e10;
            $var->{size} = $presize if $pretype eq 'CTX';
        }
        if (defined $preori && $pretype eq 'CTX' && $preori ne $var->{orientation}) {
            $var->{score} = 0;  #CTX must confirmed in the predicted orientation
        }
        next if $var->{score} <= 0;
        
        if (!defined $bestvar ||
            (defined $pretype && $sametype && ($mindist>$dist #same type and larger in size
            || $mindist==$dist && $bestvar->{score}<$var->{score})) ||
            (!defined $pretype && ($bestvar->{score}<=$var->{score}||$bestvar->{n_sub}>$var->{n_sub}||$bestvar->{n_indel}>$var->{n_indel}||$bestvar->{nbp_indel}>$var->{nbp_indel}||$bestvar->{size}<$var->{size}))) {   #same score but larger size
            $bestvar  = $var;
            $maxscore = $var->{score};
            $mindist  = $dist;
        }
    }

    if (defined $bestvar) {
  #Print out
        if ($bestvar->{type} =~ /INS/) {
#            $refpos1++;
        }

        #Microhomology Standardization
        if ($bestvar->{size} >= $self->min_indel_size && $bestvar->{score} > 0) {
            if (defined $self->ref_start_pos) {
                if (defined $bestvar->{chr1} && ($bestvar->{chr1}=~/^(\S+)\:/ || $bestvar->{chr1}=~/chromosome\.(\S+)\./)) { 
                    $bestvar->{chr1} = $1;
                }
                      
                if (defined $bestvar->{chr2} && ($bestvar->{chr2}=~/^(\S+)\:/ || $bestvar->{chr2}=~/chromosome\.(\S+)\./)) {
	                $bestvar->{chr2}=$1;
                }
                my ($pos1, $pos2);
                if ($self->ref_cat_pos > 0){
	                $pos1 = (($bestvar->{refpos1} <= $self->ref_cat_pos) ? $refpos1 : $refpos2) + $bestvar->{refpos1}-1;
	                $pos2 = (($bestvar->{refpos2} <= $self->ref_cat_pos) ? $refpos1 : $refpos2) + $bestvar->{refpos2}-1;
	                $bestvar->{size} += $self->variant_size if ($bestvar->{refpos1} - $self->ref_cat_pos)*($bestvar->{refpos2} - $self->ref_cat_pos) < 0;
                }
                else {
	                $pos1 = $refpos1 + $bestvar->{refpos1} - 1;
	                $pos2 = ($refpos2||$refpos1) + ($bestvar->{refpos2}||$bestvar->{refpos1}) - 1;
                }
                $out_str = sprintf("%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%s\t%d\t%s",$bestvar->{chr1}||$chr1,$pos1,$bestvar->{chr2}||$chr1,$pos2,$bestvar->{orientation}||'+-',$bestvar->{bkpos1},$bestvar->{bkpos2}||$bestvar->{bkpos1},$bestvar->{size},$bestvar->{type},$bestvar->{read},$bestvar->{score},$bestvar->{scar}||0,$bestvar->{read_len},$bestvar->{fraction_aligned},$bestvar->{n_seg},$bestvar->{n_sub},$bestvar->{n_indel},$bestvar->{nbp_indel},$bestvar->{strand},$bestvar->{microhomology},$bestvar->{alnstrs});
            }
            else {
                $out_str = sprintf("%d\t%d\t%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%s\t%d",$bestvar->{refpos1},$bestvar->{refpos2},$bestvar->{orientation}||'+-',$bestvar->{bkpos1},$bestvar->{bkpos2}||$bestvar->{bkpos1},$bestvar->{size},$bestvar->{type},$bestvar->{read},$bestvar->{score},$bestvar->{scar} || 0,$bestvar->{read_len},$bestvar->{fraction_aligned},$bestvar->{n_seg},$bestvar->{n_sub},$bestvar->{n_indel},$bestvar->{nbp_indel},$bestvar->{strand},$bestvar->{microhomology});
            }
            if (defined $self->ref_start_pos){
                $out_str .= sprintf("\t%s\n", $self->ref_start_pos);
            }
            else{
                $out_str .= "\n";
            }

            my $altseq   = $self->_GetContig($bestvar->{read});
            my @altbases = split //, $altseq;
            my @tmp      = split /\_/, $self->ref_start_pos;
            
            if ($bestvar->{strand} eq '+') {
                my $altbase1 = $altbases[$bestvar->{bkpos1}-2];
                my $altbase2 = $bestvar->{bkpos2} eq '-' ? '-' : $altbases[$bestvar->{bkpos2}];

                if (($refbases[$bestvar->{refpos1}-2] ne $altbase1 || $altbase2 ne '-' && $refbases[$bestvar->{refpos2}] ne $altbase2) && $log_fh) {
                    # Microhomology check
	                $log_fh->print("$tmp[1]\t") if defined $self->ref_start_pos;
	                $log_fh->printf("%d\t%d\t%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%s\t%d\t",$bestvar->{refpos1},$bestvar->{refpos2},$bestvar->{orientation}||'+-',$bestvar->{bkpos1},$bestvar->{bkpos2}||$bestvar->{bkpos1},$bestvar->{size},$bestvar->{type},$bestvar->{read},$bestvar->{score},$bestvar->{scar} || 0,$bestvar->{read_len},$bestvar->{fraction_aligned},$bestvar->{n_seg},$bestvar->{n_sub},$bestvar->{n_indel},$bestvar->{nbp_indel},$bestvar->{strand},$bestvar->{microhomology});
	                $log_fh->print("Discrepancy: $refbases[$bestvar->{refpos1}-2] vs $altbase1 \/ $refbases[$bestvar->{refpos2}] vs $altbase2\n");
                }
            }
            else {
                my $altbase1 = ($bestvar->{bkpos2} eq '-') ? '-' : $altbases[$bestvar->{bkpos2}+1-1];
                $altbase1 =~ tr/ACGT/TGCA/;
                my $altbase2 = $altbases[$bestvar->{bkpos1}-2];
                $altbase2 =~ tr/ACGT/TGCA/;

                if (($altbase1 ne '-' && $refbases[$bestvar->{refpos1}-2] ne $altbase1 || $refbases[$bestvar->{refpos2}] ne $altbase2) && $log_fh){
                    #       Microhomology check
	                $log_fh->print("$tmp[1]\t") if defined $self->ref_start_pos;
	                $log_fh->printf("%d\t%d\t%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%s\t%d\t",$bestvar->{refpos1},$bestvar->{refpos2},$bestvar->{orientation}||'+-',$bestvar->{bkpos1},$bestvar->{bkpos2}||$bestvar->{bkpos1},$bestvar->{size},$bestvar->{type},$bestvar->{read},$bestvar->{score},$bestvar->{scar} || 0,$bestvar->{read_len},$bestvar->{fraction_aligned},$bestvar->{n_seg},$bestvar->{n_sub},$bestvar->{n_indel},$bestvar->{nbp_indel},$bestvar->{strand},$bestvar->{microhomology});
	                $log_fh->print("Discrepancy: $refbases[$bestvar->{refpos1}-2] vs $altbase1 \/ $refbases[$bestvar->{refpos2}] vs $altbase2\n");
                }
            }
        }
    }
    $log_fh->close if $log_fh;

    # unless explicitly deleted, this rather large object will stay in the UR
    # cache forever
    $cm->delete;
    $cm = undef;

    return $out_str;
}


sub GLess{
    my @chroms = @_;
    #resolve chromosome naming ambiguity
    for (my $i=0;$i<=$#chroms;$i++) {
        my $chrom = $chroms[$i];
        if ($chrom =~ /\.([\w\d]+)\.fa*/){
            $chroms[$i] = $1;
            $chroms[$i] = 23 if $chroms[$i] eq 'X';
            $chroms[$i] = 24 if $chroms[$i] eq 'Y';
        }
    }
    my $less = $chroms[0] < $chroms[1] ? 1 : 0;
    return $less;
}


sub AlnStrand{
    my @alns = @_;
    my $strand;
    if ($alns[$#alns]->{ref_end} > $alns[0]->{ref_start}){
        $strand = '+';
    }
    else {
        $strand = '-';
    }
    return $strand;
}


sub _ReadStats {
    my ($self, $read, $cm) = @_;
    
    my @alns = @{$cm->_align->{$read}->{aln}};
    my $readlen = $alns[$#alns]->{r_end} + $alns[$#alns]->{r_rest} || 0;
    my $trimmed_readlen = $alns[$#alns]->{r_end} - $alns[0]->{r_start} + 1;
    my $strand = &AlnStrand(@alns);
    
    #Add Tail
    $trimmed_readlen += ($alns[$#alns]->{r_rest} < $alns[$#alns]->{ref_rest}) ? $alns[$#alns]->{r_rest} : $alns[$#alns]->{ref_rest};

    #Add Head
    if ($alns[0]->{orientation} eq 'U'){
        $trimmed_readlen += ($alns[0]->{r_start} < $alns[0]->{ref_start}) ? $alns[0]->{r_start}-1 : $alns[0]->{ref_start}-1;
    }
    else {
        $trimmed_readlen += ($alns[0]->{r_start} < $alns[0]->{ref_end}) ? $alns[0]->{r_start}-1 : $alns[0]->{ref_end}-1;
    }
    my ($nbp_aligned,$n_sub,$n_indel,$nbp_indel) = &AlignStats(@alns);
    return ($readlen,$trimmed_readlen,$nbp_aligned,$#alns+1,$n_sub,$n_indel,$nbp_indel,$strand);
}


sub AlignStats{
    my @alns = @_;
    my ($nbp_aligned,$n_sub,$n_indel) = (0,0,0);
    my $nbp_indel = 0;
    my $overlap   = 0;
    for (my $i=0;$i<=$#alns;$i++) {
        my $aln = $alns[$i];
        my $overlap = $i > 0 ? $alns[$i-1]->{r_end} - $aln->{r_start} + 1 : 0;
        my $alnlen  = $aln->{r_end} - $aln->{r_start} + 1;
        $nbp_aligned += $alnlen-$overlap;
        if (defined $aln->{discrepancy}) {
            my @mismatch = @{$aln->{discrepancy}};
            for my $disc (@mismatch)    {
	            my ($type, $size) = $disc=~/([DI])\-*(\d*)/;
	            if (defined $type) {
	                $n_indel++;
	                $nbp_indel += (defined $size && $size=~/^\d+/) ? $size : 1;
                }
	            else{
	                $n_sub++;
                }
            }
        }
    }
    my $strand = &AlnStrand(@alns);
    return ($nbp_aligned,$n_sub,$n_indel,$nbp_indel,$strand);
}


sub _GetContig{
    my ($self, $contigid) = @_;
    my $fin = $self->assembly_contig_file;
    my $in  = Bio::SeqIO->newFh(-file => $fin , -format => 'Fasta');
    my $sequence;
    while ( my $seq = <$in> ) {
    # do something with $seq
        next unless $seq->id eq $contigid;
        $sequence = uc($seq->seq());
        last;
    }
    return $sequence;
}

sub FixMicrohomology{
    my ($self, $read, $overlap, $aln1, $aln2) = @_;
    my $altseq = $self->_GetContig($read);
    my $refseq = $self->_refseq;
  # fix aln1
    my $idx    = 0;
    my $offset = 0;
    my ($refpos, $rpos);
    #Adjust breakpoint position, make sure no mismatch in the first flanking base
    my ($ref_flankbase, $r_flankbase);
    do{
        $r_flankbase = substr($altseq,$aln1->{r_end}-$overlap-1,1);
        if ($aln1->{orientation} eq 'U') {
            $ref_flankbase = substr($refseq,$aln1->{ref_end}-$overlap-1,1);
        }
        else {
            $ref_flankbase = substr($refseq,$aln1->{ref_end}+$overlap-1,1);
            $r_flankbase   =~ tr/ACGT/TGCA/;
        }
        if ($ref_flankbase ne $r_flankbase) {
            $overlap++;
            $offset++;
        }
    } until ($ref_flankbase eq $r_flankbase);

    $aln2->{r_start} -= $offset;
    if ($aln2->{orientation} eq 'U') {
        $aln2->{ref_start} -= $offset;
    }
    else{
        $aln2->{ref_start} += $offset;
    }

    $offset=0;
    do{
        $r_flankbase = substr($altseq,$aln2->{r_start}+$overlap-1,1);
        if ($aln2->{orientation} eq 'U') {
            $ref_flankbase=substr($refseq,$aln2->{ref_start}+$overlap-1,1);
        }
        else{
            $ref_flankbase =  substr($refseq,$aln2->{ref_start}-$overlap-1,1);
            $r_flankbase   =~ tr/ACGT/TGCA/;
        }
        if ($ref_flankbase ne $r_flankbase) {
            $overlap++;
            $offset++;
        }
    } until ($ref_flankbase eq $r_flankbase);
    
    $aln1->{r_end} += $offset;
    if ($aln1->{orientation} eq 'U') {
        $aln1->{ref_end} += $offset;
    }
    else{
        $aln1->{ref_end} -= $offset;
    }

  #Adjust Microhomology (don't allow mismatch in Microhomology)
    if ($aln1->{orientation} eq 'U'){
        $refpos = $aln1->{ref_end}-$overlap+1;
        $rpos   = $aln1->{r_end}-$overlap+1;
        my $refbase = substr($refseq,$refpos-1+$idx,1);
        my $rbase   = substr($altseq,$rpos-1+$idx,1);
        while ($refbase eq $rbase && $idx<$overlap){
            $idx++;
            $refbase = substr($refseq,$refpos-1+$idx,1);
            $rbase   = substr($altseq,$rpos-1+$idx,1);
        }
        $aln1->{ref_end} = $refpos+$idx-1;
        $aln1->{r_end}   = $rpos+$idx-1;
    }
    else{
        my $refpos = $aln1->{ref_end}+$overlap-1;
        my $rpos   = $aln1->{r_end}-$overlap+1;

        my $refbase= substr($refseq,$refpos-1-$idx,1);
        my $rbase  = substr($altseq,$rpos-1+$idx,1); 
        $rbase =~ tr/ACGT/TGCA/;
        
        while ($refbase eq $rbase && $idx<$overlap) {
            $idx++;
            $refbase = substr($refseq,$refpos-1-$idx,1);
            $rbase   = substr($altseq,$rpos-1+$idx,1); 
            $rbase =~ tr/ACGT/TGCA/;
        }
        $aln1->{ref_end} = $refpos-$idx+1;
        $aln1->{r_end}   = $rpos+$idx-1;
    }

    #fix aln2
    $idx = 0;
    if ($aln2->{orientation} eq 'U') {
        my $refpos = $aln2->{ref_start}+$overlap-1;
        my $rpos   = $aln2->{r_start}+$overlap-1;

        my $refbase= substr($refseq,$refpos-1-$idx,1);
        my $rbase  = substr($altseq,$rpos-1-$idx,1);
        while ($refbase eq $rbase && $idx < $overlap) {
            $idx++;
            $refbase = substr($refseq,$refpos-1-$idx,1);
            $rbase   = substr($altseq,$rpos-1-$idx,1);
        }
        $aln2->{ref_start} = $refpos-$idx+1;
        $aln2->{r_start}   = $rpos-$idx+1;
    }
    else {
        my $refpos = $aln2->{ref_start}-$overlap+1;
        my $rpos   = $aln2->{r_start}+$overlap-1;

        my $refbase= substr($refseq,$refpos-1+$idx,1);
        my $rbase  = substr($altseq,$rpos-1-$idx,1); 
        $rbase =~ tr/ACGT/TGCA/;
        
        while ($refbase eq $rbase && $idx < $overlap){
            $idx++;
            $refbase = substr($refseq,$refpos-1+$idx,1);
            $rbase   = substr($altseq,$rpos-1-$idx,1); 
            $rbase =~ tr/ACGT/TGCA/;
        }
        $aln2->{ref_start} = $refpos+$idx-1;
        $aln2->{r_start}   = $rpos-$idx+1;
    }
    return ($aln1,$aln2);
}

1;
