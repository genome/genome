package Genome::Model::Tools::Parse::Crossmatch;

use strict;
use warnings;

use Genome;
use Genome::Model::Tools::Crossmatch::Fisher;
use Bio::SeqIO;

class Genome::Model::Tools::Parse::Crossmatch {
    is => 'Command',
    has => [
      chr_pos =>
      {
          type => 'String',
          is_optional=>1,
          doc =>"chr_pos string written like chr_pos",
      },
      min_indel_size => 
      {
          type => 'Integer',
          is_optional=>1,
          default=>1,
          doc =>"minimum indel size to report"
      },
      sub_rate => 
      {
          type => 'Integer',
          is_optional=>1,
          default=>2,
          doc=>'sub rate for flanking alignment in percent',
      },
      max_homology_size_for_nhej =>
      {
          type => 'Integer',
          is_optional=>1,
          default=>100,
          doc=>'Non homologous end joining. FOR REALZ',
      },
      min_unique_sequence => 
      {
          type => 'Integer',
          is_optional=>1,
          default=>30,
          doc=>'seems pretty obvious holmes',
      },
      crossmatch_file => {
          type => 'String',
          is_optional=>0,
          doc=>'the output of crossmatch that you would like to parse.',
      },
      indels_only => {
          type => 'Boolean',
          default=>1,
          is_optional=>1,
          doc => 'indels only? default yes',
      },    
    ],
};

sub help_detail {
    return "
    Usage:   getCrossMatchIndel.pl <crossmatch alignment result>\n
    Output:  breakpoint1 (chr/pos), breakpoint2 (chr/pos), contig breakpoint, Size, Type, Contig, AlnScore
    Options:
    -x STRING  start position of the reference [chr_pos]
    -i         Only detect Indels
    -s INT     minimal size of the indels to report
    -m INT     percent substitution rate in the flanking alignment [2%]
    -h INT     Maximal microhomology size for Non-homologous end joining [100]
    -u INT     minimal length of unique sequence in the alignment [30]
    ";
}

sub execute{
    my $self = shift;
    my @vars;
    my $cm = Genome::Model::Tools::Crossmatch::Fisher->new(fin=>$self->crossmatch_file);
    my @DCposes=keys %{$cm->{dcpos}}; #discrepant position

    #examining gapped indels
    foreach my $refpos(@DCposes){
        next if(!defined $cm->{dcpos}{$refpos});
        #my $refbase=substr $refseq, $refpos-1, 1;
        my @dcs=@{$cm->{dcpos}{$refpos}};
        foreach my $dc(@dcs){
            my ($type,$size)=($dc->{type}=~/([DI])\-*(\d*)/);
            if(defined $size && $size=~/\S+/){
                my $var;
                if($type eq 'D'){
                    $type='DEL';
                    $var->{refpos1}=$refpos;
                    $var->{refpos2}=$refpos+$size;
                    $var->{bkpos1}=$dc->{rpos};
                    $var->{bkpos2}=$dc->{rpos};
                }
                else{
                    $type='INS';
                    $var->{refpos1}=$refpos;
                    $var->{refpos2}=$refpos;
                    $var->{bkpos1}=$dc->{rpos};
                    $var->{bkpos2}=$dc->{rpos}+$size;
                }
                ($var->{type},$var->{size},$var->{read},$var->{score})=($type,$size,$dc->{read},$dc->{score});
                $var->{orientation}='+-';
                push @vars,$var;
            }
        }
    }

    #examining split-reads alignment for intra-chromosomal rearrangement
    foreach my $read(keys %{$cm->{align}}){
        my @alns=@{$cm->{align}{$read}->{aln}};
        for(my $i=0;$i<$#alns;$i++){
            my $aln1=$alns[$i];
            for(my $j=$i+1;$j<=$#alns;$j++){
                my $aln2=$alns[$j];
                my ($perc_mismatch1,$perc_mismatch2)=(0,0);
                foreach my $mistype('perc_sub','perc_del','perc_ins'){
                    $perc_mismatch1+=$aln1->{$mistype};
                    $perc_mismatch2+=$aln2->{$mistype};
                }
                next unless($perc_mismatch1<$self->sub_rate && $perc_mismatch2<$self->sub_rate);
                my ($refpos1,$refpos2);
                my $type;
                my $score=$aln1->{score}+$aln2->{score};
                my $size;
                my $orientation;
                my $uniq_query1=$aln2->{r_start}-$aln1->{r_start};
                my $uniq_query2=$aln2->{r_end}-$aln1->{r_end};
                next unless ($uniq_query1>=$self->min_unique_sequence && $uniq_query2>=$self->min_unique_sequence);

                my $overlap;
                my $bkpos1=$aln2->{r_start};
                my $bkpos2=$aln1->{r_end};
                if($aln1->{refseq} eq $aln2->{refseq}){  # intra-chromosomal
                    if($aln1->{orientation}eq$aln2->{orientation}){   #indel
                        my $readspan=$aln2->{r_start}-$aln1->{r_start};
                        my $refspan=$aln2->{ref_start}-$aln1->{ref_start};
                        next if($aln1->{orientation} eq 'U' && $refspan<$self->min_unique_sequence ||   #ignore cross-mapping
                        $aln1->{orientation} eq 'C' && $refspan>-$self->min_unique_sequence);
                        $size=abs($refspan)-$readspan;
                        $orientation='+-';
                        if($size>0){
                            $type='DEL';
                            $overlap=$aln1->{r_end}-$aln2->{r_start}+1;
                            next if($overlap<-1);
                            if($aln1->{orientation} eq 'U'){
                                $refpos1=$aln1->{ref_end}-$overlap+1;
                                $refpos2=$aln2->{ref_start};
                            }
                            else{
                                $refpos1=$aln2->{ref_start}-$overlap+1;
                                $refpos2=$aln1->{ref_end};
                            }
                        }
                        else{
                            $type='INS';
                            $overlap=$aln1->{ref_end}-$aln2->{ref_start}+1;
                            if($aln1->{orientation} eq 'U'){
                                next if($overlap<-1);
                                $refpos1=$aln2->{ref_start};
                                $refpos2=$aln1->{ref_end};
                                $bkpos1=$aln1->{r_end}-$overlap+1;
                                $bkpos2=$aln2->{r_start};
                            }
                            else{
                                next if($overlap>1);
                                $refpos1=$aln1->{ref_end};
                                $refpos2=$aln2->{ref_start};
                                $bkpos1=$aln1->{r_end}-$overlap+1;
                                $bkpos2=$aln2->{r_start};
                            }
                        }
                    }
                    else{  #inversion
                        $overlap=$aln1->{r_end}-$aln2->{r_start}+1;
                        next if($overlap<-1);
                        my $readspan=abs($aln2->{r_start}-$aln1->{r_end});
                        my $refspan=abs($aln2->{ref_start}-$aln1->{ref_end});
                        $size=abs($refspan-$readspan);
                        $type='INV';
                        if($aln1->{orientation} eq 'U'){
                            $refpos1=$aln1->{ref_end}-$overlap+1;
                            #$refpos2=$aln2->{ref_start}+$overlap+1;
                            $refpos2=$aln2->{ref_start}-1;
                            $orientation='++';
                        }
                        else{
                            #$refpos2=$aln1->{ref_end}+$overlap+1;
                            $refpos2=$aln1->{ref_end}-1;
                            $refpos1=$aln2->{ref_start}-$overlap+1;
                            $orientation='--';
                        }
                    }
                }
                else{  #inter-chromosomal
                    $overlap=$aln1->{r_end}-$aln2->{r_start}+1;
                    next if($overlap<-1);
                    $refpos1=$aln1->{ref_end};
                    $refpos2=$aln2->{ref_start};
                    $size=1;
                    $type='CTX';
                    if(&GLess($aln1->{refseq}, $aln2->{refseq})){ #keep the repeat in the lower chromosome
                        if($aln1->{orientation} eq 'U'){
                            if(($aln2->{orientation} eq 'U')){
                                $orientation='+-';
                                #$refpos2=$refpos2+($overlap-1)+1;
                            }
                            else{
                                $orientation='++';
                                #$refpos2=$refpos2-($overlap-1)+1;
                            }
                            $refpos1=$refpos1-($overlap-1);
                        }
                        else{
                            if(($aln2->{orientation} eq 'U')){
                                $orientation='--';
                                #$refpos2=$refpos2+($overlap-1)+1;
                            }
                            else{
                                $orientation='-+';
                                #$refpos2=$refpos2-($overlap-1)+1;
                            }
                            $refpos1=$refpos1+$overlap-1;
                        }
                    }
                    else{
                        if($aln1->{orientation} eq 'U'){
                            if($aln2->{orientation} eq 'U'){
                                $orientation='-+';
                                #$refpos2=$refpos2+($overlap-1)+1;
                            }
                            else{
                                $orientation='++';
                                #$refpos2=$refpos2-($overlap-1)+1;
                            }
                            $refpos1=$refpos1-($overlap-1);
                        }
                        else{
                            if($aln2->{orientation} eq 'U'){
                                $orientation='--';
                                #$refpos2=$refpos2+($overlap-1)+1;
                            }
                            else{
                                $orientation='+-';
                                #$refpos2=$refpos2-($overlap-1)+1;
                            }
                            $refpos1=$refpos1+($overlap-1);
                        }
                        my $tmp=$refpos1;$refpos1=$refpos2;$refpos2=$tmp;
                        $tmp=$aln1;$aln1=$aln2;$aln2=$tmp;
                    }
                }

                if($bkpos1>$bkpos2){
                    my $tmp=$bkpos1;$bkpos1=$bkpos2;$bkpos2=$tmp;
                }
                my $var;
                ($var->{type},$var->{size},$var->{chr1},$var->{refpos1},$var->{chr2},$var->{refpos2},$var->{orientation},$var->{bkpos1},$var->{bkpos2},$var->{read},$var->{score})=($type,abs($size),$aln1->{refseq},$refpos1,$aln2->{refseq},$refpos2,$orientation,$bkpos1,$bkpos2,$read,$score);
                push @vars,$var if($score>0);
            }
        }
    }

    #Find best answer (highest score)
    my $seq;
    my $bestvar;
    my $maxscore=0;
    my $mindist=1e10;
    my ($chr1,$refpos1,$chr2,$refpos2,$pretype,$presize,$preori)=split /\_/, $self->chr_pos if(defined $self->chr_pos);
    foreach my $var(@vars){
        next if($var->{type} eq 'INV' && $self->indels_only);
        my $dist;
        if(defined $pretype && defined $presize){
            $dist=($var->{type}eq$pretype)?abs($presize-$var->{size}):1e10;
            $var->{size}=$presize if($pretype eq 'CTX');
        }
        if(defined $preori && $pretype eq 'CTX' && $preori ne $var->{orientation} ){
            $var->{score}=0;  #CTX must confirmed in the predicted orientation
        }
        next if($var->{score}<=0);

        if(!defined $bestvar ||
        (defined $pretype && $var->{type}eq$pretype && $mindist>$dist) ||  #same type and closer in size
        (!defined $pretype && ($bestvar->{score}<=$var->{score} && $bestvar->{size}<$var->{size}))  # same score but larger size
    ){
        $bestvar=$var;
        $maxscore=$var->{score};
        $mindist=$dist;
    }
    }

    if(defined $bestvar){
        #Print out
        if($bestvar->{type}=~/INS/){
            $refpos1++;
        }

        if($bestvar->{size}>=$self->min_indel_size && $bestvar->{score}>0){
            if(defined $self->chr_pos){
                if(defined $bestvar->{chr1} && $bestvar->{chr1}=~/chromosome\.(\S+)\./){ $bestvar->{chr1}=$1;}
                if(defined $bestvar->{chr2} && $bestvar->{chr2}=~/chromosome\.(\S+)\./){ $bestvar->{chr2}=$1;}
                printf "%s\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%s\t%s\t%d",$bestvar->{chr1}||$chr1,$refpos1+$bestvar->{refpos1}-1,$bestvar->{chr2}||$chr1,($refpos2||$refpos1)+($bestvar->{refpos2}||$bestvar->{refpos1})-1,$bestvar->{orientation}||'+-',$bestvar->{bkpos1},$bestvar->{bkpos2}||$bestvar->{bkpos1},$bestvar->{size},$bestvar->{type},$bestvar->{read},$bestvar->{score};
            }
            else{
                printf "%d\t%d\t%s\t%d\t%d\t%d\t%s\t%s\t%d",$bestvar->{refpos1},$bestvar->{refpos2},$bestvar->{orientation}||'+-',$bestvar->{bkpos1},$bestvar->{bkpos2}||$bestvar->{bkpos1},$bestvar->{size},$bestvar->{type},$bestvar->{read},$bestvar->{score};
            }
            if(defined $self->chr_pos){
                printf "\t%s\n",$self->chr_pos;
            }
            else{
                print "\n";
            }
        }
    }
}
sub GLess{
    my @chroms=@_;
    #resolve chromosome naming ambiguity
    for(my $i=0;$i<=$#chroms;$i++){
        my $chrom=$chroms[$i];
        if($chrom=~/\.([\w\d]+)\.fa/){
            $chroms[$i]=$1;
            $chroms[$i]=23 if($chroms[$i] eq 'X');
            $chroms[$i]=24 if($chroms[$i] eq 'Y');
        }
    }
    my $less=($chroms[0]<$chroms[1])?1:0;
    return $less;
}
