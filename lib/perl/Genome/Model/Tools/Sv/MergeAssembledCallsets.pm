package Genome::Model::Tools::Sv::MergeAssembledCallsets;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sv::MergeAssembledCallsets {
    is => 'Command::V2',
    has => [
        index_file => { #ARGV[0]
            is => 'Text',
            doc => 'the result index file (three columns: set name, .csv file, .fasta file',
        },
        output_file => {
            is => 'Text',
            default_value => '-',
            doc => 'Where to write the output (defaults to STDOUT)',
        },
        max_merge_distance => { #-d
            is => 'Integer',
            default_value => 20,
            doc => 'do not merge SVs unless they differ less than this many basepairs in start and end positions, subject to 1 duplication',
        },
        max_merge_size_difference => { #-p
            is => 'Integer',
            default_value => 20,
            doc => 'do not merge SVs unless they differ less than this many basepairs in size',
        },
        min_size => { #-s
            is => 'Integer',
            default_value => 10,
            doc => 'minimal basepair size to trust',
        },
        max_size_difference => { #-z
            is => 'Float',
            default_value => '0.1',
            doc => 'ignore SVs with more than this size difference from the expected size',
        },
        min_weighted_assembly_score => { #-w
            is => 'Integer',
            default_value => 50,
            doc => 'filter out prediction with weighted assembly score smaller than this',
        },
        cna_threshold => { #-C
            is => 'Float',
            default_value => '0.5',
            doc => 'always include SVs with CNA > this',
        },
        max_location_difference => { #-n
            is => 'Integer',
            default_value => '1000',
            doc => 'ignore SVs with assembled breapoint more this many bp off the predicted breakpoint',
        },
        output_fasta => { #-f
            is => 'Text',
            doc => 'dump breakpoint sequences to this fasta file',
        },
        use_hq_alignment_filter => { #-h
            is => 'Boolean',
            default_value => 0,
            doc => 'apply a high quality alignment filter on the calls',
        },
        sort_by_assembly_scores => { #-c
            is => 'Boolean',
            default_value => 0,
            doc => 'sort result by assembly score',
        },
    ],
    has_optional => [
        exclude_file => { #-e
            is => 'Text',
            doc => 'exclude assembled SV overlap with Breakdancer SVs in this file',
            is_optional => 1,
        },
        amibguity_allowed => { #-l
            is => 'Integer',
            default_value => 200,
            doc => 'used with "exclude_file" option: ambiguity allowed at breakpoints, this many basepairs',
        },
        overlap_allowed => { #-r
            is => 'Float',
            default_value => '0.8',
            doc => 'used with "exclude_file" option: fraction of overlap between two events',
        },
        name_filter => { #-L
            is => 'Text',
            doc => 'ignore SVs that are detected in sets with this string in the name',
        },
    ],
};

sub execute {
    my $self = shift;

    my %f_csv;
    my %f_fasta;
    my %diff_size;
    my %diff_loc;
    my $index_file = $self->index_file;
    open(INX,"<$index_file") || die "unable to open $index_file";
    while(<INX>){
        chomp;
        next if(/^\#/);
        my ($set,$f_indel,$f_fasta,$diff_size,$diff_loc)=split;
        push @{$f_csv{$set}},$f_indel;
        #push @{$f_fasta{$set}},$f_fasta;
        $f_fasta{$set}=$f_fasta;
        $diff_size{$set}=$diff_size if(defined $diff_size);
        $diff_loc{$set}=$diff_loc if(defined $diff_loc);
    }

    my $output_fh;
    if($self->output_file and $self->output_file ne '-') {
        $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    } else {
        $output_fh = \*STDOUT;
    }

    my %SVs;
    foreach my $set(keys %f_csv){
        my @files=@{$f_csv{$set}};
        foreach my $file(@files){
            $self->_add_svs(\%diff_size, \%diff_loc, \%SVs, $set,$file);
        }
    }

    my %EXSVs;
    if($self->exclude_file){  #exclude SVs
        my $exclude_file = $self->exclude_file;
        open(SVIN,"<$exclude_file") || die "unable to open $exclude_file\n";
        while(<SVIN>){
            next if(/^\#/);
            chomp;
            my $sv;
            my @extra;
            ($sv->{chr1},$sv->{start},$sv->{ori1},$sv->{chr2},$sv->{end},$sv->{ori2},$sv->{type},$sv->{size},$sv->{qual},$sv->{nrps},$sv->{rp_sample},$sv->{allele_frequency},@extra)=split /\s+/;
            next unless(defined $sv->{type} && defined $sv->{chr1} && defined $sv->{chr2});
            push @{$EXSVs{$sv->{type}}{$sv->{chr1}}},$sv;
        }
        close(SVIN);
        foreach my $chr(sort keys %SVs){
            my $idx=1;
            foreach my $chr2(sort keys %{$SVs{$chr}}){
                foreach my $pos(sort {$a<=>$b} keys %{$SVs{$chr}{$chr2}}){
                    my $SV=$SVs{$chr}{$chr2}{$pos};
                    my $hit= $self->_overlap(\%EXSVs, $chr,$chr2,$pos,$SV);
                    if($hit){
#	  $output_fh->printf( "exclude: %s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n",$chr,$SV->{OuterStart},$SV->{InnerStart},$chr2,$SV->{InnerEnd},$SV->{OuterEnd},$SV->{type},join(',',@{$SV->{ori}}),$SV->{minsize},$SV->{maxsize},join(',',@{$SV->{group}}),join(',',@{$SV->{wAsmscore}}));
                        delete $SVs{$chr}{$chr2}{$pos};
                    }
                }
            }
        }
    }

#print table
#$output_fh->printf( "#ID\tCHR\tPOS\tTYPE\tSIZE\tGROUPS\tHET\n");
#$output_fh->printf( "#ID\tCHR1\tOUTER_START\tINNER_START\tCHR2\tINNER_END\tOUTER_END\tTYPE\tORIENTATION\tMINSIZE\tMAXSIZE\tSOURCE\tSCORES\tCopy_Number\tGene\tKnown\n");
    $output_fh->printf( "#ID\tCHR1\tOUTER_START\tINNER_START\tCHR2\tINNER_END\tOUTER_END\tTYPE\tORIENTATION\tMINSIZE\tMAXSIZE\tSOURCE\tSCORES\tCopy_Number\n");

    my $fout;
    if($self->output_fasta){
        my $output_fasta = $self->output_fasta;
        if( -s $output_fasta){
            die("Output FASTA file $output_fasta already exists.");
        }
        $fout = Bio::SeqIO->new(-file => ">>$output_fasta" , '-format' => 'Fasta');
    }

    my %SVscores;
    my $filter = $self->name_filter;
    foreach my $chr(sort keys %SVs){
        my $idx=1;
        foreach my $chr2(sort keys %{$SVs{$chr}}){
            foreach my $pos(sort {$a<=>$b} keys %{$SVs{$chr}{$chr2}}){
                my $SV=$SVs{$chr}{$chr2}{$pos};
                my $cid="$chr.$idx";
                #my $pos=($SV->{OuterStart}==$SV->{InnerStart})?$SV->{OuterStart}:join(',',$SV->{OuterStart},$SV->{InnerStart});
                my $size=($SV->{minsize}==$SV->{maxsize})?$SV->{minsize}:join('-',$SV->{minsize},$SV->{maxsize});
                #$output_fh->printf( "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$cid,$chr,$pos,"Deletion",$size,join(',',@{$SV->{group}}),join(',',@{$SV->{het}}));
                my $maxscore=-1;
                my $winningGroup;
                my $winningContig;
                my $ignore=0;
                for(my $i=0;$i<=$#{$SV->{wAsmscore}};$i++){
                    if($maxscore<${$SV->{wAsmscore}}[$i]){
                        $maxscore=${$SV->{wAsmscore}}[$i];
                        $winningGroup=${$SV->{group}}[$i];
                        $winningContig=${$SV->{prestr}}[$i];
                    }
                    $ignore=1 if(defined $filter && ${$SV->{group}}[$i]=~ /$filter/i);
                }
                next if($ignore);

                my $keys=join('|',$cid,$chr,$chr2,$pos,$winningGroup,$winningContig);
                $SVscores{$keys}=$maxscore;
                $idx++;
            }
        }
    }

    my @SVkeys;
    if($self->sort_by_assembly_scores){
        @SVkeys=sort {$SVscores{$b}<=>$SVscores{$a}} keys %SVscores;
    }
    else{
        @SVkeys=sort keys %SVscores;
    }

    foreach my $key(@SVkeys){
        my ($cid,$chr,$chr2,$pos,$winningGroup,$winningContig)=split /\|/,$key;
        my $SV=$SVs{$chr}{$chr2}{$pos};
        $output_fh->printf( "%s\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n",$cid,$chr,$SV->{OuterStart},$SV->{InnerStart},$chr2,$SV->{InnerEnd},$SV->{OuterEnd},$SV->{type},join(',',@{$SV->{ori}}),$SV->{minsize},$SV->{maxsize},join(',',@{$SV->{group}}),join(',',@{$SV->{wAsmscore}}),$SV->{extra});

        if($self->output_fasta){
            my $seqobj;
            my $contigfile=$f_fasta{$winningGroup};
            my $in  = Bio::SeqIO->newFh(-file => "$contigfile" , '-format' => 'Fasta');
            while ( my $seq = <$in> ) {
                # do something with $seq
                my ($id,$var,$ins,$strand,$score)=split /\,/,$seq->id;
                $id=~s/ID\://;
                next unless($id eq $winningContig);
                my @u=split /\,/,$seq->id;
                shift @u;
                my $sequence=$seq->seq();
                my $dispid=join(',',$cid,$winningGroup,@u);
                $seqobj = Bio::Seq->new( -display_id => "$dispid",
                    -seq => $sequence);
                last;
            }
            if(defined $seqobj){
                $fout->write_seq($seqobj);
            }
            else{
                print STDERR "Contig $winningContig not found in $contigfile\n";
            }
        }

    }

    return 1;
}

sub _add_svs {
    my $self = shift;
    my $diff_size = shift;
    my $diff_loc = shift;
    my $SVs = shift;
    my ($set,$fin)=@_;
    open(FIN,"<$fin") || die "unable to open $fin\n";
    while(<FIN>){
        chomp;
        my ($chr,$start,$chr2,$end,$ori,$size,$type,$het,$wAsmscore,$read_len,$perc_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$microhomology,$scarstr,$prestr,$asm_parm,@extra)=split /\s+/;
        next unless(defined $prestr && $size=~/^\d/);
        my ($pchr1,$pre_pos1,$pchr2,$pre_pos2,$pre_type,$pre_size,$pre_ori)=split /\./,$prestr;
        my $size_diff_cutoff=$diff_size->{$set} || $self->max_size_difference;
        $size_diff_cutoff=1e10 if($type && $type=~/ctx/i);
        my $loc_diff_cutoff=$diff_loc->{$set} || $self->max_location_difference;
        $start=~s/\(\d+\)//g;
        $end=~s/\(\d+\)//g;
        $size=~s/\(\d+\)//g;
        $type=~s/\(\w+\)//g;
        my $extra_info=join("\t",@extra);
        my $f_cna=0;

        if(@extra && $extra[0]=~/\.bam/){
            my @cnstr=split /\,/,$extra[0];
            my ($ncn)=($cnstr[0]=~/\:(\S+)$/);
            if($ncn!~/NA/i){
                for(my $i=1;$i<=$#cnstr;$i++){
                    my ($cn)=($cnstr[$i]=~/\:(\S+)$/);
                    $f_cna=$cn-$ncn if($cn!~/NA/i);
                }
            }
        }

        #filtering criteria

        $perc_aligned=~s/\%//;
        my $flanksize=$read_len*$perc_aligned/100;
        my $sub_rate=$n_sub/$flanksize;
        my $subindel_rate=$n_indel/$flanksize;

        if($self->use_hq_alignment_filter &&  ($n_seg>2
                || $n_seg==2 && (($sub_rate>0.006 && $n_sub>2) || ($subindel_rate> 0.002 && $n_indel>1)|| $nbp_indel>20)
                || $n_seg==1 && ($sub_rate>0.006 || $subindel_rate>0.001 || $nbp_indel>5)
                #|| $n_seg==2 && (($sub_rate>0.006 && $n_sub>2) || ($subindel_rate> 0.002 && $n_indel>1)|| $nbp_indel>20) && ($size<=99999) && (abs($f_cna)<$self->cna_threshold)
                #|| $n_seg==1 && ($sub_rate>0.006 || $subindel_rate>0.001 || $nbp_indel>5)&& ($size<=99999) && (abs($f_cna)<$self->cna_threshold)
                || ($type ne 'INS' ) && ($perc_aligned<80)
                #|| ($type eq 'DEL' ) && ($size>99999) && ($f_cna>-$self->cna_threshold)
                #|| ($type eq 'ITX' ) && ($size>99999) && ($f_cna<$self->cna_threshold)
                #|| ($type eq 'INV' ) && ($size>99999)
                #|| $microhomology>=200
                || $size<$self->min_size || $type ne $pre_type || ($pre_size>0 && abs($size-$pre_size)/$pre_size>$self->max_size_difference) || abs($start-$pre_pos1)>$loc_diff_cutoff
                || $wAsmscore >0 && $wAsmscore<$self->min_weighted_assembly_score
            )
        ){
            printf STDERR  "%s\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d%%\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",$set,$chr,$start,$chr2,$end,$ori,$size,$type,$het,$wAsmscore,$read_len,$perc_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$microhomology,$scarstr,$prestr,$asm_parm,$extra_info;
            next;
        }

        #$output_fh->printf( "%s\t%d\t%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d%%\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n", $chr,$start,$chr2,$end,$ori,$size,$type,$het,$wAsmscore,$read_len,$perc_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$microhomology,$scarstr,$prestr,$asm_parm);

        my $found;
        my @tds=($type && $type=~/ctx/i)?(0):(-$size,0,$size);
        foreach my $td(@tds){
            my $locrange=int($size*0.01);
            my $sizerange=$locrange;
            #$locrange=($locrange>$self->max_merge_distance)?$locrange:$self->max_merge_distance;
            #$sizerange=($sizerange>$self->max_merge_size_difference)?$sizerange:$self->max_merge_size_difference;
            $locrange=$self->max_merge_distance;
            $sizerange=$self->max_merge_size_difference;
            for (my $i=-$locrange+$td;$i<=$locrange+$td;$i++){
                my $SV=$SVs->{$chr}{$chr2}{$start+$i};
                if(defined $SV
                    && ($SV->{maxsize}<=$size+$sizerange && $SV->{minsize}>=$size-$sizerange || $type=~/ctx/i)
                    && $end>=$SV->{OuterEnd}-$locrange && $end<=$SV->{InnerEnd}+$locrange
                ){
                    $found=$i;
                    last;
                }
            }
        }

        my $SVregion;

        if(defined $found){
            $SVregion=$SVs->{$chr}{$chr2}{$start+$found};
            $SVregion->{minsize}=($SVregion->{minsize}>$size)?$size:$SVregion->{minsize};
            $SVregion->{maxsize}=($SVregion->{maxsize}<$size)?$size:$SVregion->{maxsize};
            $SVregion->{OuterStart}=($SVregion->{OuterStart}>$start)?$start:$SVregion->{OuterStart};
            $SVregion->{InnerStart}=($SVregion->{InnerStart}<$start)?$start:$SVregion->{InnerStart};
            $SVregion->{OuterEnd}=($SVregion->{OuterEnd}<$end)?$end:$SVregion->{OuterEnd};
            $SVregion->{InnerEnd}=($SVregion->{InnerEnd}>$end)?$end:$SVregion->{InnerEnd};
        }
        else{
            $SVregion->{minsize}=$size;
            $SVregion->{maxsize}=$size;
            $SVregion->{OuterStart}=$start;
            $SVregion->{InnerStart}=$start;
            $SVregion->{OuterEnd}=$end;
            $SVregion->{InnerEnd}=$end;
            $SVregion->{type}=$type;
            $SVregion->{extra}=$extra_info;
        }

        my $exist=0;
        for(my $ki=0;$ki<=$#{$SVregion->{group}};$ki++){
            my $group=${$SVregion->{group}}[$ki];
            if($group eq $set){
                $exist=1;
                if($wAsmscore>${$SVregion->{wAsmscore}}[$ki]){  #select the one with higher assembly score
                    ${$SVregion->{wAsmscore}}[$ki]=$wAsmscore;
                    ${$SVregion->{prestr}}[$ki]=$prestr;
                    ${$SVregion->{het}}[$ki]=$het;
                    ${$SVregion->{ori}}[$ki]=$ori;
                }
                last;
            }
        }
        if(! $exist){
            push @{$SVregion->{group}},$set;
            push @{$SVregion->{het}},$het;
            push @{$SVregion->{wAsmscore}},$wAsmscore;
            push @{$SVregion->{prestr}},$prestr;
            push @{$SVregion->{ori}},$ori;
            $SVs->{$chr}{$chr2}{$start+($found||0)}=$SVregion;
        }
    }
    close(FIN);
}

sub _overlap{
    my $self = shift;
    my $EXSVs = shift;
    my ($chr1,$chr2,$pos,$SV)=@_;
    my ($pt,$pc1,$ps,$pc2,$pe)=($SV->{type},$chr1,($SV->{OuterStart}+$SV->{InnerStart})/2,$chr2,($SV->{InnerEnd}+$SV->{OuterEnd})/2);
    my $hit=0;
    if(defined $EXSVs->{$SV->{type}}{$chr1}){
        my @bsvs=@{$EXSVs->{$SV->{type}}{$chr1}};
        foreach my $bsv(@bsvs){
            my ($t,$c1,$s,$c2,$e)=($bsv->{type},$bsv->{chr1},$bsv->{start},$bsv->{chr2},$bsv->{end});
            if(defined $pe && defined $e && $pt eq $t){ # same SV type
                my $breakpoints_overlap=($pc1 eq $c1 && $pc2 eq $c2 && abs($ps-$s)<=$self->ambiguity_allowed && abs($pe-$e)<=$self->ambiguity_allowed)?1:0;
                my $interval_overlap=0;
                if($c1 eq $c2){  #same chromosome
                    ($ps,$pe)=sort {$a<=>$b} ($ps,$pe);
                    ($s,$e)=sort {$a<=>$b} ($s,$e);
                    my $length=$pe-$ps+1;
                    my $size=$e-$s+1;
                    my $overlap=(($e<$pe)?$e:$pe)-(($s>$ps)?$s:$ps)+1;
                    my $perc1=$overlap/$size;
                    my $perc2=$overlap/$length;
                    $interval_overlap=(defined $self->overlap_allowed && $perc1>=$self->overlap_allowed && $perc2>=$self->overlap_allowed)?1:0;
                }
                if($breakpoints_overlap || $interval_overlap){
                    $hit=1;
                    last;
                }
            }
        }
    }

    return $hit;
}

1;
