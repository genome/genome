# Annotate Assembled SVs based UCSC gene table
package Genome::Model::Tools::Sv::SvAnnot;

use strict;
use warnings;
use Genome;
use Bio::SeqIO;
use DBI;


class Genome::Model::Tools::Sv::SvAnnot {
    is => 'Genome::Model::Tools::Sv',
    has => [
        sv_file => {
            type => 'String',
            doc  => 'input sv file or files, separate by , if multiple files',
        },
        output_file  => {
            type => 'String',
            doc  => 'output sv annotation file',
        },
 
    ],
    has_optional => [
        repeat_mask => {
            type => 'Boolean',
            doc  => 'run UCSC repeat masker annotation using mysql DB access for human build36 only',
            default => 0,
        },
        sv_format => {
            type => 'String',
            doc  => 'The format of input sv file, chr1, pos1, chr2, pos2. default is standard: 1,2,4,5',
            valid_values  => ['standard', 'merged'],
            default_value => 'standard',
        },
        annot_build => {
            type => 'String',
            doc  => 'human build number, 36, 37, mouse_37 ...',
            valid_values => ["36","37","mouse_37"],
        },
        specify_chr => {
            type => 'String',
            doc  => 'Specify a single chromosome',
        },
        length_to_annot => {
            type => 'Integer',
            doc  => 'Proximity of breakpoints to gene annotations (200 bp)',
            default_value => 200,
        },
        length_to_segdup => {
            type => 'Integer',
            doc  => 'Proximity of breakpoints to segmental duplication (50 bp)',
            default_value => 50,
        },
        length_to_repeat => {
            type => 'Integer',
            doc  => 'Look for repeat annotation within +- flanking bp of breakpoint (200 bp)',
            default_value => 200,
        },
        overlap_repeat_size => {
            type => 'Integer',
            doc  => 'Overlapped size of +- flanking breakpoint and masked repeat (5 bp)',
            default_value => 5,
        },
        masked_repeat_size => {
            type => 'Integer',
            doc  => 'size of masked repeat used for overlap_repeat_size (100 bp)',
            default_value => 100,
        },
        overlap_fraction => {
            type => 'Number',
            doc  => 'Fraction of overlap (reciprocal) required to hit a dbSNP/dbVar SVs (0.5)',
            default_value => 0.5,
        },
    ],
};


sub _get_column_index {
    my $self   = shift;
    my $format = $self->sv_format;

    my %column_index = (
        standard => [1, 2, 4, 5, 7],
        merged   => [2, 4, 5, 6, 8],
    );
    unless ($column_index{$format}) {
        die $self->error_message("Input sv file format: $format is not valid");
    }

    return $column_index{$format};
}


sub _get_annot_files {
    my $self = shift;
    my $annot_num = $self->annot_build;
    my $b36_dir = Genome::Sys->dbpath("tgi/misc-annotation/human","build36-20130113")."/BreakAnnot";
    my $b37_dir = Genome::Sys->dbpath("tgi/misc-annotation/human","build37-20130113")."/BreakAnnot";
    my $m37_dir = Genome::Sys->dbpath("tgi/misc-annotation/mouse","mm9-20130310")."/BreakAnnot";

    my %annot_files;
    if ($annot_num eq '36'){
      %annot_files = (
        36 => {
            ref_gene => $b36_dir . '/Human.Mar2006.RefSeqGenes.tab',
            seg_dup  => $b36_dir . '/Human.Mar2006.SegDups.tab',
            dbsnp    => $b36_dir . '/dbsnp130.indel.named.csv',
            dbvar    => $b36_dir . '/ncbi36_submitted.gff',
            cancer_gene => $b36_dir . '/Cancer_genes.csv',
        }
      );
    }elsif($annot_num eq '37'){
      %annot_files = (
        37 => {
            ref_gene => $b37_dir . '/Human.Feb2009.RefSeqGenes.tab',
            seg_dup  => $b37_dir . '/Human.Feb2009.SegDups.tab',
            dbsnp    => $b37_dir . '/dbsnp132.indel.named.csv',
            dbvar    => $b37_dir . '/GRCh37.remap.all.germline.ucsc.gff',
            cancer_gene => $b37_dir . '/Cancer_genes.csv',
        }
      );
    }elsif($annot_num eq 'mouse_37'){
      %annot_files = (  
        mouse_37 => {
            ref_gene    => $m37_dir . '/Mouse.July2007.RefSeqgene.tab',
            cancer_gene => $m37_dir . '/Mouse_Cancer_genes.csv',
        }
      );
    }

    unless ($annot_files{$annot_num}) {
        die $self->error_message("human build number: $annot_num is not valid");
    }

    return $annot_files{$annot_num};
}


sub execute {
    my $self = shift;
    
    if ($self->repeat_mask) {
        unless ($self->annot_build eq '36') {
            $self->error_message('Currently the repeat mask mysql UCSC DB only set for human build 36');
            die;
        }
    }

    my $col_index   = $self->_get_column_index;
    my $annot_files = $self->_get_annot_files;
    my $specify_chr = $self->specify_chr;
    my $output_file = $self->output_file;

    my (%BK1s,%BK2s,%ROIs,%chrs);
    my @SVs;
    my $header;
    my @files = split /\,/, $self->sv_file;

    my $out_fh = Genome::Sys->open_file_for_writing($output_file) or die "Failed to open $output_file for writing";

    for my $sv_file (@files) {
        my $fh = Genome::Sys->open_file_for_reading($sv_file) or die "Failed to open sv_file $sv_file\n";

        while (<$fh>) {
            chomp;
            if (/^\#/) {
                $header = $_ if(/^\#/ && /Chr/i);
                next;
            }
            chomp;
            my @u = split /\s+/;
            my ($chr, $start, $chr2, $end, $type) = map{$u[$_ - 1]}@$col_index;
            next if $specify_chr and $chr ne $specify_chr;
            
            $start =~ s/\D.*//gi;
            $end   =~ s/\D.*//gi;
            next unless $start =~ /^\d+/ && $end =~ /^\d+/ && defined $chr;

            push @SVs, $_;
            $BK1s{$chr}{$start}++;
            $BK2s{$chr2}{$end}++;
            $ROIs{$chr}{$start} = $end if $chr eq $chr2;
            $chrs{$chr}++; 
            $chrs{$chr2}++;
        }
        $fh->close;
    }

    my (%AROIs, %ABKs, %SBKs, %dbBK2s, %dbVarBK2s, %RPMK);

    for my $chr (keys %chrs) {
        next if $specify_chr && $chr ne $specify_chr;
        my $annot = ReadUCSCGeneAnnotation($annot_files->{ref_gene}, $chr);

        #Select all transcripts overlapping with ROI
        my @starts = sort {$a<=>$b} keys %{$ROIs{$chr}};
        my @txEnds = sort {$a<=>$b} keys %{$$annot{$chr}};
        #push @{$Annot{$e->{chrom}}{$e->{txEnd}}{$e->{txStart}}},$e;
        for my $start (@starts) {
            my $end = $ROIs{$chr}{$start};
            while (@txEnds>0 && $start>$txEnds[0]) {
                shift @txEnds;
            }
            next unless @txEnds>0;
            for my $txEnd(@txEnds){
                my $hit=0;
                for my $txStart (keys %{$$annot{$chr}{$txEnd}}) {
	                if ($end>=$txStart) {
	                    for my $transcript (@{$$annot{$chr}{$txEnd}{$txStart}}) {
	                        push @{$AROIs{$chr}{$start}{$end}},$transcript;
                        }
	                    $hit=1;
                    }
                }
                last if $hit == 0;
            }
        }

        #Select all transcripts that overlap with the breakpoint
        my @poses = sort {$a<=>$b} (keys %{$BK1s{$chr}}, keys %{$BK2s{$chr}});
        @txEnds   = sort {$a<=>$b} keys %{$$annot{$chr}};

        my $annot_length = $self->length_to_annot;

        for my $pos (@poses) {
            while (@txEnds>0 && $pos>$txEnds[0]+$annot_length) {
                shift @txEnds;
            }
            next unless @txEnds>0;
            for my $txEnd(@txEnds) {
                for my $start (keys %{$$annot{$chr}{$txEnd}}) {
	                if ($pos>$start-$annot_length) {
	                    for my $transcript (@{$$annot{$chr}{$txEnd}{$start}}) {
	                        push @{$ABKs{$chr}{$pos}},$transcript;
                        }
                    }
                }
            }
        }

        next if $self->annot_build =~ /^mouse/;  #For mouse, only apply gene annotation and ignore everything else

        #Select all SegDup that overlaps the breakpoint
        my $segdup = ReadUCSCSegDupAnnotation($annot_files->{seg_dup}, $chr);
        my @segEnds= sort {$a<=>$b} keys %{$$segdup{$chr}};
        for my $pos(@poses) {
            while (@segEnds>0 && $pos>$segEnds[0]+$annot_length) {
                shift @segEnds;
            }
            next unless @segEnds>0;
            for my $start (keys %{$$segdup{$chr}{$segEnds[0]}}) {
                if ($pos>$start-$annot_length) {
	                for my $dup(@{$$segdup{$chr}{$segEnds[0]}{$start}}){
	                    push @{$SBKs{$chr}{$pos}},$dup;
                    }
                }
            }
        }

        #Obtain dbSNP annotation
        @poses    = sort {$a<=>$b} (keys %{$BK2s{$chr}});
        my $dbSNP = Read_dbSNPAnnotation($annot_files->{dbsnp}, $chr);
        my @chromEnds = sort {$a<=>$b} keys %{$$dbSNP{$chr}};
        for my $pos (@poses) {
            while (@chromEnds>0 && $pos>$chromEnds[0]+$annot_length) {
                shift @chromEnds;
            }
            next unless @chromEnds>0;
            for my $start (keys %{$$dbSNP{$chr}{$chromEnds[0]}}) {
                if ($pos>=$start-$annot_length) {
	                for my $var (@{$$dbSNP{$chr}{$chromEnds[0]}{$start}}){
	                    push @{$dbBK2s{$chr}{$pos}},$var;
                    }
                }
            }
        }

        my $dbVar = Read_dbVarAnnotation($annot_files->{dbvar}, $chr);
        @chromEnds=sort {$a<=>$b} keys %{$$dbVar{$chr}};
        for my $pos (@poses) {
            while (@chromEnds>0 && $pos>$chromEnds[0]+$annot_length) {
                shift @chromEnds;
            }
            next unless @chromEnds>0;
            for my $start(keys %{$$dbVar{$chr}{$chromEnds[0]}}) {
                if ($pos>=$start-$annot_length) {
	                for my $var (@{$$dbVar{$chr}{$chromEnds[0]}{$start}}) {
	                    push @{$dbVarBK2s{$chr}{$pos}},$var;
                    }
                }
            }
        }
        #Prepare repeatMasker table
        if ($self->repeat_mask) {
            my $db   = "ucsc";
            my $user = "mgg_admin";
            my $password = "c\@nc3r";
            my $dataBase = "DBI:mysql:$db:mysql2";
            my $dbh = DBI->connect($dataBase, $user, $password) || die "ERROR: Could not connect to database: $! \n";

            my $table = "chr$chr"."_rmsk";
            my $query = 
                "SELECT genoStart, genoEnd, repClass
                FROM $table
                WHERE genoEnd >= ? && genoStart <= ?
                ORDER BY genoStart";
            $RPMK{$chr} = $dbh->prepare($query) || die "Could not prepare statement '$query': $DBI::errstr \n";
        }
    }

    if (defined $header) {
        if ($self->repeat_mask) {
            $out_fh->print("$header\tRefseqGene\tDataBases\tSegDup\tRepeat\tShortIndex\n");
        }
        else {
            $out_fh->print("$header\tRefseqGene\tDataBases\tSegDup\tShortIndex\n");
        }
    }

    for my $sv (@SVs) {
        my @u   = split /\s+/, $sv;
        my @pos = map{$u[$_ - 1]}@$col_index; #($chr, $start, $chr2, $end, $type)

        my $geneAnnot = $self->GetGeneAnnotation(@pos, \%AROIs, \%ABKs);
        $out_fh->print("$sv\t$geneAnnot");

        unless ($self->annot_build =~ /^mouse/) {  #For human, apply the following crap
            pop @pos; # $type is not needed to pass on
            my $dbSNPAnnot    = $self->GetVarAnnotation(@pos, \%dbBK2s);
            my $dbVarAnnot    = $self->GetVarAnnotation(@pos, \%dbVarBK2s);
            my $dbSegDupAnnot = $self->GetSegDupAnnotation(@pos, \%SBKs);
            my $repeatAnnot   = $self->GetRepeatMaskerAnnotation(@pos, \%RPMK) if $self->repeat_mask;

            $out_fh->printf("\t%s", join(',', $dbSNPAnnot, $dbVarAnnot));
            $out_fh->printf("\t%s", $dbSegDupAnnot);
            $out_fh->printf("\t%s", $repeatAnnot) if $self->repeat_mask;
            $out_fh->printf("\tchr%s\:%d\-%d,chr%s\:%d\-%d",$pos[0],$pos[1]-500,$pos[1]+500,$pos[2],$pos[3]-500,$pos[3]+500);

            if ($pos[0] eq $pos[2]) { # $chr eq $chr2
                $out_fh->print(',chr'.$pos[0].':'.$pos[1].'-'.$pos[3]);
            }
        }
        $out_fh->print("\n");
    }
    $self->status_message('SV annotation is done.');
    return 1;
}


sub GetGeneAnnotation {
    my ($self,$chr,$start,$chr2,$end,$type,$AROI,$ABK) = @_;
    my $annot_files   = $self->_get_annot_files;
    my $cancergenelst = ReadCancerGenelst($annot_files->{cancer_gene});

    #Overlapping genes
    my %Cancergenes;
    my %OverlapGenes;

    for my $e (@{$AROI->{$chr}->{$start}->{$end}}) {
        $OverlapGenes{$e->{name2}}++;
        if (defined $cancergenelst->{uc($e->{name2})}) {
            $Cancergenes{$e->{name2}}++;
        }
    }

    #Annotate Breakpoint
    my (@e1s,@e2s);
    if (defined $ABK->{$chr}->{$start}) {
        @e1s = @{$ABK->{$chr}->{$start}};
    }
    if (defined $ABK->{$chr2}->{$end}) {
        @e2s = @{$ABK->{$chr2}->{$end}};
    }

    #select the transcript (from multiple ones)
    my ($e1,$e2,$struct1,$struct2,$annot1,$annot2);
    for (my $i=0;$i<=$#e1s;$i++) {
        for(my $j=0;$j<=$#e2s;$j++){
            if ($e1s[$i]->{name} eq $e2s[$j]->{name}) {  #same carrying transcript
	            $e1 = $e1s[$i];
	            $e2 = $e2s[$j];
            }
        }
    }
    if (!defined $e1) {
        for (my $i=0;$i<=$#e1s;$i++) {
            if (!defined $e1 || $e1->{txSize} < $e1s[$i]->{txSize}) {  #longest isoform
            	$e1 = $e1s[$i];
            }
        }
    }
    if (!defined $e2) {
        for (my $i=0;$i<=$#e2s;$i++) {
            if (!defined $e2 || $e2->{txSize} < $e2s[$i]->{txSize}) {  #longest isoform
	            $e2 = $e2s[$i];
            }
        }
    }
    $OverlapGenes{$e1->{name2}}++ if defined $e1;
    $OverlapGenes{$e2->{name2}}++ if defined $e2;

    $struct1 = AnnotExon($start,$e1) if defined $e1;
    $annot1  = defined $struct1 ? sprintf "%s%d:%d\/%d",$struct1->{unit},$struct1->{id},$struct1->{pos},$struct1->{size} : 'NA';

    $struct2 = AnnotExon($end, $e2)  if defined $e2;
    $annot2  = defined $struct2 ? sprintf "%s%d:%d\/%d",$struct2->{unit},$struct2->{id},$struct2->{pos},$struct2->{size} : 'NA';

    my @overlapgenes = keys %OverlapGenes;
    my $gene = @overlapgenes ? 'Gene:'.join('|',@overlapgenes) : '-';

    if (defined $e1 && defined $e2) {
        if ($e1->{name} eq $e2->{name}){  #same transcript
            $gene .= sprintf ",%s\|%s\:%s\-%s",$e1->{name}||'NA',$e1->{name2}||'NA',$annot1||'NA',$annot2||'NA';
            if (defined $struct1 && defined $struct2 && (abs($struct2->{id}-$struct1->{id})>0 || $struct1->{unit}=~/exon/i)) {
	            $gene.=',AffectCoding';
            }
            if (defined $struct1 && defined $struct2 && (abs($struct2->{id}-$struct1->{id})>1)) {
	            $gene.=',novelSplice';
            }
        }
        else {
            $gene .= sprintf ",%s\|%s\:%s\-%s\|%s\:%s",$e1->{name}||'NA',$e1->{name2}||'NA',$annot1||'NA',$e2->{name}||'NA',$e2->{name2}||'NA',$annot2||'NA';
            $gene .= ',AffectCoding,Fusion';
        }
    }
    elsif (defined $e1 || defined $e2) {
        $gene .= sprintf ",%s\|%s\:%s\-%s\|%s\:%s",$e1->{name}||'NA',$e1->{name2}||'NA',$annot1||'NA',$e2->{name}||'NA',$e2->{name2}||'NA',$annot2||'NA';
        $gene .= ',AffectCoding';
    }
    else{
        # neither $e1 nor $e2 is defined
    }

    my @cgenes = keys %Cancergenes;
    $gene .= ',Cancer:'.join('|',@cgenes) if $#cgenes >= 0;
    return $gene;
}


sub GetRepeatMaskerAnnotation {
    my ($self,$chr1,$pos1,$chr2,$pos2,$RPMK)=@_;
    my ($nbrpt1,$rep1) = $self->GetBKRepeatMaskerAnnotation($chr1, $pos1, $RPMK);
    my ($nbrpt2,$rep2) = $self->GetBKRepeatMaskerAnnotation($chr2, $pos2, $RPMK);
    my $repeatannot = '-';
    my $rm_size = $self->masked_repeat_size;

    if ($nbrpt1>$rm_size || $nbrpt2>$rm_size) {
        $repeatannot = sprintf "Repeat:%s-%s", $rep1 || 'NA',$rep2 ||'NA';
    }
    return $repeatannot;
}


sub GetBKRepeatMaskerAnnotation {
    my ($self, $chr1, $pos, $RPMK) = @_;
    if ($chr1 =~ /[MN]/) {#MT, N_xxxx?
        return (0, undef);
    }

    my $rp_length  = $self->length_to_repeat;
    my $rp_overlap = $self->overlap_repeat_size;

    my $start = $pos - $rp_length;
    my $stop  = $pos + $rp_length;
    my $db_prepare = $RPMK->{$chr1};
    $db_prepare->execute($start, $stop) || die "Could not execute statement for repeat masker table with (chr$chr1, $start, $stop): $DBI::errstr \n";
    
    my %repCount;
    while (my ($chrStart, $chrStop, $repClass) = $db_prepare->fetchrow_array()) {
        my $start_last = $chrStart > $start ? $chrStart : $start;
        my $stop_last  = $chrStop < $stop   ? $chrStop  : $stop;
        $repCount{$repClass} = $stop_last-$start_last+1 if $start_last-$rp_overlap <= $pos && $pos <= $stop_last+$rp_overlap;
    }
    
    my @sortedrepClass;
    my $maxClass;
    my $maxClassCount=0;
    for (keys %repCount) {
        if (defined $repCount{$_} && $repCount{$_}>$maxClassCount) {
            $maxClass = $_;
            $maxClassCount = $repCount{$_};
        }
    }
    return ($maxClassCount, $maxClass);
}


sub GetSegDupAnnotation {
    my ($self,$chr,$start,$chr2,$end,$SBK) = @_;
    my (@e1s, @e2s);

    if (defined $SBK->{$chr}->{$start}) {
        @e1s = @{$SBK->{$chr}->{$start}};
    }
    if (defined $SBK->{$chr2}->{$end}) {
        @e2s = @{$SBK->{$chr2}->{$end}};
    }

    #select the transcript (from multiple ones)
    my ($e1,$e2,$struct1,$struct2,$annot1,$annot2);

    for (my $i=0;$i<=$#e1s;$i++) {
        for(my $j=0;$j<=$#e2s;$j++) {
            if($e1s[$i]->{name} eq $e2s[$j]->{name}){  #same carrying transcript
	            $e1 = $e1s[$i];
	            $e2 = $e2s[$j];
            }
        }
    }
    if (!defined $e1) {
        for (my $i=0;$i<=$#e1s;$i++) {
            if(!defined $e1 || $e1->{Size} < $e1s[$i]->{Size}){
	            $e1 = $e1s[$i];
            }
        }
    }
    if(!defined $e2){
        for (my $i=0;$i<=$#e2s;$i++) {
            if (!defined $e2 || $e2->{Size} < $e2s[$i]->{Size}) {
	            $e2=$e2s[$i];
            }
        }
    }

    my $segdup = '-';
    if (defined $e1 || defined $e2) {
        $segdup = sprintf "SegDup:%s\-%s",$e1->{name}||'NA',$e2->{name}||'NA';
    }
    return $segdup;
}


sub GetVarAnnotation {
    my ($self,$chr,$start,$chr2,$end,$db) = @_;
    return '-' if $chr ne $chr2;
    my @vars;
    if (defined $$db{$chr2}{$end}) {
        @vars=@{$$db{$chr2}{$end}};
    }
    my $bestvar;
    my ($maxratio1,$maxratio2) = (0,0);
    my $frac = $self->overlap_fraction;

    for my $var(@vars) {
        my $pos1 = $end<$var->{chromEnd} ? $end : $var->{chromEnd};
        my $pos2 = $start>$var->{chromStart} ? $start : $var->{chromStart};
        my $overlap = $pos1-$pos2+1;
        my $ratio1 = $overlap/(abs($end-$start)+1);
        my $ratio2 = $overlap/(abs($var->{chromEnd}-$var->{chromStart})+1);
        if ($ratio1 >= $frac && $ratio2 >= $frac && ($ratio1 >= $maxratio1 || $ratio2 >= $maxratio2)) {
            $bestvar   = $var;
            $maxratio1 = $ratio1;
            $maxratio2 = $ratio2;
        }
    }
    my $varreport = '-';
    if (defined $bestvar) {
        $varreport = $bestvar->{name};
    }
    return $varreport;
}


sub AnnotExon {
    my ($pos, $e) = @_;
    my @exonStarts = split /\,/, $e->{exonStarts};
    my @exonEnds   = split /\,/, $e->{exonEnds};

    my $report;
    for (my $i=0;$i<=$#exonStarts;$i++) {
        if ($exonStarts[$i]<=$pos && $pos<=$exonEnds[$i]) {
            my $id=$i+1;
            my $rpos=$pos-$exonStarts[$i]+1;
            if ($e->{strand} eq '-') {
	            $id = $#exonStarts-$i+1;
	            $rpos=$exonEnds[$i]-$pos+1;
            }
            ($report->{unit},$report->{id},$report->{pos},$report->{size}) = ('Exon',$id,$rpos,$exonEnds[$i]-$exonStarts[$i]+1);
        }
        elsif ($i>0 && $exonEnds[$i-1]<$pos && $pos<$exonStarts[$i]) { 
            my $id   = $i;
            my $rpos = $pos-$exonEnds[$i-1]+1;
            if ($e->{strand} eq '-') {
	            $id   = $#exonStarts-$i+1;
	            $rpos = $exonStarts[$i]-$pos+1;
            }
            ($report->{unit},$report->{id},$report->{pos},$report->{size}) = ('Intron',$id,$rpos,$exonStarts[$i]-$exonEnds[$i-1]+1);
        }
    }
    return $report;
}


sub ReadUCSCGeneAnnotation {
    my ($file, $chr) = @_;
    my %Annot;
    open (AN, "<$file") || die "Unable to open UCSC gene annotation: $file\n";
    while (<AN>) {
        chomp;
        next if /^\#/;
        my $e;
        ($e->{bin},$e->{name},$e->{chrom},$e->{strand},$e->{txStart},$e->{txEnd},$e->{cdsStart},$e->{cdsEnd},$e->{exonCount},$e->{exonStarts},$e->{exonEnds},$e->{id},$e->{name2},$e->{cdsStartStat},$e->{cdsEndStat},$e->{exonFrames}) = split;
        $e->{chrom} =~ s/chr//;
        next unless $e->{chrom} eq $chr;
        $e->{txSize} = abs($e->{txStart}-$e->{txEnd}+1);
        push @{$Annot{$e->{chrom}}{$e->{txEnd}}{$e->{txStart}}}, $e;
    }
    close AN;
    return \%Annot;
}


sub ReadUCSCSegDupAnnotation {
    my ($file, $chr) = @_;
    my %SegDup;
    open (SEGDUP,"<$file") || die "Unable to open UCSC seg dup annotation: $file\n";
    while (<SEGDUP>) {
        chomp;
        next if /^\#/;
        my $e;
        my @extra;
        ($e->{bin},$e->{chrom},$e->{Start},$e->{End},$e->{name},$e->{score},$e->{strand},@extra) = split;
        $e->{chrom} =~ s/chr//;
        next unless $e->{chrom} eq $chr;
        $e->{Size} = abs($e->{Start}-$e->{End}+1);
        push @{$SegDup{$e->{chrom}}{$e->{End}}{$e->{Start}}}, $e;
    }
    close SEGDUP;
    return \%SegDup;
}


sub Read_dbSNPAnnotation{
    my ($file, $chr) = @_;
    my %dbSNP;
    open (DBSNP, "<$file") || die "Unable to open dbSNP annotation: $file\n";
    while (<DBSNP>) {
        chomp;
        next if /^\#/;
        my $p;
        ($p->{bin},$p->{chrom},$p->{chromStart},$p->{chromEnd},$p->{name},$p->{score},$p->{strand},$p->{refNCBI},$p->{refUCSC},$p->{observed},$p->{molType},$p->{class},$p->{valid},$p->{avHet},$p->{avHetSE},$p->{func},$p->{locType},$p->{weight}) = split /\t+/;
        $p->{chrom} =~ s/chr//;
        next unless $p->{chrom} eq $chr;
        push @{$dbSNP{$p->{chrom}}{$p->{chromEnd}}{$p->{chromStart}}}, $p;
    }
    close DBSNP;
    return \%dbSNP;
}


sub Read_dbVarAnnotation{
    my ($file, $chr) = @_;
    my %dbVar;
    open (DBVAR,"<$file") || die "Unable to open dbVar annotation: $file\n";
    while (<DBVAR>) {
        chomp;
        next if /^\#/;
        my ($p,$db,$tmp1,$tmp2,$tmp3);
        ($p->{chrom},$db,$p->{var},$p->{chromStart},$p->{chromEnd},$tmp1,$tmp2,$tmp3,$p->{name}) = split /\t+/;
        $p->{chrom} =~ s/chr//;
        next unless $p->{chrom} eq $chr;
        push @{$dbVar{$p->{chrom}}{$p->{chromEnd}}{$p->{chromStart}}}, $p;
    }
    close DBVAR;
    return \%dbVar;
}


sub ReadCancerGenelst{
    my @fin = @_;
    my %genes;
    for my $f (@fin) {
        open (FIN, "<$f") || die "unable to open cancer gene list: $f\n";
        $_ = <FIN>;
        while (<FIN>) {
            next if /^\#/;
            chomp;
            my @t = split;
            $genes{$t[0]}++;
        }
    }
    return \%genes;
}

1;
