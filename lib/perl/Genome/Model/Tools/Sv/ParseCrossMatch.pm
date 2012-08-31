package Genome::Model::Tools::Sv::ParseCrossMatch;


use strict;
use warnings;
use Genome::Model::Tools::Sv;
require Genome::Sys;

our $VERSION = $Genome::Model::Tools::Sv::VERSION;

class Genome::Model::Tools::Sv::ParseCrossMatch {
    is  => 'Genome::Model::Tools::Sv',
    has => [
       input_file    => {
           type => 'String',
           doc  => 'Input cross_match output file',
       },
       min_base_qual => {
           type => 'Number',
           doc  => 'minimum base quality',
           default_value => 15,
           is_optional   => 1,
       },
       homopolymer_indel_size => {
           type => 'Number',
           doc  => 'homo polymer indel size',
           default_value => 2,
           is_optional   => 1,
       },
       _align => {
           type => 'HASH',
           is_optional => 1,
       },
       dcpos => {
           type => 'HASH',
           is_optional => 1,
       },
    ],
};

sub create{
    my $class = shift;
    my $self  = $class->SUPER::create(@_) or return;

    unless (-s $self->input_file) {
        $self->error_message('Input cross_match file: '.$self->input_file.' is invalid');
        return;
    }
    $self->Read;
    return $self;
}


sub Read {
    my $self = shift;
    
    my %DCReads;
    my %DCPoses;
    my @discrepancies;

    my ($aln, $paln_readname, $tmpstr, $score, $perc_sub, $perc_del, $perc_ins, $readname, $r1, $r2, $r3, $comp, $refseq, $g1, $g2, $g3);

    my $fh = Genome::Sys->open_file_for_reading($self->input_file) or return;

    while (my $line = $fh->getline) {
        chomp $line;
        if ($line =~ /^ALIGNMENT/){
            if (defined $aln) {
	            if($#discrepancies >= 0) {
	                my @disc = @discrepancies;
	                $aln->{discrepancy} = \@disc;
                }
	            push @{$DCReads{$paln_readname}->{aln}}, $aln;
	            undef @discrepancies;
	            undef $aln;
            }

            my @u = split /\s+/, $line;
            $aln->{line} = join '|', @u;
            pop @u if $u[$#u] eq '*';  #include both higher and lower scoring match

            if ($#u == 12) {
	            ($tmpstr,$score,$perc_sub,$perc_del,$perc_ins,$readname,$r1,$r2,$r3,$refseq,$g1,$g2,$g3) = @u;
	            undef $comp;
            }
            elsif ($#u >= 13) {
	            ($tmpstr,$score,$perc_sub,$perc_del,$perc_ins,$readname,$r1,$r2,$r3,$comp,$refseq,$g1,$g2,$g3) = @u;
            }
            else {
            }
            $paln_readname = $readname;

            my ($r_start, $r_end, $r_rest);
            
            if ($r1=~/\(/) {
	            ($r_start, $r_end) = ($r2, $r3);
	            ($r_rest) = ($r1 =~ /\((\d+)\)/);
            }
            else {
	            ($r_start, $r_end) = ($r1, $r2);
	            ($r_rest) = ($r3 =~ /\((\d+)\)/);
            }

            my $rlen = $r_end + $r_rest;
            my ($ref_start, $ref_end, $ref_rest);
            
            if ($g1 =~ /\(/){
	            ($ref_start, $ref_end) = ($g2, $g3);
	            ($ref_rest) = ($g1 =~ /\((\d+)\)/);
            }
            else {
	            ($ref_start, $ref_end) = ($g1, $g2);
	            ($ref_rest) = ($g3 =~ /\((\d+)\)/);
            }

            $aln = {
                %$aln,
                score => $score,
                perc_sub    => $perc_sub,
                perc_del    => $perc_del,
                perc_ins    => $perc_ins,
                r_start     => $r_start,
                r_end       => $r_end,
                r_rest      => $r_rest,
                refseq      => $refseq,
                ref_start   => $ref_start,
                ref_end     => $ref_end,
                ref_rest    => $ref_rest,
                orientation => $comp || 'U',
            };
        }
        elsif ($line =~ /DISCREPANCY/) {
            my @u = split /\s+/, $line;
            if ($u[1] =~ /([SDI])\-*(\d*)/) {  #substitutions and indels
	            my ($dctype, $dcsize) = ($1,$2);
	            $dcsize = 1 if !defined $dcsize || length($dcsize) <= 0;
	            shift @u;
	            push @discrepancies, join"\t", @u;
	
	            my $dcpos = $u[3];
            	my ($base, $qual) = ($u[2] =~ /(\S)\((\d+)\)/);
	            if (defined $comp && $comp eq 'C'){  #reverse complement
	                $base =~ tr/ACGT/TGCA/;
	                #complement indel position
	                if($dctype eq 'D'){
	                    #$dcpos=$dcpos-$dcsize+1;
	                    #$dcpos--;
                    }
#	                elsif($dctype eq 'I'){
#	                    $dcpos=$dcpos-1;
#	                }
#	                else{}
                }

	            my $dcinfo = {
                    type       => $u[0],
                    rpos       => $u[1],
                    base       => $base,
                    qual       => $qual,
                    read       => $readname,
                    score      => $score,
                    aln_orient => $comp || 'U',
                };
                
	            my $maxlen_homopolymer = _Indel_454_homopolymerErr($u[$#u], $comp);
                
	            if ($qual >= $self->min_base_qual && #Require minimal quality score to be trusted
	                (
	                    ($dctype eq 'S') ||
	                    ($dctype=~/[DI]/ && $dcsize > 1) ||
	                    ($dctype=~/[DI]/ && $dcsize == 1 && $maxlen_homopolymer < $self->homopolymer_indel_size)
                    )
                ) {
	                my @DClist;
	                if (defined $DCPoses{$dcpos}) {
	                    @DClist = @{$DCPoses{$dcpos}};
                    }
	                push @DClist, $dcinfo;
	                $DCPoses{$dcpos} = \@DClist;
                }
            }
        }
        elsif ($line =~ /SCORE_HISTOGRAM/) {
            my @u = split /\s+/;
            my ($bestscore, $nbest) = ($u[$#u] =~ /(\d+)\((\d+)\)/);
            my ($sec_bestscore, $nsec_best) = ($u[$#u-1] =~ /(\d+)\((\d+)\)/);
            
            $DCReads{$u[1]} = {
                best      => $bestscore,
                nbest     => $nbest,
                sec_best  => $sec_bestscore || 0,
                nsec_best => $nsec_best || 0,
                associate => $readname,
            }
        }
        else{}
    }
    if (defined $aln) {
        $aln->{discrepancy} = \@discrepancies if $#discrepancies >= 0;
        push @{$DCReads{$paln_readname}->{aln}}, $aln;
    }
    
    $fh->close;
    
    $self->_align(\%DCReads);
    $self->dcpos(\%DCPoses);
    
    return (\%DCReads, \%DCPoses);
}


sub _Indel_454_homopolymerErr{
  #Test find if a 1bp indel is caused by upstream homopolymer errors
    my ($string, $comp) = @_;
    my @bases  = split //, $string;
    my $midpos = int($#bases/2);
    my $maxlen_homopolymer = 0;

    for my $start ($midpos, $midpos+1) {
        my $pbase = $bases[$start];
        my $end;
        for ($end=$start+1;$end<=$#bases;$end++) {
            last if $bases[$end] ne $pbase;
        }
        my $len = $end - $start;
        $maxlen_homopolymer = ($maxlen_homopolymer < $len) ? $len : $maxlen_homopolymer;
    }

    for my $start ($midpos, $midpos-1) {
        my $pbase = $bases[$start];
        my $end;
        for ($end=$start-1; $end>=0; $end--){
            last if $bases[$end] ne $pbase;
        }
        my $len = $start - $end;
        $maxlen_homopolymer = ($maxlen_homopolymer<$len) ? $len : $maxlen_homopolymer;
    }

    return $maxlen_homopolymer;
}


sub GetAlleleInfo{
    my ($self, $refseq, $f_qual, $refpos, $Min_Reads) = @_;
    
    my $refbase      = substr $refseq, $refpos-1, 1;
    my $DCreadsCount = 0;
    my %Allele_info;
    
    if (defined $self->dcpos->{$refpos}) {
        my @DCreads = @{$self->dcpos->{$refpos}};
        $DCreadsCount = $#DCreads + 1;
        if ($DCreadsCount >= $Min_Reads) {
            for my $read (@DCreads){
	            my $base   = $read->{base} || $refbase;
	            my $DCtype = $read->{type};
	            $base = '-' if $DCtype=~/[DI]/;
	            $Allele_info{$DCtype}{$base}{$read->{aln_orient}}{SumQual} += $read->{qual};
	            if (!defined $Allele_info{$DCtype}{$base}{$read->{aln_orient}}{MaxQual} ||
	                $Allele_info{$DCtype}{$base}{$read->{aln_orient}}{MaxQual} < $read->{qual}
                ) {
	                $Allele_info{$DCtype}{$base}{$read->{aln_orient}}{MaxQual} = $read->{qual};
                }
	            $Allele_info{$DCtype}{$base}{$read->{aln_orient}}{Rpos}{$read->{rpos}}++;
            }
        }
    }

    my %readpos;
    #Map reference positions to unpadded read positions
    for my $read (keys %{$self->_align}) {  #all the reads that are aligned
        my $align = $self->_align->{$read};
        if ($align->{ref_start} <= $refpos && $refpos <= $align->{ref_end}) {
            my $rpos;
            if ($align->{orientation} eq 'U') {  #forward alignment
	            $rpos = $refpos - $align->{ref_start} + $align->{r_start};
            }
            else {  #reverse alignment
	            $rpos = $align->{ref_end} - $refpos + $align->{r_start};
            }
            for my $dc (@{$align->{discrepancy}}) {
	            my ($dctype,$pos) = ($dc=~/(\S+)\s+(\d+)/);
	            next if $pos > $rpos;  #doesn't affect alignment coordinates
	            if ($dctype =~ /D\-*(\d*)/) {
	                $rpos -= $1 || 1;
                }
	            elsif ($dctype =~ /I\-*(\d*)/) {
	                $rpos += $1 || 1;
                }   
	            else{
                }
            }
            $readpos{$read} = $rpos;
        }
    }

    #Get Quality from file
    my $fh = Genome::Sys->open_file_for_reading($f_qual) or return;
    my $line;
    
    do{$line = $fh->getline;} until $line =~ /^>(\S+)\s/ || $fh->eof;
    my $readname = $1;
    
    while (!$fh->eof){
        my $qualstr = '';
        do{
            $line = $fh->getline; 
            chomp $line;
            $qualstr = join ' ', $qualstr, $line if $line !~ /^>(\S+)\s/;
        } until $line =~ /^>(\S+)\s/ || $fh->eof;
        
        if (!$fh->eof) {
            my @quals = split /\s+/, $qualstr;
            if (defined $readpos{$readname}) {
	            my $align = $self->_align->{$readname};
	            my $rpos  = $readpos{$readname};
	            $self->warning_message("rpos=$rpos out of scope $#quals in getAlleleInfo in ParseCrossMatch.pm") 
                    if $rpos < 1 || $rpos > $#quals;
                    
        	    my $qual = $quals[$rpos];
	            if ($qual > $self->min_base_qual) {
	                my $base   = $refbase;
	                my $DCtype = 'W';
	                $Allele_info{$DCtype}{$base}{$align->{orientation}}{SumQual} += $qual;
	                if (!defined $Allele_info{$DCtype}{$base}{$align->{orientation}}{MaxQual} ||
	                    $Allele_info{$DCtype}{$base}{$align->{orientation}}{MaxQual} < $qual
                    ) {
	                    $Allele_info{$DCtype}{$base}{$align->{orientation}}{MaxQual}=$qual;
                    }
	                $Allele_info{$DCtype}{$base}{$align->{orientation}}{Rpos}{$rpos}++;
                }
            }
            ($readname) = ($_ =~ /^\>(\S+)\s/);
            #print "$readname\n";
        }
    }
    $fh->close;
    return \%Allele_info;
}

1;
