
# review jlolofie
#
# 1. it was claimed that this tool doesnt correctly handle the case when
#    chromosomes are out of order; as a result, a new tool was created:
#    G:M:T:Snp:IntersectChromPos; Ideally we'd only have one tool that
#    does this stuff; Both are used elsewhere in the codebase

package Genome::Model::Tools::Snp::Intersect;

use IO::File;
use strict;
use warnings;
use UR;

class Genome::Model::Tools::Snp::Intersect {
    is => 'Command',
    has => [
        file1               => { is => 'FileName', shell_args_position => 1,
                                    doc => 'first file to intersect' },
        file2               => { is => 'FileName', shell_args_position => 2, is_optional => 1,
                                    doc => 'second file to intersect' },
        intersect_output    => { is => 'FileName', is_optional => 1, 
                                    doc => 'instead of stdout, direct the intersection to this file' },
        f1_only_output      => { is => 'FileName', is_optional => 1, 
                                    doc => 'items present only in the first input should be dumped here' },
        f2_only_output      => { is => 'FileName', is_optional => 1, 
                                    doc => 'items present only in the second input should be dumped here' },
        format              => { is => 'Text', is_optional => 1, default_value => 'default', 
                                    doc => 'combine: show content of both, compare: show both plus snp comparison details.' },
        delimiter1          => { is => 'Text', default_value => '\s+', },
        delimiter2          => { is => 'Text', default_value => '\s+', },
        headers1            => { is => 'Integer', default_value => 0, is_optional=>1, doc => 'file 1 has n header lines' },
        headers2            => { is => 'Integer', default_value => 0, is_optional=>1, doc => 'file 2 has n header lines' },
        _only_in_f1         => { is => 'Integer', default_value => 0, is_optional=>1, doc => 'number only in file 1'},
        _only_in_f2         => { is => 'Integer', default_value => 0, is_optional=>1, doc => 'number only in file 2'},
        _intersection       => { is => 'Integer', default_value => 0, is_optional=>1, doc => 'number intersected'},
    ],
    doc => "intersect two snp lists by position, with optional genotype overlap detail"
};

sub help_synopsis {
    my $self = shift;
    return qq/gmt snp intersect list1.snps list2.snps 

gmt snp intersect list1.snps list2.snps -i intersect.out -f1 f1.out -f2 f2.out

maq cns2view 1.cns | gmt snp intersect mypositions  | less
/
}

sub help_detail{
    return <<EOS;
Intersect two SNP lists and produce a report on the intersection.

The lists should be sorted by chromosome, then position.  Numeric chromosomes should come first, in numeric order, followed by alphabetic chromosomes in alpha order.  Positions should be in numeric order.

The counts of the intersection and the f1/f2 only positions are displayed to STDERR after execution.
You can re-create this output with a simple shell command if f1_only_output and f2_only_output are specified.
    wc -l intersection.out f1.only f2.only

By default the intersection file contains all data from file 1 at the intersected positions.  The detail output shows data from both files.

Detail Output:

If the "detail" option is specified, the intersection file will contain data from both files, and an additional field which describes the genotype difference match/miss, and the transition from homozygous to heterozygous.

This tool expects the 3rd and 4th columns of the inputs to be the reference call and genotype call (in IUB format).  This is the same format maq uses for cns2snp and cns2view output.

The detail output will show: chromosome, position, ref genotype, f1 genotype, f2 genotype, comparison tag.
After these columsn all of the f1 remaining fields will appear, then all of the f2 fields, with a ':' field between.

The comparison tag follows the pattern (match|mismach)-(het|hom)-(het|hom)-\\d+-base-overlap.  After running a detail

You can re-create the detail summary by comparison with a shell command later w/o re-running:
    cat intersection.out | columns 5 | sort |  uniq -c 

EOS
}

my %iub = (
    qw/
        A   A
        C   C
        G   G
        T   T
        R   AG  
        Y   CT  
        K   GT  
        M   AC  
        S   GC  
        W   AT  
        B   CGT 
        D   AGT 
        H   ACT 
        V   ACG 
        N   AGCT 
    /
);

my %iub_overlap;
for my $c1 (sort keys %iub) {
    for my $c2 (sort keys %iub) {
        my $n1 = $iub{$c1};
        my $n2 = $iub{$c2};
        my $c = 0;
        for $b (qw/A C G T/) {
            if ((index($n1,$b) != -1) and (index($n2,$b) != -1)) {
                $c++;
            }
        }
        $iub_overlap{$c1}{$c2} = $c;
    }
}

sub execute {
    my $self = shift;

    # setting a terrible example by using 2 and 3 letter variable names...   
    
    my $f1 = $self->file1;
    
    my $f2 = $self->file2;
    unless (defined $f2) {
        # when only one file is specified, STDIN becomes the "1st" file
        $f2 = $f1;
        $f1 = '-';
    }

    unless ($f1 eq '-' or -e $f1) {
        $self->error_message("File not found: $f1");
        return;
    }
    my $h1 = ($f1 eq '-' ? 'STDIN' : IO::File->new($f1));
    unless ($h1) {
        $self->error_message("Failed to open file $f1: $!");
        return;
    } 

    unless (-e $f2) {
        $self->error_message("File not found: $f2");
        return;
    }
    my $h2 = IO::File->new($f2);
    unless ($h2) {
        $self->error_message("Failed to open file $f2: $!");
        return;
    } 

    my $fi  = $self->intersect_output; 
    my $xi;
    if (my $fi = $self->intersect_output) {
        unless ($fi eq '/dev/null') {
            $xi = IO::File->new(">$fi") or die "Failed to open $fi: $!\n";
        }
    }
    else {
        $xi = 'STDOUT';
    }

    my ($x1,$x2);
    my $f1o = $self->f1_only_output;
    if ($f1o) {
        $x1 = IO::File->new(">$f1o");
        unless ($x1) {
            $self->error_message("Failed to open file $x1: $!");
            return;
        } 
    }

    my $f2o = $self->f2_only_output; 
    if ($f2o) {
        $x2 = IO::File->new(">$f2o");
        unless ($x2) {
            $self->error_message("Failed to open file $x2: $!");
            return;
        } 
    }

    my $n1 = 0;
    my $n2 = 0;
    my $ni = 0;    
    my %intersect_groups;

    my $delimiter1 = $self->delimiter1;
    my $delimiter2 = $self->delimiter2;

    no warnings;
    my ($v1,$c1,$p1,$r1,$g1,@t1);
    my $getf1 = sub {
        $v1 = $h1->getline;
        chomp($v1);
        if (defined $v1) {
            ($c1,$p1,$r1,$g1,@t1) = split(/$delimiter1/,$v1);
            #print "f1: $c1 : $p1 : $g1\n";
        }
        else {
            $c1 = undef;
        }
    };

    my ($v2,$c2,$p2,$r2,$g2,@t2);
    my $getf2 = sub {
        $v2 = $h2->getline;
        chomp($v2);
        if (defined $v2) {
            ($c2,$p2,$r2,$g2,@t2) = split(/$delimiter2/,$v2);
            #print "f2: $c2 : $p2 : $g2\n";
        }
        else {
            $c2 = undef;
        }
    };
    use warnings;

    my $format = $self->format || 'default';
    my $printer = sub {
        no warnings;
        my ($h,$c1,$p1,$r1,$g1,$t1,$r2,$g2,$t2)=@_;
        my $g1_het = ($g1 eq $iub{$g1} ? 'hom' : 'het');
        if ($g2) {
            if ($format eq 'compare') {
                my $g2_het = ($g2 eq $iub{$g2} ? 'hom' : 'het');
                my $m = $iub_overlap{$g1}{$g2};
                unless (defined $m) {
                    print "no value for >$g1< >$g2<\n";
                }
                my $desc = ($g1 eq $g2 ? 'match' : 'miss').'-'.$g1_het.'-'.$g2_het.'-'.$m.'-base-overlap';
                $intersect_groups{$desc}++; 
                $h->print(join("\t",$c1,$p1,$r1,$g1,$g2,$desc,@$t1,":",@$t2),"\n");
            }
            elsif ($format eq 'combine') {
                $h->print(join("\t",$c1,$p1,$r1,$g1,@$t1,':',$r2,$g2,@$t2),"\n");
            }
            else {
                $h->print(join("\t",$c1,$p1,$r1,$g1,@$t1),"\n");
            }
        }
        else {
            if ($format eq 'compare') {
                $h->print(join("\t",$c1,$p1,$r1,$g1,$g1_het,@$t1,),"\n");
            } 
            else {
                $h->print(join("\t",$c1,$p1,$r1,$g1,@$t1),"\n");
            }
        }
    };

    my $headers1 = $self->headers1 || 0;
    my $headers2 = $self->headers2 || 0;
    my $headers_count = ($headers1 >= $headers2 ? $headers1 : $headers2);
    while ($headers_count) {
        if ($headers1) {
            $getf1->();
            $headers1--;
            #print "got1\n";
        }
        else {
            $c1 = $p1 = $r1 = $g1 = '';
            @t1 = ();
        }
        if ($headers2) {
            $getf2->();
            $headers2--;
            #print "got2\n";
        }
        else {
            $c2 = $p2 = $r2 = $g2 = '';
            @t2 = ();
        }
        $printer->($xi,$c1,$p1,$r1,$g1,\@t1,$r2,$g2,\@t2) if $xi;
        $printer->($x1,$c1,$p1,$r1,$g1,\@t1) if $x1;
        $printer->($x2,$c2,$p2,$r2,$g2,\@t2) if $x2;
        $headers_count--;
    }

    $getf1->();
    $getf2->();
    while ($v1 or $v2) {
        # compare chromosomes
        no warnings;
        my $cc;
        if ((not defined $c1) or ($c1 == 0 and $c2 > 0)) {
            $cc = 1;
        }
        elsif ((not defined $c2) or ($c2 == 0 and $c1 > 0)) {
            $cc = -1;
        }
        elsif ($c1 == 0 and $c2 == 0) {
            $cc = ($c1 cmp $c2) 
        }
        else {
            $cc = ($c1 <=> $c2);
        }

        if (($cc == -1) or ($cc == 0 and $p1 < $p2)) {
            $n1++;
            $printer->($x1,$c1,$p1,$r1,$g1,\@t1) if $x1;
            $getf1->();
        }
        elsif ($cc == 1 or ($cc == 0 and $p2 < $p1)) {
            $n2++;
            $printer->($x2,$c2,$p2,$r2,$g2,\@t2) if $x2;
            $getf2->();
        }
        elsif ($cc == 0 and $p1 == $p2) {
            $ni++;
            $printer->($xi,$c1,$p1,$r1,$g1,\@t1,$r2,$g2,\@t2) if $xi;
            $getf1->();
            $getf2->();
        }
        else {
            die "$v1\n$v2\n";
        }
    }

    $h1->close;
    $h2->close;

    $self->_only_in_f1($n1);
    $self->_only_in_f2($n2);
    $self->_intersection($ni);
    print STDERR "$f1 only:\t$n1\n";
    print STDERR "$f2 only:\t$n2\n";
    print STDERR "intersection:\t$ni\n";
    for my $g (sort keys %intersect_groups) {
        print STDERR "\t", $g, ": ", $intersect_groups{$g}, "\n";
    }

    1;
}

sub chr_cmp {
    my ($a,$b) = @_;
    no warnings;
    if ($a == 0 and $b > 0) {
        return 1;
    }
    elsif ($b == 0 and $a > 0) {
        return -1;
    }
    elsif ($a == 0 and $b == 0) {
        return ($a cmp $b) 
    }
    else {
        return ($a <=> $b);
    }
}


1;

