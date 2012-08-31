package Genome::Model::Tools::Snp::IntersectChromPos;

use IO::File;
use strict;
use warnings;
use UR;

class Genome::Model::Tools::Snp::IntersectChromPos {
    is => 'Command',
    has => [
        file1               => { is => 'filename', is_optional=>0,},
        file2               => { is => 'filename', is_optional=>0,},
        intersect_output    => { is => 'FileName', is_optional => 1, 
                                doc => 'instead of stdout, direct the intersection to this file' },
        f1_only_output      => { is => 'FileName', is_optional => 1, 
                                doc => 'items present only in the first input should be dumped here' },
        f2_only_output      => { is => 'FileName', is_optional => 1, 
                                doc => 'items present only in the second input should be dumped here' },
        delimiter1          => { is => 'Text', default_value => '\s+', },
        delimiter2          => { is => 'Text', default_value => '\s+', },
        headers1            => { is => 'Integer', default_value => 0, is_optional=>1, doc => 'file 1 has n header lines' },
        headers2            => { is => 'Integer', default_value => 0, is_optional=>1, doc => 'file 2 has n header lines' },
        consider_genotype   => { is_optional=>1, default=>0, doc=>'thingy for dave larson he knows what it means' },
        ignore_sorting      => { is => 'Boolean', is_optional=>1, default=>0, 
                                 doc=>'Setting this to true causes the tool to ignore any sorting or missing chromosome issues. Use this flag ONLY if you are SURE both files are sorted the same by chromosome and position.' },
    ],
};

sub execute {
    my $self = shift;
    my $file1_fh=IO::File->new($self->file1);
    my $file2_fh=IO::File->new($self->file2);

    unless($file1_fh && $file2_fh) {
        $self->error_message("Could not open file1 or file2");
    }


    my $intersect_fh;
    my $f1_only_fh;
    my $f2_only_fh;

    if($self->f1_only_output) {
        $f1_only_fh=IO::File->new(">" . $self->f1_only_output);
    }
    if($self->f2_only_output) {
        $f2_only_fh=IO::File->new(">" . $self->f2_only_output);
    }
    if($self->intersect_output) {
        $intersect_fh=IO::File->new(">" . $self->intersect_output);
    }
    unless($f1_only_fh && $f2_only_fh && $intersect_fh) {
        $self->error_message("Could not open some output files! $!");
    }     


    $self->warning_message("this tool assumes your files have chromosomes in the same order and the positions are ascending. it should detect any other condition and barf, though."); 
    #both files exist  
    
    my $header_lines_f1 = $self->headers1;
    my $header_lines_f2 = $self->headers2;
    
    my $line2 = $file2_fh->getline;
    for my $skip (1..$header_lines_f2) { 
        $line2 = $file2_fh->getline; 
    }
    my ($chr2, $pos2, $ref1, $genotype1) = split ($self->delimiter2, $line2);
    
    my $line1 = $file1_fh->getline;
    for my $skip (1..$header_lines_f1) { 
        $line1 = $file1_fh->getline; 
    }
    my ($chr1, $pos1, $ref2, $genotype2) = split ($self->delimiter1, $line1);

    my ($prev_chrom1, $prev_chrom2); 
    while(defined $line1 && defined $line2) {
        if($chr1 eq $chr2) {
            if ($pos1 < $pos2) {
                $f1_only_fh->print($line1);
                $line1 = $file1_fh->getline;
                $prev_chrom1=$chr1;
                ($chr1, $pos1) = split ($self->delimiter1, $line1);

            } 
            elsif ($pos1 == $pos2) { 
                if($self->consider_genotype) {
                    if($genotype1 eq $genotype2) {
                        $intersect_fh->print($line1);
                    }
                    else {
                        $f1_only_fh->print($line1);
                        $f2_only_fh->print($line2);
                    }
                }
                else {
                    $intersect_fh->print($line1);
                }
                $line1 = $file1_fh->getline;
                $prev_chrom1=$chr1;
                ($chr1, $pos1) = split ($self->delimiter1, $line1);
                $line2 = $file2_fh->getline;
                $prev_chrom2=$chr2;
                ($chr2, $pos2) = split ($self->delimiter2, $line2);
            }
            elsif ($pos1 > $pos2) { 
                $f2_only_fh->print($line2);
                $line2 = $file2_fh->getline;
                $prev_chrom2=$chr2;
                ($chr2, $pos2) = split ($self->delimiter2, $line2);
            }
        }
        elsif($chr1 ne $chr2) {
            if($chr2 eq $prev_chrom1) {
                #file 2 is lagging
                $f2_only_fh->print($line2);
                $line2 = $file2_fh->getline;
                $prev_chrom2=$chr2;
                ($chr2, $pos2) = split ($self->delimiter2, $line2);
            }
            elsif($chr1 eq $prev_chrom2) {
                #file 1 is lagging
                $f1_only_fh->print($line1);
                $line1 = $file1_fh->getline;
                $prev_chrom1=$chr1;
                ($chr1, $pos1) = split ($self->delimiter1, $line1);
            }
            # If both current chromosomes are NT, or if either previous chromosome is NT, the else will make us die (and shouldnt), so just do a cmp
            # OR if we are ignoring sorting, just do a cmp and continue
            elsif ( ($self->ignore_sorting) || 
                    (($chr1 =~ m/NT_/ && $chr2 =~ m/NT_/)||($prev_chrom1 =~ m/NT_/)||($prev_chrom1 =~ m/NT_/))
                  ) {
                if (($chr1 cmp $chr2) < 0) {
                    $f1_only_fh->print($line1);
                    $line1 = $file1_fh->getline;
                    $prev_chrom1=$chr1;
                    ($chr1, $pos1) = split ($self->delimiter1, $line1)
                }
                if (($chr1 cmp $chr2) > 0) {
                    $f2_only_fh->print($line2);
                    $line2 = $file2_fh->getline;
                    $prev_chrom2=$chr2;
                    ($chr2, $pos2) = split ($self->delimiter2, $line2)
                }
            }
            else {
                # the files are whacked out and not sorted the same?
                $self->error_message("Line 1: $line1");
                $self->error_message("Line 2: $line2");
                $self->error_message("previous chromosome1: $prev_chrom1 previous chromosome2: $prev_chrom2");
                $self->error_message("NO NO! YOU GIVE ME FILES SORTED BY CHROMOSOME, POSITION!");
                die;
            }
        }
    }
    if(!defined $line1) {
        $f2_only_fh->print($file2_fh->getlines);
    }
    elsif(!defined $line2) {
        $f1_only_fh->print($file1_fh->getlines);
    }

}


1;

