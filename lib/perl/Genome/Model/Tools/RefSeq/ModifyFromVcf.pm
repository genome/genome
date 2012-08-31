package Genome::Model::Tools::RefSeq::ModifyFromVcf;

use warnings;
use strict;
use File::stat;
use Genome::Info::IUB;
use Data::Dumper;
use IO::File;
use Math::Random;

class Genome::Model::Tools::RefSeq::ModifyFromVcf {
    is => 'Command',
    has => [
    reference_fasta => {
        type => 'String',
        is_optional => 0,
        default=>"/gscmnt/sata839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa",
        doc=>'reference seqeunce fasta. Will be used with samtools to create a new reference',
    },
    mutation_list => {
        type => 'String',
        is_optional => 0,
        doc => 'list of mutations in vcf format',
    },
    min_mutation_size => {
        type => 'String',
        doc => "Optionally only mutate if the mutation is bigger than a certain size.",
    },
    fasta_name_string => {
        type => 'String',
        doc => 'String to append to fasta names, if any. example 1 --> 1_somatic if i supply "somatic"',
    },
    haplotype_1_output_file => { 
        type => 'String',
        doc => 'output fasta name',
    },
    haplotype_2_output_file => { 
        type => 'String',
        doc => 'output fasta name',
    },
    ],
};

sub help_brief {
    "Generates two haplotypes from an input sequence and a VCF file, taking phasing into account"
}

sub help_detail {
    "This will read in a VCF file and apply the mutations listed in the first data column to the reference. If phased " .
    "the mutations will be applied to the appropriate haplotype (hap1 for leftmost, hap2 for rightmost) and if unphased will" .
    " be randomly assigned."
}

sub execute {
    my $self=shift;
    my $fasta = $self->fasta_dir;

    unless (-s $self->mutation_list) {
        $self->error_message("mutation_list has no size?");
        return;
    }
    my $fh = IO::File->new($self->mutation_list,"r");
    unless($fh) {
        $self->error_message("Unable to open " . $self->mutation_list);
        return;
    }

    my @out = (IO::File->new($self->haplotype_1_output_file, "w"),IO::File->new($self->haplotype_2_output_file,"w"));

    my $prev_chr = '';
    my @pos = (1,1);
    my @new_refseq = ('','');
    my $refseq;

    while(my $haps = $self->read_mutation_list($fh)) { 
        if($prev_chr ne $haps->[0]->{chr}) {
            #finish up last chromosome
            if($prev_chr) {
                for(my $i = 0; $i < 2; $i++) {
                    my $start=$pos[$i]-1;
                    my $newrefseq = $new_refseq[$i];
                    my $len=length($refseq)-$start;
                    my $subseq=substr($refseq,$start,$len);
                    $newrefseq.=$subseq;
                    my $desc="mutated according to " . $self->mutation_list;

                    my $append_to_fasta_name = '';
                    if($self->fasta_name_string) {
                        $append_to_fasta_name = "_" . $self->fasta_name_string;
                    }
                    $out[$i]->print(">$prev_chr${append_to_fasta_name}_hap$i $desc\n");
                    $newrefseq =~ s/(.{50})/$1\n/g;
                    $out[$i]->print($newrefseq, "\n");
                }
            }
            @pos = (1,1);
            $prev_chr = $haps->[0]->{chr};
            @new_refseq = ('','');
            $refseq = '';
            unless(open(FAIDX,"samtools faidx $fasta $prev_chr |")) {
                $self->error_message("Couldn't pipe samtools faidx");
                return;
            }
            while(my $fasta_line = <FAIDX>) {
                next if $fasta_line =~ /^>/;
                chomp $fasta_line;
                $refseq .= $fasta_line;
            }
            unless(close(FAIDX)) {
                $self->error_message("Error reading from samtools faidx pipe");
                return;
            }
        }

        for(my $i = 0; $i < 2; $i++) {
            my $mut = $haps->[$i];
            print "applying";
            print Dumper $mut;
            my $start = $pos[$i]-1;
            my $len = $mut->{start}-$pos[$i];

            my $subseq=substr($refseq,$start,$len);
            die "Variant appears to be off of the reference sequence\n" if(length($subseq)!=$len);

            $new_refseq[$i].=$subseq;
            $pos[$i] += $len;
            if($mut->{type}=~/DEL/i){
                $pos[$i]=$mut->{stop}+1;
            }
            elsif($mut->{type}=~/INS/i){
                #note that we are storing the variant start for insertions as the base right after the insertion, this makes this construction easier, but possibly confusing
                $new_refseq[$i].=$mut->{variant};
            }
            elsif($mut->{type}=~/SNP/i){
                $new_refseq[$i].=$mut->{variant};   
                $pos[$i]++;
            }
            else{
                printf "%s\t%d\t%d\t%d\t%s\n",$mut->{chr},$mut->{start},$mut->{end},$mut->{size},$mut->{type};
                die;
            }
        } 
    }
    return 1;
}

#This will read in the VCF file and add mutations to the two haplotypes, returning two arrays.
sub read_mutation_list{
    my ($self, $fh) = @_;
    while(<$fh>){
        chomp;
        next if /^##/;  #skip meta information lines
        if(/^#CHROM/) {   #stop at header
            #8 fixed columns then format
            my @header = split "\t";
            unless($header[8] eq "FORMAT" && defined($header[9])) {
                #check that there is a format field and at least one sample. We will use the first sample (column 9)
                die "VCF must contain genotype information\n";
            }
        }
        my @fields = split /\t/;
        my @alleles = ($fields[3],split /,/, $fields[4]);

        my ($genotype) = split /:/, $fields[9]; #genotype is always the first 
        my (@haplotype, $phasing_type);
        ($haplotype[0], $phasing_type, $haplotype[1]) = $genotype =~ /(\d)([\/\|\\])(\d)/;

        next unless(defined$haplotype[0] && defined $haplotype[1]); #should handle the case where we don't have a genotype call ie ".|."
        if($phasing_type eq '/') {
            #genotype is unphased, randomly flip haplotypes
            if(random_uniform() < 0.5) {
                ($haplotype[0], $haplotype[1]) = ($haplotype[1], $haplotype[0]);
            }
        }
        my @mutations;
        for(my $i = 0; $i < 2; $i++) {

            my $mut;
            $mut->{chr} = $fields[0];
            $mut->{start}=$fields[1];
            if($alleles[$haplotype[$i]] =~ /^([ID])(.+)/) {
                if($1 eq 'I') {
                    $mut->{start} = $fields[1] + 1; #altering the "start" coordinate to make assembly of the new reference easier
                    $mut->{variant} = $2;
                    $mut->{type} = 'INS';
                    $mut->{reference} = 0;
                }
                else {
                    $mut->{stop} = $fields[1] + $2 - 1;
                    $mut->{variant} = 0;
                    $mut->{type} = 'DEL';
                    $mut->{reference} = '';
                }
            }
            else {
                $mut->{stop} = $fields[1];
                $mut->{variant} = $alleles[$haplotype[$i]];
                $mut->{type} = 'SNP';
                $mut->{reference} = $alleles[0];
            }

            push @mutations,$mut;
        }
        return \@mutations;
    }
    return; #done so return
}


1;
