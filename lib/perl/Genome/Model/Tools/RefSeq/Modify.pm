package Genome::Model::Tools::RefSeq::Modify;

use warnings;
use strict;
use File::stat;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Genome::Info::IUB;
use Data::Dumper;
use IO::File;




class Genome::Model::Tools::RefSeq::Modify {
    is => 'Command',
    has => [
    fasta_dir => {
        type => 'String',
        is_optional => 0,
        default=> Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/',
        doc=>'right now this tool only works with directories that BIO DB FASTA works with',
    },
    mutation_list => {
        type => 'String',
        is_optional => 0,
        doc => 'list of mutations in annotation input format',
    },
    min_mutation_size => {
        type => 'String',
        doc => "Optionally only mutate if the mutation is bigger than a certain size.",
    },
    fasta_name_string => {
        type => 'String',
        doc => 'String to append to fasta names, if any. example 1 --> 1_somatic if i supply "somatic"',
    },
    output_file => { 
    type => 'String',
    doc => 'output fasta name',
},
    ],
};

sub help_brief {
    "Prepares a fasta file to be used as a new refseq in processing profiles"
}

sub help_detail {
    "Copies a fasta file out to the reference path, and then schedules jobs which will " .
    "create appropriate BWA, Maq, and Samtools index files."
}

sub execute {
    my $self=shift;
    my $fasta = $self->fasta_dir;
    my $db=Bio::DB::Fasta->new($fasta);

    unless (-s $self->mutation_list) {
        $self->error_message("mutation_list has no size?");
        return 0;
    }

    my @muts=$self->read_mutation_list($self->mutation_list);
    unless(@muts) {
        $self->error_message("Able to read mutation list, but returned no valid mutations.");
        return 0;
    }


    my $out=IO::File->new($self->output_file, ">");
    for my $x(1..22,'X','Y') {
        my $refseq=$db->seq($x);
        my $newrefseq='';
        my $pos=1;
        foreach my $mut(@muts){
            if($x eq $mut->{chr}) {
                print "applying";
                print Dumper $mut;
                my $start=$pos-1;
                my $len=$mut->{start}-$pos;

                #my $subseq=unpack("x$start a$len", $refseq);
                my $subseq=substr($refseq,$start,$len);
                next if(length($subseq)!=$len);

                $newrefseq.=$subseq;
                $pos+=$len;
                if($mut->{type}=~/DEL/i){
                    $pos=$mut->{stop}+1;
                }
                elsif($mut->{type}=~/INS/i){
                    $newrefseq.=$mut->{variant};
                }
                elsif($mut->{type}=~/SNP/i){
                    $newrefseq.=$mut->{variant};   
                    $pos++;
                }
                else{

                    printf "%s\t%d\t%d\t%d\t%s\n",$mut->{chr},$mut->{start},$mut->{end},$mut->{size},$mut->{type};
                    die;
                }
            }
        }
        my $start=$pos-1;
        my $len=length($refseq)-$start;
#my $subseq=unpack("x$start a$len", $refseq);
my $subseq=substr($refseq,$start,$len);
$newrefseq.=$subseq;
my $desc="mutated according to " . $self->mutation_list;

my $append_to_fasta_name = '';
if($self->fasta_name_string) {
    $append_to_fasta_name = "_" . $self->fasta_name_string;
}
$out->print(">$x$append_to_fasta_name $desc\n");
$newrefseq =~ s/(.{50})/$1\n/g;
$out->print($newrefseq, "\n");
}
}


sub read_mutation_list{
    my ($self, $mutlist) = @_;
    open(MUT,"<$mutlist") || die "unable to open $mutlist\n";
    my @muts;
    while(<MUT>){
        chomp;
        my $mut;
        my ($chr,$start,$stop,$ref,$var)=split /\s+/;
        $mut->{chr}=$chr;
        $mut->{start}=$start;
        $mut->{stop}=$stop;
        $mut->{reference}=$ref;
        if(length($var)==1 &&  $var ne '0' && $ref ne '0') {
            my ($variant_allele) = Genome::Info::IUB->variant_alleles_for_iub($mut->{reference}, $var);
            $mut->{variant}=$variant_allele;
        }      
        else {
            $mut->{variant}=$var;
        }  
        my $type = $self->infer_variant_type($mut);
        $mut->{type}=$type;
        if($type eq 'INS') {
            $mut->{start} += 1;
        }
#    print STDERR "$_\n";
push @muts,$mut;
}
return @muts;
}

sub infer_variant_type {
    my ($self,$variant) = @_;

    # If the start and stop are the same, and ref and variant are defined its a SNP
    if (($variant->{stop} == $variant->{start})&&
    ($variant->{reference} ne '-')&&($variant->{reference} ne '0')&&
    ($variant->{variant} ne '-')&&($variant->{variant} ne '0')) {
        return 'SNP';
        # If start and stop are 1 off, and ref and variant are defined its a DNP
    } elsif (($variant->{stop} - $variant->{start} == 1)&&
    ($variant->{reference} ne '-')&&($variant->{reference} ne '0')&&
    ($variant->{variant} ne '-')&&($variant->{variant} ne '0')) {
        return 'DNP';
        # If reference is a dash, we have an insertion
    } elsif (($variant->{reference} eq '-')||($variant->{reference} eq '0')) {
        return 'INS';
    } elsif (($variant->{variant} eq '-')||($variant->{variant} eq '0')) {
        return 'DEL';
    } else {
        die("Could not determine variant type from variant:");
    }
}


1;
