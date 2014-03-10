package Genome::Model::Tools::Simulation::AlterReferenceSequence;

use strict;
use warnings;
use Data::Dumper;
use Genome;           
use Genome::Info::IUB;
our $VERSION = '0.01';

class Genome::Model::Tools::Simulation::AlterReferenceSequence {
    is => 'Command',
    has_optional_input => [
    ref_fasta => {
        type => 'String',
        is_optional => 0,
        is_input=>1,
    },
    mutation_bed => {
        type => 'String',
        is_optional => 0,
        doc => 'list of mutations in bed format',
    },
    output_file => { 
        type => 'String',
        doc => 'output fasta name.' ,
        is_optional=>1,
    },
    region=> {
        type=>'String',
        is_optional =>0,
        doc =>"restrict the output to a subsection of the input ref",
    },
    limit_regions => {
        type =>'String',
        is_optional=>1,
        doc =>'restrict the output by a target bed. may conflict with region. use intelligently',
    },
    ],
};

sub help_brief {
    "Mutates a reference in preparation for read simulation"
}

sub help_detail {
}

sub execute {
    my $self=shift;
    $DB::single=1;

    unless (-s $self->mutation_bed) {
        $self->error_message("mutation_list has no size?");
        return 0;
    }
    my $output_file_name = $self->mutation_bed;
    $output_file_name =~ s/\.bed/\.fasta/g;
    $self->output_file($output_file_name);

    my @muts=$self->read_mutation_list($self->mutation_bed);
    unless(@muts) {
        $self->error_message("Able to read mutation list, but returned no valid mutations.");
        return 0;
    }
    my $region = $self->region;

    my ($refseq, $offset)=$self->read_region_of_fasta($self->region, $self->ref_fasta);
    my $refseq_length= length($refseq);
    my ($out1, $out1_name) =Genome::Sys->create_temp_file();
    #my $out1=IO::File->new($self->output_file, ">");
    my ($out2, $out2_name) = Genome::Sys->create_temp_file();
#   my $out2=IO::File->new($self->output_file . "2", ">");
    my $pos = $offset;
    my $newrefseq1='';
    my $newrefseq2='';
    foreach my $mut(@muts){
        print "applying";
        print Dumper $mut;
        my $start=$pos-$offset;
        my $len=$mut->{start}-$pos;
        if($mut->{type}=~/SNP/i){
            $len++;
        }
        #my $subseq=unpack("x$start a$len", $refseq);
        my $subseq=substr($refseq,$start,$len);
        next if(length($subseq)!=$len);

        $newrefseq1.=$subseq;
        $newrefseq2.=$subseq;
        $pos+=$len;
#       if($mut->{type}=~/DEL/i){
#           $pos=$mut->{stop};
#       }
#       elsif($mut->{type}=~/INS/i){
#           $newrefseq.=$mut->{variant};
#       }
        if($mut->{type}=~/SNP/i){
            $newrefseq1.=$mut->{variant};
            $newrefseq2.=$mut->{reference};
            $pos++;
        }
        else{

            printf "%s\t%d\t%d\t%d\t%s\n",$mut->{chr},$mut->{start},$mut->{end},$mut->{size},$mut->{type};
            die;
        }
#        $newrefseq1=$self->print_and_flush($newrefseq1, $out1,0);
#        $newrefseq2=$self->print_and_flush($newrefseq2, $out2,0);
    }
    my $start=$pos-$offset;
    my $len=length($refseq)-$start;
    my $subseq=substr($refseq,$start,$len);
    $newrefseq1.=$subseq;
    $newrefseq2.=$subseq;
    if($self->limit_regions) {
        $self->dump_some_regions($newrefseq1, $out1, $offset, $self->limit_regions, "A");
        $self->dump_some_regions($newrefseq2, $out2, $offset, $self->limit_regions, "B");
    }
    else {
        my $desc="mutated according to " . $self->mutation_bed;
        $out1->print(">$region-A\t$refseq_length\t$desc\n");
        $out2->print(">$region-B\t$refseq_length\t$desc\n");
        $self->print_and_flush($newrefseq1, $out1);
        $self->print_and_flush($newrefseq2, $out2);
    }

    $out1->close;
    $out2->close;
    my $final_output = $self->output_file;
    if(Genome::Sys->shellcmd(cmd=>"cat $out1_name $out2_name > $final_output")) {
        return 1;
    }
    else {
        return 0;
    }

}

sub write_region_limited_bed_file {
    my $self = shift;
    my $target_regions_padded_bed = shift;
    my $mutation_bed = $self->mutation_bed;
    my ($basename,$dirname,$suffix) = File::Basename::fileparse($mutation_bed,qw/.bed/);
    my $output_bed = "$dirname/$basename.region_limited.bed";
    my $cmd = "intersectBed -wa -u -a $mutation_bed -b $target_regions_padded_bed > $output_bed";
    Genome::Sys->shellcmd(cmd=>$cmd);
}
    



sub dump_some_regions {
    my($self, $ref, $out_fh, $offset, $bed_file, $hap) = @_;
    my $bed_fh = IO::File->new($bed_file);
    my ($temp_bed, $temp_bed_path) = Genome::Sys->create_temp_file();
    $DB::single=1;
    my $merged_padded_bed = Genome::Sys->create_temp_file_path();
    while(my $line = $bed_fh->getline) {
        chomp($line);
        my ($chr, $start, $stop, undef) = split /\t/, $line;
        if($self->region && $self->in_region($chr, $start, $stop)) {
            $start-=600; 
            $stop +=600;
            $temp_bed->print("$chr\t$start\t$stop\n");
        }
    }
    $temp_bed->close;
    my $cmd = "mergeBed -i $temp_bed_path > $merged_padded_bed";
    Genome::Sys->shellcmd(cmd=>$cmd);
    $bed_fh->close;
    $bed_fh = IO::File->new($merged_padded_bed);
    $self->write_region_limited_bed_file($merged_padded_bed); #for analysis convenience
    while(my $line = $bed_fh->getline) {
        chomp($line);
        my ($chr, $bed_start, $bed_stop, undef) = split /\t/, $line;
        my $stop = $bed_stop - $offset;
        my $start = $bed_start - $offset;
        if($stop > length($ref)) {
            $stop = length($ref);
        }
        if($start < 0) {
            $start =0;
        }
        my $length = $stop - $start;
        if($self->region) {
            if($self->in_region($chr,$bed_start, $bed_stop)) {
                $out_fh->print(">$chr:$bed_start-$bed_stop:$hap\t$length\n"); 
                my $fasta_seq = substr($ref, $start, $length);
                $self->print_and_flush($fasta_seq, $out_fh);
                $self->debug_message("Dumped 22:$bed_start-$bed_stop to fasta\n");
            }
        }
    }
}

sub in_region {
    my ($self, $chr,$start, $stop) = @_;
    my $region = $self->region;
    my ($reg_chr, $start_stop) = split ":", $region;
    my ($reg_start, $reg_stop)= split "-", $start_stop;
    if($chr eq $reg_chr) {
        if(($reg_start <= $start) && ($stop <= $reg_stop)) {
            return 1;
        }
    }
    return 0;
}




sub print_and_flush {
    my $self = shift;
    my $refseq = shift;
    my $out_fh=shift;
    my $chomp=0;
    if(length($refseq) % 60==0) {
        $chomp=1;
    }
    $refseq =~ s/(.{60})/$1\n/g;
    if($chomp) {
        chomp($refseq);
    }
    $out_fh->print($refseq . "\n");
}


sub read_region_of_fasta {
    my $self = shift;
    my $region = shift;
    my $fasta = shift;
    my ($chr, $start_stop) = split ":", $region;
    my $cmd = "samtools faidx $fasta $region";  #FIXME use the stupid auto samtools version decider at some point
    $self->debug_message("ref command: $cmd");
    chomp(my @ref_string = `$cmd`);    
    shift @ref_string;
    my $return_string = join("", @ref_string);
    if($start_stop) {
        my ($start, $stop) = split "-", $start_stop;
        return($return_string, $start);
    }
    return ($return_string, 0);
}

sub read_mutation_list{
    my ($self, $mutlist) = @_;
    my $mut_fh = IO::File->new($mutlist);
    my @muts;
    while(my $line = $mut_fh->getline){
        chomp($line);
        my $mut;
        my ($chr,$start,$stop,$ref_var)=split /\s+/, $line;
        my ($ref, $var) = split "/", $ref_var;
        if($self->region && $self->in_region($chr,$start, $stop)) {
            $mut->{chr}=$chr;
            $mut->{start}=$start;
            $mut->{stop}=$stop;
            $mut->{reference}=$ref;
            $mut->{variant}=$var;
            my $type = $self->infer_variant_type($mut);
            $mut->{type}=$type;

#    print STDERR "$_\n";
            push @muts,$mut;
        }
    }
    return @muts;
}

sub infer_variant_type {
    my ($self,$variant) = @_;

# If the start and stop are the same, and ref and variant are defined its a SNP
    if (($variant->{stop} == $variant->{start}+1)&&
        ($variant->{reference} ne '-')&&($variant->{reference} ne '0')&&
        ($variant->{variant} ne '-')&&($variant->{variant} ne '0')) {
        return 'SNP';
# If start and stop are 1 off, and ref and variant are defined its a DNP
    } elsif (($variant->{reference} eq '-')||($variant->{reference} eq '0')) {
        return 'INS';
    } elsif (($variant->{variant} eq '-')||($variant->{variant} eq '0')) {
        return 'DEL';
    } else {
        die("Could not determine variant type from variant:");
    }
}


1;
