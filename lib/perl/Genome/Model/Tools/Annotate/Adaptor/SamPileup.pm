package Genome::Model::Tools::Annotate::Adaptor::SamPileup;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Annotate::Adaptor::SamPileup {
    is => 'Genome::Model::Tools::Annotate',
    has => [
    pileup_file => {
        is  => 'String',
        is_input  => 1,
        doc => 'The somatic file output from samtools pileup -c to be adapted',
    },
    output_file => {
        is => 'Text',
        is_input => 1,
        is_output => 1,
        doc => "Store output in the specified file instead of sending it to STDOUT."
    },
    ],
};

sub help_brief {
    "Converts samtools pileup -c output into a gmt annotate transcript-variants friendly input format",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt annotate adaptor sam-pileup  ...    
EOS
}

sub help_detail {                           
    return <<EOS 
Converts samtools pileup -c output into a gmt annotate transcript-variants friendly input format
EOS
}

sub execute {
    my $self = shift;

    unless (-s $self->pileup_file) {
        $self->error_message("pileup file must exist");
        die;
    }
    my $pileup_fh = IO::File->new($self->pileup_file);

    # establish the output handle for the transcript variants
    my $output_fh;
    if (my $output_file = $self->output_file) {
        $output_fh = $self->_create_file($output_file);
    }
    else {
        $output_fh = 'STDOUT';
    }

    my @return;   
    while (my $line=$pileup_fh->getline) {
        chomp $line;
        my @fields = split("\t", $line);
        
        if($fields[2] eq '*') {
            #INDELOL!!!
            @return = $self->parse_indel_line($line);
            for my $ret (@return) {
                $output_fh->print(join("\t",@{$ret}) . "\n");
            }
        }
        else { #SNP or Reference. We only care about SNPs
            my ($chr, $start, $ref_base, $variant_iub, $consensus_quality, $snp_quality, $max_map_q, $depth_tumor, $depth_normal) = @fields;
            if($ref_base ne $variant_iub) {
                #It's a SNP
                $output_fh->print("$chr\t$start\t$start\t$ref_base\t$variant_iub\tSNP\t$consensus_quality\t$snp_quality\t$max_map_q\t$depth_tumor\t$depth_normal\n");
            }
        }
    }
}


sub parse_indel_line {
    my $self=shift;
    my ($line) = @_;
    my @return;
    #$self->status_message("(SEARCHING FOR: $line)");

    my ($chr,
        $start_pos,
        $star,
        $consensus,
        @rest_of_fields
    ) = split /\s+/, $line; 
    my ($indel1, $indel2) = split /\//, $consensus;
    my @indels;
    push(@indels, $self->parse_samtools_indel_allele($indel1));
    push(@indels, $self->parse_samtools_indel_allele($indel2));
    for my $indel(@indels) {

        if ($indel->{'sequence'} eq '*') { next; }
        my $hash;
        my $stop_pos;
        my $start;
        if($indel->{'length'} < 0) {
            #it's a deletion!
            $hash->{variation_type}='DEL';
            $start= $start_pos+1;
            $stop_pos = $start_pos + abs($indel->{'length'});
            $hash->{reference}=$indel->{'sequence'};
            $hash->{variant}=0;
        }
        else {
            #it's an insertion
            $hash->{variation_type}='INS';
            $start=$start_pos;
            $stop_pos = $start_pos+1;
            $hash->{reference}=0;
            $hash->{variant}=$indel->{'sequence'};

        }

        $hash->{chromosome}=$chr;
        $hash->{start}=$start;
        $hash->{stop}=$stop_pos;
        #$hash->{num_reads}=$num_reads_across;
        push @return, [@$hash{"chromosome","start","stop","reference","variant","variation_type"}];
    }
    return @return;
}

sub parse_samtools_indel_allele {
    my ($self,$string) = @_;
    my %return;
    if($string eq '*') {
        $return{sequence} = '*';
        $return{length} = 0;
        return \%return;
    }
    my ($sign, $sequence) = $string =~ /([\+\-])([ACGTacgt]+)/;
    unless(defined($sign) && defined($sequence)) {
        $self->error_message("Unexpected indel allele $string");
    }
    $return{sequence} = $sequence;
    if($sign eq '-') {
        $return{length} = - length($sequence);
    }
    else {
        $return{length} = length($sequence);
    }
    return \%return;
}
1;

