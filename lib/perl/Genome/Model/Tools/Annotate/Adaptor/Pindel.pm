package Genome::Model::Tools::Annotate::Adaptor::Pindel;

use warnings;
use strict;

use Genome;
use Workflow;

class Genome::Model::Tools::Annotate::Adaptor::Pindel {
    is => ['Command'],
    has => [
        input_file => {
            is  => 'String',
            is_input => '1',
            is_output => '1',
            doc => 'The input file. This must be a file in pindel output format. At present this only works for insertions and deletions, not di.',
        },
        output_file => {
            is  => 'String',
            is_input => '1',
            is_output => '1',
            doc => 'The output file. This will be a file in annotation input format representing the pindel calls from the input file',
        },
        reference_fasta => {
            is => 'String',
            default=> Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa',
        },
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
        },
        # Make workflow choose 64 bit blades, this is needed for samtools faidx
        lsf_resource => {
            is_param => 1,
            default_value => 'rusage[mem=4000] select[type==LINUX64] span[hosts=1] -M 4000000',
        },
        lsf_queue => {
            is_param => 1,
            default_value => 'long'
        }, 
       ],
};

sub help_brief {
    "Transforms a file from pindel output format to annotation input format in preparation for transcript annotation";
}

sub help_synopsis {
    return <<"EOS"
gmt annotate adaptor pindel --input-file pindel.outfile --output-file pindel.adapted 
gmt annotate adaptor pindel --in pindel.outfile --out pindel.adapted 
EOS
}

sub help_detail {                           
    return <<EOS 
Transforms a file from pindel output format to annotation input format in preparation for transcript annotation
EOS
}

sub execute {
    my $self = shift;

    # test architecture to make sure we can run (needed for samtools faidx)
    unless (`uname -a` =~ /x86_64/) {
       $self->error_message("Must run on a 64 bit machine");
       die;
    }
    
    my $output_fh = IO::File->new($self->output_file, ">");
    my $big_fh = IO::File->new($self->output_file . ".big_deletions" , ">");

    unless($output_fh && $big_fh) {
        $self->error_message("Unable to open file handles for output");
        return;
    }

    my $input_fh = IO::File->new($self->input_file);
    unless($input_fh) {
        $self->error_message("Unable to open input fh for " . $self->input_file);
        return;
    }
    while(my $line = $input_fh->getline) {
        next unless($line =~ m/^#+$/);
        my $call = $input_fh->getline;
        my $reference = $input_fh->getline;
        my $first_read = $input_fh->getline;
        $self->parse_and_print($call, $reference, $first_read, $output_fh, $big_fh);
    }
}

sub parse_and_print {
    my $self=shift;
    my ($call, $reference, $first_read, $output_fh, $big_fh) = @_;
    #parse out call bullshit
    chomp $call;
    my @call_fields = split /\s+/, $call;    
    my $type = $call_fields[1];
    my $size = $call_fields[2];
    my $chr = $call_fields[4];
    my $start = $call_fields[6];
    my $stop = $call_fields[7];
    my $support = $call_fields[12];
    my ($ref, $var);
    if($type =~ m/D/) {
        $var =0;
        ###Make pindels coordinates(which seem to be last undeleted base and first undeleted base) 
        ###conform to our annotators requirements

        ###also deletions which don't contain their full sequence should be dumped to separate file
        $start = $start + 1;
        $stop = $stop - 1;
        my $allele_string;
        if(($reference =~ m/</) && ($size < 100)) { #pindel has abridged our allele string. bastard program 
            my $sam_default = Genome::Model::Tools::Sam->path_for_samtools_version;
            my $faidx_cmd = "$sam_default faidx " . $self->reference_fasta . " $chr:$start-$stop";
            my @faidx_return= `$faidx_cmd`;
            shift(@faidx_return);
            chomp @faidx_return;
            $allele_string = join("",@faidx_return);
        }
        else {
            ($allele_string) = ($reference =~ m/^[ACGTN]+([acgt]+)[ACGTN]+$/);
        }
        $ref = $allele_string;
    }
    elsif($type =~ m/I/) {
        $ref=0;
        my ($letters_until_space) =   ($reference =~ m/^([ACGTN]+) /);
        my $offset_into_first_read = length($letters_until_space);
        $var = substr($first_read, $offset_into_first_read, $size);
    }
    if($size < 100) {
        $output_fh->print("$chr\t$start\t$stop\t$ref\t$var\t$support\n");
    }
    else {
        $big_fh->print("$chr\t$start\t$stop\t$size\t$support\n");
    }
}









    
    
