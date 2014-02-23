package Genome::Model::Tools::Bsmap::MethRatioWorkflow;

use strict;
use warnings;
use Genome;

my $DEFAULT_VERSION = '2.6';
my $METHRATIO_COMMAND = 'methratio.py';

my %METHRATIO_VERSIONS = (
    '2.6' => '/gscuser/cmiller/usr/src/bsmap-2.6/' . $METHRATIO_COMMAND,
);


class Genome::Model::Tools::Bsmap::MethRatioWorkflow {
    is => 'Command',
    has => [
        bam_file => {
            is => 'Text',
            doc => 'The bam file to do the counting on (must be a product of bsmap alignment',
            is_input => 1,
        },
        output_directory => {
            is => 'Text',
            doc => 'Where to output the methyl counts',
            is_output => 1,
            is_input => 1,
        },
        reference => {
            is => 'Text',
            doc => '36, 37, or a path to the reference fasta',
            is_input => 1,
        },
        version => {
            is => 'Version',
            is_optional => 1,
            is_input => 1,
            default_value => $DEFAULT_VERSION,
            doc => "Version of methratio to use",
        },        
    ],
    has_optional => [
        chromosome => {
            is => 'Text',
            doc => 'process only this chromosome',
            is_input => 1,
        },
        output_zeros => {
            is => 'Boolean',
            doc => 'report loci with zero methylation ratios',
            default => 1,
        },
        no_header => {
            is => 'Boolean',
            doc => 'do not put a header on the file',
            default => 0,
        },

# other options not exposed:
#   -u, --unique          process only unique mappings/pairs.
#   -p, --pair            process only properly paired mappings.
#   -q, --quiet           don't print progress on stderr.
#   -r, --remove-duplicate
#                         remove duplicated reads.
#   -t N, --trim-fillin=N
#                         trim N end-repairing fill-in nucleotides. [default: 2]
#   -g, --combine-CpG     combine CpG methylaion ratios on both strands.
#   -m FOLD, --min-depth=FOLD
#                         report loci with sequencing depth>=FOLD. [default: 1]

    ],
};

sub execute {
    my $self = shift;
    my $fasta;

    if ($self->reference eq "36") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "NCBI-human-build36");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    }
    elsif ($self->reference eq "37") {
        my $reference_build_fasta_object = Genome::Model::Build::ReferenceSequence->get(name => "GRCh37-lite-build37");
        $fasta = $reference_build_fasta_object->cached_full_consensus_path('fa');
    } else { #path to fasta
        if( -s $self->reference){
            $fasta = $self->reference;
        } else {
            $self->error_message('reference must be either "36", "37", or the path to a valid reference file');
            return 0;
        }
    }

    if ($self->chromosome) {
        mkdir $self->output_directory . '/' . $self->chromosome;
        $self->output_directory($self->output_directory . '/' . $self->chromosome);
    }

    if(!(defined($self->version))){
        $self->error_message('methratio version "$self->version" not found.');
        return 0;
    }

    my $cmd = "python " . $METHRATIO_VERSIONS{$self->version};
    $cmd .= " -o ". $self->output_directory ."/snvs.hq";
    $cmd .= " -d " . $fasta;
    if($self->output_zeros){
        $cmd .= " -z";
    }
    if($self->chromosome){
        $cmd .= " -c " . $self->chromosome;
    }
    if($self->no_header){
        $cmd .= " -n";
    }
    $cmd .= " " . $self->bam_file;
    
    $self->debug_message("Running command: $cmd");

    my $return = Genome::Sys->shellcmd(
        cmd => "$cmd",
        );
    unless($return) {
        $self->error_message("Failed to execute: Returned $return");
        die $self->error_message;
    }    

    return 1;
}

1;

sub available_methratio_versions {
    my $self = shift;
    return keys %METHRATIO_VERSIONS;
}
