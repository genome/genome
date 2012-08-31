package Genome::Model::Tools::DetectVariants2::GatkGermlineIndelUnifiedGenotyper;

use strict;
use warnings;

use Cwd;

use Genome;

class Genome::Model::Tools::DetectVariants2::GatkGermlineIndelUnifiedGenotyper{
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    has_constant => [
        detect_snvs => {},
        detect_svs => {},
        detect_indels => { value => 1 },
    ],
    has => [
        mb_of_ram => {
            is => 'Text',
            doc => 'Amount of memory to allow GATK to use',
            default => 5000,
        },
    ],
    has_param => [
         lsf_resource => {
             default_value => "-R 'rusage[mem=6000] select[type==LINUX64 && model != Opteron250 && mem>6000 && maxtmp>100000] span[hosts=1]' -M 6000000",
         },
     ],
};

sub _detect_variants {
    my $self = shift;
    my $refseq = $self->reference_sequence_input;
    $refseq =~ s/\/opt\/fscache//;
    my $gatk_cmd = Genome::Model::Tools::Gatk::GermlineIndelUnifiedGenotyper->create( 
        bam_file => $self->aligned_reads_input, 
        vcf_output_file => $self->_temp_staging_directory."/indels.hq",
        mb_of_ram => $self->mb_of_ram,
        reference_fasta => $refseq,
        version => $self->version,
    );
    unless($gatk_cmd->execute){
        $self->error_message("Failed to run GATK command.");
        die $self->error_message;
    }
    unless(-s $self->_temp_staging_directory."/indels.hq"){
        my $filename = $self->_temp_staging_directory."/indels.hq";
        Genome::Sys->write_file($filename, '');
    }
    return 1;
}

sub has_version {
    my $self = shift;

    return Genome::Model::Tools::Gatk->has_version(@_);
}

sub parse_line_for_bed_intersection {
    my $class = shift;
    my $line = shift;

    unless ($line) {
        die $class->error_message("No line provided to parse_line_for_bed_intersection");
    }

    # Skip header lines
    if ($line =~ /^#/) {
        return;
    }

    my ($chromosome, $start, $id, $reference, $variant) = split("\t", $line);
    my $stop;
    if(length($reference) == 1 and length($variant) == 1) {
        #SNV case
        $stop = $start;
        $start -= 1; #convert to 0-based coordinate
    } elsif (length($reference) == 1 and length($variant) > 1) {
        #insertion case
        $stop = $start; #VCF uses 1-based position of base before the insertion (which is the same as 0-based position of first inserted base), insertions have no length
        $reference = '*';
        $variant = substr($variant, 1);
    } elsif (length($reference) > 1 and length($variant) == 1) {
        #deletion case
        $reference = substr($reference, 1);
        $stop = $start + length($reference);
        $variant = '*';
    } else {
        die $class->error_message("Unhandled variant type encountered for line: $line");
    }


    unless (defined $chromosome && defined $stop && defined $reference && defined $variant) {
        die $class->error_message("Could not get chromosome, position, reference, or variant for line: $line");
    }

    return [$chromosome, $stop, $reference, $variant];
}

1;
