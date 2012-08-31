package Genome::Model::Tools::DetectVariants2::GatkGermlineIndel;

use strict;
use warnings;

use Cwd;

use Genome;

class Genome::Model::Tools::DetectVariants2::GatkGermlineIndel{
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
             default_value => "-M 8000000 -R 'select[type==LINUX64 && mem>8000] rusage[mem=8000]'",
         },
     ],
};

sub _detect_variants {
    my $self = shift;
    my $refseq = $self->reference_sequence_input;
    $refseq =~ s/\/opt\/fscache//;
    my $gatk_cmd = Genome::Model::Tools::Gatk::GermlineIndel->create( 
        bam_file => $self->aligned_reads_input, 
        output_file => $self->_temp_staging_directory."/gatk_output_file",
        bed_output_file => $self->_temp_staging_directory."/indels.hq",
        mb_of_ram => $self->mb_of_ram,
        reference => $refseq,
        gatk_params => ' -T IndelGenotyperV2 --window_size 300 -et NO_ET',
        version => $self->version,
    );
    unless($gatk_cmd->execute){
        $self->error_message("Failed to run GATK command.");
        die $self->error_message;
    }
    unless(-s $self->_temp_staging_directory."/indels.hq"){
        my $filename = $self->_temp_staging_directory."/indels.hq";
        my $output = `touch $filename`;
        unless($output){
            $self->error_message("creating an empty indels.hq file failed.");
            die $self->error_message;
        }
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

    my ($chromosome, $start, $stop, $refvar) = split "\t",  $line;

    my ($reference, $variant, $type);
    ($refvar) = split ":", $refvar;
    if ($refvar =~ m/\-/) {
        $type = '-';
    }
    else {
        $type = '+';
    }
    $refvar =~ s/[\+,\-]//;
    if ($type eq '+') {
        $reference = '0';
        $variant = $refvar;
    }
    else {
        $variant = '0';
        $reference = $refvar;
    }

    unless (defined $chromosome && defined $stop && defined $reference && defined $variant) {
        die $class->error_message("Could not get chromosome, position, reference, or variant for line: $line");
    }

    return [$chromosome, $stop, $reference, $variant];
}
1;
