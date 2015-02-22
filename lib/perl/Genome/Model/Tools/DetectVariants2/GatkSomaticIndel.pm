package Genome::Model::Tools::DetectVariants2::GatkSomaticIndel;

use strict;
use warnings;

use Cwd;

use Genome;

class Genome::Model::Tools::DetectVariants2::GatkSomaticIndel{
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    has => [
        mb_of_ram => {
            is => 'Text',
            doc => 'Amount of memory to allow GATK to use',
            default => 5000,
        },
    ],
    has_param => [
         lsf_resource => {
             default_value => "-M 8000000 -R 'select[mem>8000] rusage[mem=8000]'",
         },
     ],
};

sub _detect_variants {
    my $self = shift;
    my $refseq = $self->reference_sequence_input;
    $refseq =~ s/\/opt\/fscache//;
    my $gatk_raw_output = $self->_temp_staging_directory."/gatk_output_file";
    my $gatk_somatic_output = $self->_temp_staging_directory."/indels.hq";

    my %parameters = (
        tumor_bam   => $self->aligned_reads_input, 
        normal_bam  => $self->control_aligned_reads_input,
        output_file => $gatk_raw_output,
        mb_of_ram   => $self->mb_of_ram,
        reference   => $refseq,
        version     => $self->version,
    );
    $parameters{gatk_params} = $self->params if $self->params;

    my $gatk_cmd = Genome::Model::Tools::Gatk::SomaticIndel->create(%parameters); 
    unless($gatk_cmd->execute){
        $self->error_message("Failed to run GATK command.");
        die $self->error_message;
    }

    # Gatk outputs different types of variants. We want just the somatic stuff for now. So grep out somatic results.
    my $cmd = "grep SOMATIC $gatk_raw_output > $gatk_somatic_output";
    my $return = Genome::Sys->shellcmd( 
        cmd => $cmd, 
        allow_failed_exit_code => 1
    );
    # Grep will return 1 if nothing was found... and 2 if there was some error. So just make sure the return is 0 or 1.
    unless( $return == 0 || $return == 1){
        die $self->error_message("Could not execute grep to separate germline and somatic calls in gatk");
    }
    unless(-e $self->_temp_staging_directory."/indels.hq"){
        my $filename = $self->_temp_staging_directory."/indels.hq";
        my $output = system("touch $filename");
        if($output){
            $self->error_message("creating an empty indels.hq file failed. output=".$output);
            die $self->error_message;
        }
    }

    return 1;
}

sub has_version {
    my $self    = shift;
    my $version = shift;

    unless (defined $version) {
        $version = $self->version;
    }
    return Genome::Model::Tools::Gatk::SomaticIndel->has_version($version);
}


sub parse_line_for_bed_intersection {
    my ($self, $line) = @_;
    die $self->error_message('No line provided to parse_line_for_bed_intersection for GatkSomaticIndel')
        unless $line;

    my ($chromosome, undef, $stop, $refvar) = split /\t/,  $line;
    my ($sign, $refvar_base) = $refvar =~ /^([\-\+])(\S+)$/;
    die $self->error_message("Invalid indel line: $line") unless $sign and $refvar_base;

    my ($ref, $var);
    if ($sign eq '+') {
        ($ref, $var) = (0, $refvar_base);
    }
    else {
        ($ref, $var) = ($refvar_base, 0);
    }
    
    die $self->error_message("Could not get chromosome, position, reference, variant for line: $line")
        unless defined $chromosome and defined $stop and defined $ref and defined $var;

    return [$chromosome, $stop, $ref, $var];
}

1;
