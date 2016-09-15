package Genome::Model::Tools::Vcf::Convert::Indel::GatkSomaticIndel;

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::Tools::Vcf::Convert::Indel::GatkSomaticIndel {
    is  => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from GATK somatic indel output',
};

our ($N_INDEX, $T_INDEX);


sub help_synopsis {
    <<'HELP';
    Generate a VCF file from Gatk somatic indel output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the indels.
HELP
}

sub source {
    return 'GatkSomaticIndel';
}


#Currently GATK output a off-formatted vcf file, need some
#modification to be standard vcf.
sub initialize_filehandles {
    my $self = shift;

    if ($self->_input_fh || $self->_output_fh) {
        return 1; #Already initialized
    }

    my $input  = $self->input_file;
    my $output = $self->output_file;

    #GATK specific VCF output. For now modify this file. Hacky way
    my $dir     = dirname $input;
    my $raw_vcf = $dir . '/gatk_output_file.vcf'; 

    unless (-s $raw_vcf) {
        die $self->error_message("gatk_output_file.vcf is not available under the input directory");
    }

    my $input_fh  = Genome::Sys->open_file_for_reading($raw_vcf)
        or die "Failed to open $raw_vcf for reading\n";
    my $output_fh = Genome::Sys->open_gzip_file_for_writing($output) 
        or die "Failed to open $output for writing\n";
    
    $self->_input_fh($input_fh);
    $self->_output_fh($output_fh);

    return 1;
}


sub parse_line {
    my ($self, $line) = @_;
    
    if ($line =~ /^#CHROM\s/) { #header
        my @headers = split /\s+/, $line;
        my %index = (
            $headers[9]  => 9,
            $headers[10] => 10,
        );
        ($N_INDEX, $T_INDEX) = map{$index{$self->$_}}qw(control_aligned_reads_sample aligned_reads_sample);
    
        unless ($N_INDEX and $T_INDEX) {
            $self->fatal_message("Failed to get correct normal/tumor sample header for: $line");
        }
    }

    return if $line =~ /^#/;  #skip the other vcf headers
    return unless $line =~ /SOMATIC;/; #skip non-somatic events

    my @columns = split /\s+/, $line;
    my @info    = split /;/, $columns[7];

    unless (@info == 15) {
        $self->fatal_message("line: $line got invalid info field");
    }

    my ($n_dp)  = $info[1]  =~ /N_DP=(\S+)$/;
    my ($n_bq)  = $info[4]  =~ /N_NQSBQ=(\S+?),/;
    my ($n_dp4) = $info[6]  =~ /N_SC=(\S+)$/;
    my ($t_dp)  = $info[9]  =~ /T_DP=(\S+)$/;
    my ($t_bq)  = $info[12] =~ /T_NQSBQ=(\S+?),/;
    my ($t_dp4) = $info[14] =~ /T_SC=(\S+)$/;

    $n_dp4 = _rearrange($n_dp4);
    $t_dp4 = _rearrange($t_dp4);

    unless ($n_dp4 and $t_dp4) {
        $self->fatal_message("Failed to get either normal dp4 or tumor dp4 for line: $line");
    }

    my ($t_gt) = $columns[$T_INDEX] =~ /^(\S+):/;
    my ($n_gt) = $columns[$N_INDEX] =~ /^(\S+):/;

    #Now construct new line
    $columns[6]  = 'PASS';
    $columns[7]  = '.';
    $columns[8]  = 'GT:DP:DP4:BQ:SS';
    $columns[9]  = join ':', $n_gt, $n_dp, $n_dp4, int($n_bq), '.';
    $columns[10] = join ':', $t_gt, $t_dp, $t_dp4, int($t_bq), 2;  #this file should contain only somatic indels

    my $new_line = join "\t", @columns;
    return $new_line;
}


sub _rearrange {
    my $dp4 = shift;
    my @ct  = split /,/, $dp4;

    return unless @ct == 4;
    return join ',', $ct[2], $ct[3], $ct[0], $ct[1];
}


1;

