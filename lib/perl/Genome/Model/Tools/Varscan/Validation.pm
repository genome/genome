
package Genome::Model::Tools::Varscan::Validation;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Varscan::Somatic    Runs Varscan somatic pipeline on Normal/Tumor BAM files
#
#    AUTHOR:     Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#    CREATED:    12/09/2009 by D.K.
#    MODIFIED:   12/29/2009 by D.K.
#
#    NOTES:
#
#####################################################################################################################################

use strict;
use warnings;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Varscan::Validation {
    is => 'Genome::Model::Tools::Varscan',

    has_input => [                                # specify the command's single-value properties (parameters) <---
        normal_bam       => { is => 'Text', doc => "Path to Normal BAM file", is_optional => 0, },
        tumor_bam        => { is => 'Text', doc => "Path to Tumor BAM file", is_optional => 0, },
        samtools_path	=> { is => 'Text', doc => "Path to SAMtools executable", is_optional => 0, is_input => 1, default => "samtools" },
        output           => { is => 'Text', doc => "Path to Tumor BAM file", is_optional => 1, is_output => 1 },
        output_snp       => { is => 'Text', doc => "Basename for SNP output, eg. varscan_out/varscan.status.snp" , is_optional => 1, is_output => 1},
        output_indel     => { is => 'Text', doc => "Basename for indel output, eg. varscan_out/varscan.status.indel" , is_optional => 1, is_output => 1},
        output_validation=> { is => 'Text', doc => 'Basename for validation output, eg. varscan_out/varscan.status.validation', is_optional => 1, is_output => 1, },
        reference        => { is => 'Text', doc => "Reference FASTA file for BAMs" , is_optional => 1, default_value => (Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa')},
        skip_if_output_present => { is => 'Text', doc => "If set to 1, skip execution if output files exist", is_optional => 1, },
        varscan_params   => { is => 'Text', doc => "Parameters to pass to Varscan" , is_optional => 1, default_value => '--min-var-freq 0.08 --p-value 0.10 --somatic-p-value 0.01 --validation 1 --min-coverage 8'},
        samtools_version => { is => 'Text', doc => 'Version of samtools to use', default=> 'r544' },
    ],

    has_param => [
        lsf_resource => { default_value => 'select[model!=Opteron250 && type==LINUX64 && tmp>1000] rusage[mem=4000,tmp=1000]'},
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Run the Varscan somatic variant detection"                 
}

sub help_synopsis {
    return <<EOS
Runs Varscan from BAM files
EXAMPLE:    gmt varscan somatic --normal-bam [Normal.bam] --tumor-bam [Tumor.bam] --output varscan_out/Patient.status ...
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;

    ## Get required parameters ##
    my $normal_bam = $self->normal_bam;
    my $tumor_bam = $self->tumor_bam;

    ## Get output directive ##
    my $output = my $output_snp = my $output_indel = my $output_validation = "";

    if($self->output) {
        $output = $self->output;
        $output_snp = $output . ".snp";
        $output_indel = $output . ".indel";
        $output_validation = $output . ".validation";
    } elsif($self->output_snp && $self->output_indel) {
        $output_snp = $self->output_snp;
        $output = $output_snp;
        $output_indel = $self->output_indel;
        $output_validation = $self->output_validation || $output_snp . '.validation';
    } else {
        die "Please provide an output basename (--output) or output files for SNPs (--output-snp) and indels (--output-indel)\n";
    }

    my $reference = $self->reference;
    my $varscan_params = $self->varscan_params;

    ## Check skip if output present ##

    if($self->skip_if_output_present) {
        if(-e $output_snp) {
            my $snp_len = `cat $output_snp | wc -l`;
            chomp($snp_len);
            if($snp_len > 1) {
                return 1;
            }
        }
    }

    my $temp_dir = Genome::Sys->create_temp_directory();
    my $temp_snp = join('/', $temp_dir, 'output.snp');
    my $temp_indel = join('/', $temp_dir, 'output.indel');
    my $temp_output = join('/', $temp_dir, 'output');
    my $temp_validation = $temp_output . '.validation'; #generated name in varscan

    if(-e $normal_bam && -e $tumor_bam) {
        ## Prepare pileup commands ##

#        my $sam_pathname = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
#       my $normal_pileup = "$sam_pathname pileup -f $reference $normal_bam";
#        my $tumor_pileup = "$sam_pathname pileup -f $reference $tumor_bam";
        
        my $normal_pileup = $self->samtools_path . " mpileup -f $reference $normal_bam";
        my $tumor_pileup = $self->samtools_path . " mpileup -f $reference $tumor_bam";
        my $varscan_path = Genome::Model::Tools::Varscan->path_for_version($self->version);

        my $cmd = $self->java_command_line("somatic <\($normal_pileup\) <\($tumor_pileup\) $temp_output --output-snp $temp_snp --output-indel $temp_indel $varscan_params");

        Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files => [$normal_bam, $tumor_bam, $reference],
            output_files => [$temp_snp, $temp_indel],
            allow_zero_size_output_files => 1,
        );

        Genome::Sys->copy_file($temp_snp, $output_snp);
        Genome::Sys->copy_file($temp_indel, $output_indel);
        Genome::Sys->copy_file($temp_validation, $output_validation) if Genome::Sys->check_for_path_existence($temp_validation); #optional file
    } else {
        die "Error: One of your BAM files doesn't exist!\n";
    }

    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


1;
