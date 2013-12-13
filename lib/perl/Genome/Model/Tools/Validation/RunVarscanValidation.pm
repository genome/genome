package Genome::Model::Tools::Validation::RunVarscanValidation;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;

class Genome::Model::Tools::Validation::RunVarscanValidation {
    is => 'Command',
    has => [
    output_dir => { 
        type => 'String',
        is_optional => 0,
        doc => 'The output directory in which to place the results',
    },

    normal_build => {
        type => 'String',
        is_optional => 1,
        doc => 'The build for the normal model to use. Used to retrieve the normal bam file',
    },

    tumor_build => {
        type => 'String',
        is_optional => 1,
        doc => 'The build for the tumor model to use. Used to retrieve the tumor bam file',
    },

    validation_build => {
        type => 'String',
        is_optional => 1,
        doc => 'The build for the somatic-validation model, which has the tumor and normal bams',
    },
        
    normal_purity => {
        type => 'Float',
        is_optional => 1,
        default => 1,
        doc => 'The purity of the matched normal of the sample. 1 by default. If you have an impure normal sample, then you will want to reduce this value (eg for 10% contamination of the normal use 0.9)',
    },

    min_var_freq => {
        type => 'Float',
        is_optional => 1,
        default => 0.08,
        doc => 'The minimum variant frequency that Varscan uses to determine if a variant exists. 0.08 by default. If you have an impure normal smaple, then you will want to increase this value (eg for 30% contamination of the normal use 0.2',
    },

    snp_target_file => {
        type => 'String',
        is_optional => 0,
        doc => 'A tab-separated file containing, minimally, chromosome then position for filtering the Varscan snp output to targeted positions',
    },
    ]
};


sub execute {
    my $self=shift;

    #we use the builds to grab the bams
    my $tumor_bam;
    my $normal_bam;

    if(defined($self->tumor_build) && defined($self->normal_build)){
        my $tumor_build = Genome::Model::Build->get($self->tumor_build);
        my $normal_build = Genome::Model::Build->get($self->normal_build);

        $tumor_bam = $tumor_build->whole_rmdup_bam_file;
        $normal_bam = $normal_build->whole_rmdup_bam_file;
    } elsif (defined($self->validation_build)){
        my $build = Genome::Model::Build->get($self->validation_build);
        my $data_dir = $build->data_directory;
        $tumor_bam = "$data_dir/alignments/tumor/*.bam";
        $normal_bam = "$data_dir/alignments/normal/*.bam";
    } else {
        print STDERR "Must specify either tumor-build and normal-build parameters OR specify the somatic-validation build";
        return 0;
    }

    

    my $min_var_freq = $self->min_var_freq;
    my $normal_purity = $self->normal_purity;

    my $target_file = $self->snp_target_file;

    my $output_dir = $self->output_dir;
    my $output_file = "$output_dir/varScan.output";
    my $stdout_log = "$output_dir/varScan.out";
    my $stderr_log = "$output_dir/varScan.err";

    my $user = Genome::Sys->username;
    
    #Run varscan in validation mode, accounting for the tumor and normal purities. When this is done, then run gmt snp match to limit to our target regions
    my $command = qq{bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -R 'select[mem>2000 && type==LINUX64] rusage[mem=2000]' -o $stdout_log -e $stderr_log  'gmt varscan validation --samtools-version r544 --output $output_file --normal-bam $normal_bam --tumor-bam $tumor_bam --varscan-params="--min-coverage 8 --min-var-freq $min_var_freq --p-value 0.10 --somatic-p-value 1.0e-02 --validation --normal-purity $normal_purity";gmt snp match -f $target_file $output_file.snp > $output_file.snp.targeted'};
    print `$command`;
        
    return 1;

}


1;

sub help_brief {
    "Runs varscan in validation mode on the tumor and normal bams of a validation run"
}

sub help_detail {
    <<'HELP';
This command is a simple wrapper for running varscan in validation mode with tumor and normal purities. It then filters the result by the passed target file of mutations.
HELP
}
