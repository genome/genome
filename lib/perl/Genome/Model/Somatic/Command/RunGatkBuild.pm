package Genome::Model::Somatic::Command::RunGatkBuild;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;

class Genome::Model::Somatic::Command::RunGatkBuild {
    is => 'Command',
    has => [
    build_id => { 
        type => 'String',
        is_optional => 0,
        doc => "build id of the somatic build to process",
    },
    analysis_dir => { 
        type => 'String',
        is_optional => 0,
        doc => "directory to drop the GATK files into",
    },
    force => {
        type => 'Boolean',
        is_optional => 1,
        default => 0,
        doc => "whether or not to launch blade jobs regardless of symlinks etc." ,
    },
    ]
};


sub execute {
    my $self=shift;
    my $force = $self->force;
    #check the passed directory
    my $dir = $self->analysis_dir;
    $dir=$dir."/";
    unless(-d $dir) {
        $self->error_message("$dir is not a directory");
        return;
    }
    
    $self->error_message( "$dir" );
    #retrieve the build et al
    my $build_id = $self->build_id;

    my $build = Genome::Model::Build->get($build_id);
    unless(defined($build)) {
        $self->error_message("Unable to find build $build_id");
        return unless $force;;
    }
    my $model = $build->model;
    unless(defined($model)) {
        $self->error_message("Somehow this build does not have a model");
        return unless $force;
    }
    unless($model->type_name eq 'somatic') {
        $self->error_message("This build must be a somatic pipeline build");
        return unless $force;
    }

    #retrieve the tumor bam
    my $tumor_bam = $build->tumor_build->whole_rmdup_bam_file;
    unless($tumor_bam) {
        $self->error_message("Couldn't determine tumor bam file from somatic model");
        return unless $force;
    }

    if(-z $tumor_bam) {
        $self->error_message("$tumor_bam is of size 0 or does not exist");
        return unless $force;
    }
    
    my $normal_bam = $build->normal_build->whole_rmdup_bam_file;

    unless($normal_bam) {
        $self->error_message("Couldn't determine normal bam file from somatic model");
        return unless $force;
    }
    if(-z $normal_bam) {
        $self->error_message("$normal_bam is of size 0 or does not exist");
        return unless $force;
    }

    my $genome_name = $build->tumor_build->model->subject->source_common_name;   #this should work
    unless($genome_name) {
        $self->error_message("Unable to retrieve sample name from tumor build");
        return unless $force;
    }
    my $user = Genome::Sys->username; 
    $dir=$dir."$genome_name/gatk/";
    
    if (-d $dir){
    	$self->error_message("$dir already exist, use force to proceed");
	return unless $force;
    }else{
    	mkdir($dir, 0777) || print $!;
    }
    print "gatk output in $dir\n" ;
    
    my $username = getlogin;

    #launching blade jobs 
    #Run GATK
    my $jobid1=`bsub -N -u $user\@genome.wustl.edu -q long -e $genome_name.gatk.err -J '$genome_name gatk' -e $dir\/$genome_name.err 'perl -I /gsc/scripts/opt/genome/current/pipeline/lib/perl/ \`which gmt\` gatk somatic-indel --normal-bam $normal_bam --tumor-bam $tumor_bam --output-file $dir\/$genome_name.GATK.indel --formatted-file $dir\/$genome_name.GATK.formatted --somatic-file $dir\/$genome_name.GATK.somatic' `;
    print "$jobid1\n";
    ($jobid1)=($jobid1=~m/<(\d+)>/);
    print "$jobid1\n";

    #Run GATK Annotation
    my $jobid2=`bsub -N -u $user\@genome.wustl.edu -q long -J '$genome_name.gatk.anno' -w 'ended($jobid1)' -e $dir\/$genome_name.anno.err 'perl -I /gsc/scripts/opt/genome/current/pipeline/lib/perl/ \`which gmt\` annotate transcript-variants --variant-file $dir\/$genome_name.GATK.somatic --output-file $dir\/$genome_name.GATK.somatic.anno --annotation-filter top'`;
    print "$jobid2\n";
    ($jobid2)=($jobid2=~m/<(\d+)>/);
    print "$jobid2\n";

    my $jobid3=`bsub -N -u $user\@genome.wustl.edu -q long -J '$genome_name.gatk.ucsc.anno' -w 'ended($jobid1)' -e $dir\/$genome_name.ucsc.anno.err 'perl -I /gsc/scripts/opt/genome/current/pipeline/lib/perl/ \`which gmt\` somatic ucsc-annotator --input-file $dir\/$genome_name.GATK.somatic --output-file $dir\/$genome_name.GATK.somatic.ucsc.anno --unannotated-file $dir\/$genome_name.GATK.somatic.ucsc.unanno'`;
    print "$jobid3\n";
    ($jobid3)=($jobid3=~m/<(\d+)>/);
    print "$jobid3\n";

    #Run Tiering
    my $jobid4=`bsub -N -u $user\@genome.wustl.edu -q short -J '$dir\/$genome_name.gatk.anno' -w 'ended($jobid2) && ended($jobid3)' -e $dir\/$genome_name.tier.err 'perl -I /gsc/scripts/opt/genome/current/pipeline/lib/perl/ \`which gmt\` somatic tier-variants --variant-file $dir\/$genome_name.GATK.somatic --transcript-annotation-file $dir\/$genome_name.GATK.somatic.anno --ucsc-file $dir\/$genome_name.GATK.somatic.ucsc.anno --tier1-file $dir\/$genome_name.GATK.somatic.anno.tier1 --tier2-file $dir\/$genome_name.GATK.somatic.anno.tier2 --tier3-file $dir\/$genome_name.GATK.somatic.anno.tier3  --tier4-file $dir\/$genome_name.GATK.somatic.anno.tier4' `;
    print "$jobid4\n";
    ($jobid4)=($jobid4=~m/<(\d+)>/);
    print "$jobid4\n";
    return 1;
}


1;

sub help_brief {
    "Helps run gakt by create directory according to cases and submit the jobs"
}

sub help_detail {
    <<'HELP';
This script helps runs gatk indel predictions. It uses a somatic model to locate the tumor and normal bam files and then launches the gatk job.
Example: genome somatic command run-gatk-build --analysis-dir /gscmnt/sata197/info/medseq/PCGP_Analysis/Indels/ --build-id 104295099
HELP
}

