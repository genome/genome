package Genome::Model::Somatic::Command::RunFalsePositiveFilter;

use strict;
use warnings;

use Genome;
use IO::File;


class Genome::Model::Somatic::Command::RunFalsePositiveFilter {
	is => 'Command',
	has => [
	somatic_build_id => { type => 'String', is_optional => 0, doc => "the somatic build_id to process.", },
    	analysis_dir => { type => 'String', is_optional => 0, doc => "Directory where the filtered SNVs output will be", },
        ]
};

sub help_brief {
    	"Generates false positive filtered SNVs and its annotation for the last succeed builds of a somatic model"
}

sub help_detail {
    	<<'HELP';
Hopefully this script run false positive filtered SNVs and its annotation on the last succeed somatic builds in a model group for Tier1-Tier3 SNVs
HELP
}

sub execute {
	my $self=shift;
    	my $analysis_dir=$self->analysis_dir;
    	my $somatic_build_id =$self->somatic_build_id;
    	unless ($somatic_build_id){
    		$self->error_message("must have either somatic-build-id");
    		return;
    	}
    	my $build = Genome::Model::Build->get($somatic_build_id);
        unless(defined($build)) {
        $self->error_message("Unable to find build $somatic_build_id");
        return;
	}
    	my $model = $build->model;
    	unless(defined($model)) {
        	$self->error_message("Somehow this build does not have a model");
	        return;
    	}
    	unless($model->type_name eq 'somatic') {
        	$self->error_message("This build must be a somatic pipeline build");
        	return;
    	}
    	
    	#retrieve the tumor bam
    	my $tumor_bam = $build->tumor_build->whole_rmdup_bam_file;
    	unless($tumor_bam) {
        	$self->error_message("Couldn't determine tumor bam file from somatic model");
        	return;
    	}

    	if(-z $tumor_bam) {
        	$self->error_message("$tumor_bam is of size 0 or does not exist");
        	return;
    	}

	my $data_directory = $build->data_directory;
        unless(-d $data_directory) {
        	$self->error_message("$data_directory is not a directory");
        	return;
        }
        		
	my $common_name = $build->tumor_build->model->subject->source_common_name;

	my $user = getlogin || getpwuid($<); #get current user name
	my $out_dir=$analysis_dir."/".$common_name."_FPfiltered/";
# check whether the directory can be created or already exists		
	unless (-d "$out_dir"){
        	system ("mkdir -p $out_dir");
        }
        
        my $snp_transcript_annotation = "$data_directory/upload_variants_snp_1_output.out";
        my  %snv_tiers = (
        	"Tier1" => "$data_directory/tier_1_snp_file.out",
        	"hc2" => "$data_directory/tier_2_snp_high_confidence_file.out",
        	"hc3" => "$data_directory/tier_3_snp_high_confidence_file.out",
        );

	        # if not exist, check if using new files
	unless(-e $snp_transcript_annotation) {
        	$snp_transcript_annotation = "$data_directory/uv1_uploaded_tier1_snp.csv";
	       	%snv_tiers = (
	               	"Tier1" => "$data_directory/t1v_tier1_snp.csv",
	               	"hc2" => "$data_directory/hc2_tier2_snp_high_confidence.csv",
	               	"hc3" => "$data_directory/hc3_tier3_snp_high_confidence.csv",
        	);
        }

       	foreach my $tier (keys %snv_tiers) {
       		my $snv_file = $snv_tiers{$tier};
       		unless (-e $snv_file){
	       		$self->error_message("The $common_name: $snv_file doesn't exist, exit now!");
       			return;
       		}
       		
       		my $cmd ="";
        	$cmd ="bsub -N -u $user\@genome.wustl.edu -J $common_name.$tier.FP -R \'select\[type==LINUX64\]\' \'gmt somatic filter-false-positives --bam-file=$tumor_bam --variant-file=$snv_file --output-file=$out_dir/$tier.csv\'";
        	$cmd=`$cmd`;
		my ($jobid1) =($cmd=~ m/<(\d+)>/);
		print "$cmd\n$jobid1\n";
		my $cmd_anno="";
                $cmd_anno ="bsub -N -u $user\@genome.wustl.edu -J $common_name.$tier.anno -w \'ended\($jobid1\)\'  \'perl -I /gsc/scripts/opt/genome/current/pipeline/lib/perl/ \`which gmt\` annotate transcript-variants --use-version 2 --variant-file $out_dir/$tier.csv --output-file $out_dir/$tier.anno --annotation-filter top\'";
		$cmd_anno=`$cmd_anno`;
       		my ($jobid2) = ($cmd_anno=~ m/<(\d+)>/);
       		print "$cmd_anno\n$jobid2\n";
       		
       		
	}
        return 1;
}


1;


