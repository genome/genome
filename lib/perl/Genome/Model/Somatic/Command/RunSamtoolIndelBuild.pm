package Genome::Model::Somatic::Command::RunSamtoolIndelBuild;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::Model::Somatic::Command::RunSamtoolIndelBuild {
	is => 'Command',
	has => [
    	somatic_build_id => { type => 'String', is_optional => 0, doc => "somatic build id to process", },
        analysis_dir => { type => 'String', is_optional => 0, doc => "Directory where the recurrent gene output will be", },
        ]
};

sub help_brief {
    	"Generates somtools predicted somatic indels (without assembly)"
}

sub help_detail {
    	<<'HELP';
Hopefully this script will run the run tiered samtools predicted somatic indels.
HELP
}

sub execute {
	my $self=shift;
    	my $somatic_build_id=$self->somatic_build_id;
    	my $analysis_dir=$self->analysis_dir;
    	my $build = Genome::Model::Build->get($somatic_build_id);
    	
    	my %indel_tiers;
    	my %lines;
	my %tumor_bam;
	my %normal_bam;

	my ($line,$sample, $chr, $pos, $gene, $change, $key)=("","","","","","","");
	my @column; my @lines; my @name;
 	
        # find bam files of somatic build and its common name & cancer type	
       	my $tumor_wgs_bam = $build->tumor_build->whole_rmdup_bam_file;
        my $normal_wgs_bam = $build->normal_build->whole_rmdup_bam_file;
	#next unless($tumor_build->model->subject->sub_type !~ /M[13]/);

        #satisfied we should start doing stuff here
        my $data_directory = $build->data_directory . "/";

       	unless(-d $data_directory) {
       		$self->error_message("$data_directory is not a directory");
       		return;
       	}

 #       my $indel_transcript_annotation = "$data_directory/annotate_output_indel.out";

	%indel_tiers=(
		"Tier1" => "$data_directory/t1i_tier1_indel.csv",
		"Tier2" => "$data_directory/t2i_tier2_indel.csv",
		"Tier3" => "$data_directory/t3i_tier3_indel.csv",
#	"Tier4" => "$data_directory/t4i_tier4_indel.csv",		
	);
		
       	my $common_name = $build->tumor_build->model->subject->source_common_name;
	print "$common_name\n";
	my $user = getlogin || getpwuid($<); #get current user name
        	# submit samtool assembly for each tier of indels for a sample
       	
       	my $sample_dir=$analysis_dir."/".$common_name."/samtools/";
	my $normal_dir=$analysis_dir."/".$common_name."/samtools/normal";
	my $tumor_dir=$analysis_dir."/".$common_name."/samtools/tumor";
	`mkdir $sample_dir`;
	`mkdir $normal_dir`;
	`mkdir $tumor_dir`; 
       	foreach my $tier (keys %indel_tiers) {
       		my $indel_file = $indel_tiers{$tier};
       		unless (-e $indel_file){
	       		$self->error_message("The $indel_file doesn't exist, exit now!");
       			return;
       		}
       		
       		# Normal data
       		my $jobid1 ="bsub -N -u $user\@genome.wustl.edu -J $common_name.$tier.N -R 'select[type==LINUX64]' 'gmt somatic assemble-indel --assembly-indel-list=$sample_dir/$tier.normal --bam-file=$normal_wgs_bam --data-directory=$normal_dir --indel-file=$indel_file'";
		print "$jobid1\n"; $jobid1=`$jobid1`;
       		$jobid1=($jobid1=~ /<(\d+)>/);
       		
       		my $jobid2 ="bsub -N -u $user\@genome.wustl.edu -J $common_name.$tier.NT -w 'ended($jobid1)' 'grep -v \"NT\" $sample_dir/$tier.normal > $sample_dir/$tier.normal.noNT'";
		print "$jobid2\n"; $jobid2=`$jobid2`;
       		$jobid2=($jobid2=~ /<(\d+)>/);
		
		#Tumor data
           	my $jobid3 ="bsub -N -u $user\@genome.wustl.edu -J $common_name.$tier.T -R 'select[type==LINUX64]' 'gmt somatic assemble-indel --assembly-indel-list=$sample_dir/$tier.tumor --bam-file=$tumor_wgs_bam --data-directory=$tumor_dir --indel-file=$indel_file'";
		print "$jobid3\n"; $jobid3=`$jobid3`;
       		$jobid3=($jobid3=~ /<(\d+)>/);
		       		
       		my $jobid4 ="bsub -N -u $user\@genome.wustl.edu -J $common_name.$tier.NT -w 'ended($jobid3)' 'grep -v \"NT\" $sample_dir/$tier.tumor > $sample_dir/$tier.tumor.noNT'";
		print "$jobid4\n"; $jobid4=`$jobid4`;
       		$jobid4=($jobid4=~ /<(\d+)>/);
		
            	# intersect
            	my $jobid5="bsub -N -u $user\@genome.wustl.edu -J $common_name.$tier.N.somatic -R 'select\[type==LINUX64\]' -w 'ended($jobid2) && ended($jobid4)' 'gmt somatic intersect-assembled-indels --normal-indel-file=$sample_dir/$tier.normal.noNT --tumor-indel-file=$sample_dir/$tier.tumor.noNT --tumor-assembly-data-directory=$tumor_dir --germline-output-list=$sample_dir/$tier.noNT.germline  --somatic-output-list=$sample_dir/$tier.noNT.somatic\'";
		print "$jobid5\n";
            	$jobid5=($jobid5=~ /<(\d+)>/);
       		print "intersect tumor/normal $tier in job $jobid5\n";
       	}

        return 1;
}


1;


