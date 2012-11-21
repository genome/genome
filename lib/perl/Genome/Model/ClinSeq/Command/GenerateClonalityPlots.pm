package Genome::Model::ClinSeq::Command::GenerateClonalityPlots;
#Written by Malachi Griffith and Nate Dees
#Updated into tool form by Scott Smith

use strict;
use warnings;
use Genome;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use Genome::Model::ClinSeq::Util qw(:all);


class Genome::Model::ClinSeq::Command::GenerateClonalityPlots {
    is => 'Command::V2',
    has_input => [
        somatic_var_build   => { is => 'Genome::Model::Build::SomaticVariation', id_by => 'somatic_var_build_id', 
                                doc => 'Build ID for a somatic variation model' },

        common_name         => { is => 'Text', 
                                doc => 'Human readable name for the patient / sample comparison' },
        
        output_dir          => { is => 'Text', 
                                doc => 'Directory to place temp files and results' },

        read_counts         => { is => 'Text', is_optional => 1,
                                doc => 'Instead of generating read counts, use this pre-prepared readcount file' },
    ],
    has_param => [
        verbose             => { is => 'Boolean', is_optional => 1, default_value => 0,
                                doc => 'To display more output, set this flag' },
        
        limit               => { is => 'Number', is_optional => 1,
                                doc => 'Limit the number of SNVs to the first N (mostly for testing).' }

    ],
    doc => "This script attempts to automate the process of creating a 'clonality' plot"
};

sub help_synopsis {
    return <<INFO
  genome model clin-seq generate-clonality-plots --somatic-var-build=129396826  --output-dir=/tmp/output --common-name='AML54'  --verbose
INFO
}

sub status_message {
    my $self = shift;
    my $msg = shift;
    $msg = BLUE . $msg . RESET;
    return $self->SUPER::status_message($msg);
}

sub warning_message {
    my $self = shift;
    my $msg = shift;
    $msg = YELLOW . $msg . RESET;
    return $self->SUPER::status_message($msg);
}

sub error_message {
    my $self = shift;
    my $msg = shift;
    $msg = RED . $msg . RED;
    return $self->SUPER::status_message($msg);
}

sub execute {
    my $self = shift;

    #This script running a series of commands obtained from Nate Dees that results in the creation of a clonality plot (.pdf)
    my $somatic_var_build = $self->somatic_var_build; 
    my $output_dir = $self->output_dir;
    my $common_name = $self->common_name;
    my $verbose = $self->verbose;
    my $limit = $self->limit;

    if (not defined $limit) {
        $self->status_message("limit is not defined");
    }
    else {
        $self->status_message("limit is $limit");
    }

    if ($verbose){print BLUE, "\n\nCreating clonality plot for $common_name", RESET;}

    # This tool calls some scripts which have not been converted into tools
    my $script_dir = Cwd::abs_path(File::Basename::dirname(__FILE__) . '/../original-scripts/') . '/';
    unless (-d $script_dir) {
        die $self->error_message("failed to find script dir $script_dir!")
    }

    #Get somatic variation effects dir, tumor bam and normal bam from a somatic variation model ID
    my %data_paths;
    #... /genome/lib/perl/Genome/Model/Build/SomaticVariation.pm
    $data_paths{root_dir} = $somatic_var_build->data_directory ."/";
    $data_paths{effects_dir} = "$data_paths{root_dir}"."effects/";
    $data_paths{cnvs_hq} = "$data_paths{root_dir}"."variants/cnvs.hq";
    $data_paths{normal_bam} = $somatic_var_build->normal_bam;
    $data_paths{tumor_bam} = $somatic_var_build->tumor_bam;
    my $reference_build = $somatic_var_build->reference_sequence_build;
    $data_paths{reference_fasta} = $reference_build->full_consensus_path('fa');
    $data_paths{display_name} = $reference_build->__display_name__;

    my $somatic_effects_dir = $data_paths{effects_dir};

    #Make sure the specified parameters are correct
    $somatic_effects_dir = &checkDir('-dir'=>$somatic_effects_dir, '-clear'=>"no");
    $output_dir = &checkDir('-dir'=>$output_dir, '-clear'=>"no");

    #Step 1 - gather the tier 1-3 snv files from the build:
    my $tier1_snv_file = $somatic_effects_dir . "snvs.hq.novel.tier1.v2.bed";
    my $tier2_snv_file = $somatic_effects_dir . "snvs.hq.novel.tier2.v2.bed";
    my $tier3_snv_file = $somatic_effects_dir . "snvs.hq.novel.tier3.v2.bed";
    my $cp_cmd = "cp $tier1_snv_file $tier2_snv_file $tier3_snv_file $output_dir";
    if ($verbose){print YELLOW, "\n\n$cp_cmd", RESET;}
    Genome::Sys->shellcmd(cmd => $cp_cmd);

    #Step 2 - put them together in one file:
    my $snv_file = $output_dir . "allsnvs.hq.novel.tier123.v2.bed";
    my $cat_cmd;
    if (defined $limit) {
        $self->warning_message("limiting SNVs to the first $limit from the combined list!");
        $cat_cmd = "cat $output_dir"."snvs* | head -n $limit > $snv_file";
    }
    else {
        $cat_cmd = "cat $output_dir"."snvs* > $snv_file";
    }
    if ($verbose){print YELLOW, "\n\n$cat_cmd", RESET;}
    Genome::Sys->shellcmd(cmd => $cat_cmd);

    #Step 3 - take it out of bed format to be fed into bam-readcounts:
    my $adapted_file ="$output_dir"."allsnvs.hq.novel.tier123.v2.bed.adapted";
    my $awk_cmd = "awk \'{OFS=\"\\t\";FS=\"\\t\";}{print \$1,\$3,\$3,\$4}\' $snv_file | sed \'s/\\//\\t/g\' > $adapted_file";
    if ($verbose){print YELLOW, "\n\n$awk_cmd", RESET;}
    Genome::Sys->shellcmd(cmd => $awk_cmd);

    #Define the BAM files.  
    #The 'old' method supplied both Tumor & Normal coverage and both would be used to assess a minimum coverage cutoff for plotting
    #The 'new' method uses only the Tumor coverage
    my $tumor_bam = $data_paths{tumor_bam};
    my $normal_bam = $data_paths{normal_bam};

    #Step 4 - run bam readcounts and assess the particular reads for the reference and variant and print out details about the numbers of reads and the percentages for multiple bam files:
    my $readcounts_outfile;
    if ($self->read_counts) {
        $readcounts_outfile = $self->read_counts;
    }
    else {
        $readcounts_outfile = "$adapted_file".".readcounts";
        my $read_counts_cmd = "$script_dir"."borrowed/ndees/give_me_readcounts.pl  --sites_file=$adapted_file --bam_list=\"Tumor:$tumor_bam,Normal:$normal_bam\" --reference_fasta=$data_paths{reference_fasta} --output_file=$readcounts_outfile";
        if ($verbose){print YELLOW, "\n\n$read_counts_cmd", RESET;}
        Genome::Sys->shellcmd(cmd => $read_counts_cmd);
    }

    #Step 5 - create a varscan-format file from these outputs:
    #perl ~kkanchi/bin/create_pseudo_varscan.pl     allsnvs.hq.novel.tier123.v2.bed.adapted     allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts     >     allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts.varscan
    my $readcounts_varscan_file = "$readcounts_outfile".".varscan";
    my $varscan_format_cmd = "$script_dir"."borrowed/kkanchi/create_pseudo_varscan.pl $adapted_file $readcounts_outfile > $readcounts_varscan_file";
    if ($verbose){print YELLOW, "\n\n$varscan_format_cmd", RESET;}
    Genome::Sys->shellcmd(cmd => $varscan_format_cmd);

    #TODO: Replace steps 4-5 above by using the following script:
    #TODO: Once you know this is working, use $readcounts_clonality_outfile instead of $readcounts_varscan_file in the clonality commands below.  Then comment out steps 2-5 above
    #gmt validation prepare-wgs-for-clonality-plot --help
    #USAGE
    # gmt validation prepare-wgs-for-clonality-plot --output-file=? --snv-file=? [--bam-file=?]
    #    [--genome-build=?] [--min-mapping-quality=?] [--output-readcounts-file=?] [--readcounts-file=?]
    #Use the optional --bam-file input so that readcounts are generated for you.
    my $readcounts_clonality_outfile = "$output_dir"."readcounts.clonality";
    my $readcounts_formatted_outfile = "$output_dir"."readcounts.formatted";
    my $prepare_cmd = "gmt validation prepare-wgs-for-clonality-plot --output-file=$readcounts_clonality_outfile --snv-file=$adapted_file --bam-file=$tumor_bam --genome-build=$data_paths{reference_fasta} --output-readcounts-file=$readcounts_formatted_outfile";
    if ($verbose){print YELLOW, "\n\n$prepare_cmd", RESET;}
    Genome::Sys->shellcmd(cmd => $prepare_cmd);

    #Step 6 - Take the cnvs.hq file from the somatic-variation build, and run the cna-seg tool to create known regions of copy-number
    #Specify config file paths for hg19/build37
    #gmt copy-number cna-seg --copy-number-file=/gscmnt/ams1184/info/model_data/2875816457/build111674790/variants/cnvs.hq  --min-markers=4  --detect-somatic  --centromere-file=/gscmnt/sata186/info/medseq/kchen/work/SolexaCNV/scripts/centromere.hg19.csv  --gap-file=/gscmnt/sata186/info/medseq/kchen/work/SolexaCNV/scripts/hg19gaps.csv  --output-file=hg1.cnvhmm

    #Make a copy of the cnvs.hq file
    $cp_cmd = "cp $data_paths{cnvs_hq} $output_dir";
    if ($verbose){print YELLOW, "\n\n$cp_cmd", RESET;}
    Genome::Sys->shellcmd(cmd => $cp_cmd);
    my $chmod_cmd = "chmod 664 $output_dir"."cnvs.hq";
    Genome::Sys->shellcmd(cmd => $chmod_cmd);

    my $centromere_file;
    my $gap_file;
    my $clinseq_annotations_dir = "/gscmnt/sata132/techd/mgriffit/reference_annotations/";
    if ($data_paths{display_name} =~ /NCBI\-human\-build36/){
        $centromere_file = $clinseq_annotations_dir . "hg18/ideogram/centromere.hg18.csv";
        $gap_file = $clinseq_annotations_dir . "hg18/ideogram/hg18gaps.csv";
    }elsif($data_paths{display_name} =~ /GRCh37\-lite\-build37/){
        $centromere_file = $clinseq_annotations_dir . "hg19/ideogram/centromere.hg19.csv";
        $gap_file = $clinseq_annotations_dir . "hg19/ideogram/hg19gaps.csv";
    }else{
        print RED, "\n\nUnrecognized build - unable to identify centromere and gapfiles, you will need to generate these and place in the appropriate location\n\n", RESET;
        return;
    }
    my $cnvhmm_file = "$output_dir"."cnaseq.cnvhmm";
    my $cnaseg_cmd = "gmt copy-number cna-seg --copy-number-file=$data_paths{cnvs_hq}  --min-markers=4  --detect-somatic  --centromere-file=$centromere_file  --gap-file=$gap_file  --output-file=$cnvhmm_file";
    if ($verbose){print YELLOW, "\n\n$cnaseg_cmd", RESET;}
    Genome::Sys->shellcmd(cmd => $cnaseg_cmd);


    my $varscan_file = $readcounts_varscan_file;
    #my $varscan_file = $readcounts_clonality_outfile;

    #Step 7 - then, put the cna-seg and varscan-format snv file together in this clonality tool:
    #gmt validation clonality-plot     --cnvhmm-file     /gscuser/ndees/103/wgs/SV_somatic/CNV/aml103.cnvhmm     --output-image     aml103.clonality.pdf     --r-script-output-file     clonality.R     --varscan-file     allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts.varscan     --analysis-type     wgs     --sample-id     'AML103'     --positions-highlight     IL2RA_NF1_positions

    #gmt validation clonality-plot  --cnvhmm-file='/gscmnt/sata132/techd/mgriffit/hg1/clonality/hg1.cnvhmm'  --output-image hg1.clonality.pdf  --r-script-output-file clonality.R  --varscan-file allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts.varscan  --analysis-type wgs  --sample-id 'HG1'
    #Without clusters
    my $output_image_file1 = "$output_dir"."$common_name".".clonality.pdf";
    my $r_script_file = "$output_dir"."clonality.R";
    my $uc_common_name = uc($common_name);
    my $clonality_cmd1 = "gmt validation clonality-plot  --cnvhmm-file=$cnvhmm_file  --output-image=$output_image_file1  --r-script-output-file=$r_script_file  --varscan-file=$varscan_file  --analysis-type=wgs  --sample-id='$uc_common_name'";
    if ($verbose){print YELLOW, "\n\n$clonality_cmd1\n", RESET;}
    Genome::Sys->shellcmd(cmd => $clonality_cmd1);

    #With clusters
    my $clustered_data_output_file = $output_dir . $common_name . ".clustered.data.tsv";
    my $output_image_file2 = "$output_dir"."$common_name".".clonality.clusters.pdf";
    my $clonality_cmd2 = "gmt validation clonality-plot  --cnvhmm-file=$cnvhmm_file  --output-image=$output_image_file2  --r-script-output-file=$r_script_file  --varscan-file=$varscan_file  --analysis-type=wgs  --sample-id='$uc_common_name'  --plot-clusters  --clustered-data-output-file=$clustered_data_output_file";
    if ($verbose){print YELLOW, "\n\n$clonality_cmd2\n", RESET;}

    #TODO: until the --plot-clusters functionality is more stable, allow a failed exit code
    Genome::Sys->shellcmd(cmd => $clonality_cmd2, allow_failed_exit_code => 1);

    #Keep the files that were needed to run the cna-seg and clonality plot steps so that someone can rerun with different parameters 
    #Delete intermediate files though?

    if ($verbose){print "\n\n";}

    return 1;
};

1;

