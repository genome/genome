# % Last Change: Wed Apr 16 03:00 PM 2014 C
package Genome::Model::Tools::Analysis::Concordance;

use strict;
use Genome;
use IO::File;
#use Statistics::R;
use File::Basename;
use warnings;
require Genome::Sys;
use FileHandle;

class Genome::Model::Tools::Analysis::Concordance {
    is => 'Command',
    has => [
    bam_file_1 => {
        is => 'String',
        is_optional => 1,
        doc => 'Get indexed bam file #1',
    },

    bam_file_2 => {
        is => 'String',
        is_optional => 1,
        doc => 'Get indexed bam file #2',
    },

    model_id_1 => {
        is => 'String',
        is_optional => 1,
        doc => 'Get genome model ID #1',
    },

    model_id_2 => {
        is => 'String',
        is_optional => 1,
        doc => 'Get genome model ID #2',
    },

    reference_genome => {
        is => 'String',
        is_optional => 1,
        doc => 'Reference genome of indexed bam files',
    },

    snp_file => {
        is => 'String',
        is_optional => 0,
        doc => '1-based Tab-delimited file of SNP positions. Three columns requred: chr, start, end',
    },

    output_file => {
        is => 'String',
        is_optional => 1,
        doc => 'Output the numbers and percentages of matched and mismatched SNPs',
    }
    ]
};

sub help_brief {
    "check sample concordant using bam-readcount"
}

sub help_detail {
    "check sample concordant using bam-readcount"
}

####################################

sub execute {
    my $self = shift;
    my $model_id_1 = $self->model_id_1;
    my $model_id_2 = $self->model_id_2;
    my $bam_file_1 = $self->bam_file_1;
    my $bam_file_2 = $self->bam_file_2;
    my $reference_genome = $self->reference_genome;
    my $snp_file = $self->snp_file;
    my $output_file = $self->output_file;

    if( !defined($bam_file_1) && !defined($bam_file_2) && !defined($model_id_1) && !defined($model_id_2) ){
        die("Must provide indexed bam files OR model IDs.\n");
    }

    if( (defined($bam_file_1) || defined($bam_file_2)) && (defined($model_id_1) || defined($model_id_2)) ){
        die("Must provide indexed bam files OR model IDs.\n");
    }

    if( (!defined($bam_file_1) || !defined($bam_file_2)) && (!defined($model_id_1) && !defined($model_id_2)) ){
        die("Must provide 2 indexed bam files.\n");
    }

    if( (defined($bam_file_1) && defined($bam_file_2)) && !defined($reference_genome) && (!defined($model_id_1) && !defined($model_id_2)) ){
        die("Must provide reference genome of 2 bam files.\n");
    }

    if( (!defined($bam_file_1) && !defined($bam_file_2)) && (!defined($model_id_1) || !defined($model_id_2)) ){
        die("Must provide 2 model IDs.\n");
    }

    if( defined($model_id_1) && defined($model_id_2) ){
        # get model IDs
        my $model_1 = Genome::Model->get($model_id_1);
        my $model_2 = Genome::Model->get($model_id_2);

        # get the bam paths of the last succeeded build
        my $build_1 = $model_1->last_succeeded_build;
        my $bam_file_1 = $build_1->merged_alignment_result->bam_file;
        my $build_2 = $model_2->last_succeeded_build;
        my $bam_file_2 = $build_2->merged_alignment_result->bam_file;

        # get reference genome of the model
        my $ref_seq_build_id = $model_1->reference_sequence_build->build_id;
        my $ref_seq_build = Genome::Model::Build->get($ref_seq_build_id);
        my $reference_genome = $ref_seq_build->full_consensus_path('fa');
    }

    # dbSNP dbsnp_138.hg19.sort, Three columns requred: chr start end

    # create temp directory for munging
    my $inFile;
    my $outFile;
    my $tempdir = Genome::Sys->create_temp_directory();
    $tempdir or die "Unable to create temporary directory $!";
#    unless($tempdir) {
#        my $self->error_message("Unable to create temporary file $!");
#        die;
#    }

    # run bam-readcount
    my $bam_readcount_1 = Genome::Model::Tools::Sam::Readcount->execute(
        bam_file => $bam_file_1,
        reference_fasta => $reference_genome,
        region_list => $snp_file,
#        minimum_mapping_quality => $min_mapping_quality,
        output_file => "$tempdir/readcount_temp_file_1",
    );

    my $bam_readcount_2 = Genome::Model::Tools::Sam::Readcount->execute(
        bam_file => $bam_file_2,
        reference_fasta => $reference_genome,
        region_list => $snp_file,
#        minimum_mapping_quality => $min_mapping_quality,
        output_file => "$tempdir/readcount_temp_file_2",
    );

    # parse in perl # dump to temp directory
    # for i in *rct; do echo $i; cat $i | tr ":" "\t" | cut -f1-3,19,20,33,34,47,48,61,62,75,76 > $i.txt; done
    my $parse_file_1 = Genome::Sys->open_file_for_writing("$tempdir/parse_file_1");
    my $open_readcount_1 = Genome::Sys->open_file_for_reading("$tempdir/readcount_temp_file_1");
    while (my $line = $open_readcount_1->getline){
        $line =~ tr/:/\t/;
        my @field = split(/\t/,$line);
        $parse_file_1->print(join("\t", ($field[0],$field[1],$field[2],$field[18],$field[19],$field[32],$field[33],$field[46],$field[47],$field[60],$field[61],$field[74],$field[75])), "\n");
    }
    $open_readcount_1->close;
    $parse_file_1->close;

    my $parse_file_2 = Genome::Sys->open_file_for_writing("$tempdir/parse_file_2");
    my $open_readcount_2 = Genome::Sys->open_file_for_reading("$tempdir/readcount_temp_file_2");
    while (my $line = $open_readcount_2->getline){
        $line =~ tr/:/\t/;
        my @field = split(/\t/,$line);
        $parse_file_2->print(join("\t", ($field[0],$field[1],$field[2],$field[18],$field[19],$field[32],$field[33],$field[46],$field[47],$field[60],$field[61],$field[74],$field[75])), "\n");
    }
    $open_readcount_2->close;
    $parse_file_2->close;

    # run R script and output # goes in temp dir
    my $r_script_file = '/gscuser/tli/gc3018/bin/snp.R';
    my $cmd_1 = "Rscript $r_script_file '$tempdir/parse_file_1' '$tempdir/r_output_file_1'";
    my $return_1 = Genome::Sys->shellcmd(
        cmd => "$cmd_1",
    );
    unless($return_1) {
        $self->error_message("Failed to execute: $cmd_1.\n");
        die $self->error_message;
    }

    my $cmd_2 = "Rscript $r_script_file '$tempdir/parse_file_2' '$tempdir/r_output_file_2'";
    my $return_2 = Genome::Sys->shellcmd(
        cmd => "$cmd_2",
    );
    unless($return_2) {
        $self->error_message("Failed to execute: $cmd_2.\n");
        die $self->error_message;
    }

    # output using bash
    my $comp_sh = '/gscuser/tli/gc3018/bin/2way_comp.sh';
    my $normal = "$tempdir/r_output_file_1";
    my $pre = "$tempdir/r_output_file_2";
    #if (-e $output_file) {
    if (-e $normal && -e $pre) {
        my $cmd_3 = "bash $comp_sh $normal $pre $output_file";
        print "RUN: $cmd_3\n";
        system($cmd_3);
    }

#    unlink glob('readcount_temp_file_*');
#    unlink glob('parse_file_*');
#    unlink glob('r_output_file_*');
}

1;

