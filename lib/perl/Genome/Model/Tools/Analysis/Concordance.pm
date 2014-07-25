# % Last Change: Fri Jul 25 02:44:09 PM 2014 CDT
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
        doc => '1-based Tab-delimited file of SNP positions. Three columns required: chr, start, end',
    },

    output_file => {
        is => 'String',
        is_optional => 1,
        doc => 'Output the numbers and percentages of matched and mismatched SNPs',
    }
    ]
};

sub help_synopsis {
  return <<EOS
    gmt analysis concordance -snp-file=/gscuser/tli/gc3018/db/dbsnp_138.hg19.sort.uniq.chr21 --bam-file-1=/gscuser/tli/168_norm_chr21.10k.bam --bam-file-2=/gscuser/tli/168_pre_chr21.10k.bam --output-file=/gscuser/tli/168_comp.chr21 --reference-genome=/gscuser/tli/21.fa

    gmt analysis concordance -snp-file=/gscuser/tli/gc3018/db/dbsnp_138.hg19.sort.uniq.chr21 --model-id-1=9b7fa0a9e3c643379c26a3367de8ace0 --model-id-2=a1f870cb759548339c841393ef823b18 --output-file=model.chr21
EOS
}

sub help_brief {
    "check sample concordant using bam-readcount"
}

sub help_detail {
    "check sample concordant using bam-readcount"
}

####################################

sub execute {
    my $self = shift;
    my $bam_file_1 = $self->bam_file_1;
    my $bam_file_2 = $self->bam_file_2;
    my $model_id_1 = $self->model_id_1;
    my $model_id_2 = $self->model_id_2;
    my $reference_genome = $self->reference_genome;
    my $snp_file = $self->snp_file;
    my $output_file = $self->output_file;

    if (!defined($bam_file_1) && !defined($bam_file_2) && !defined($model_id_1) && !defined($model_id_2)){
        die("Must provide indexed bam files OR model IDs.\n");
    }

    if ((defined($bam_file_1) || defined($bam_file_2)) && (defined($model_id_1) || defined($model_id_2))){
        die("Must provide indexed bam files OR model IDs.\n");
    }

    if ((!defined($bam_file_1) || !defined($bam_file_2)) && (!defined($model_id_1) && !defined($model_id_2))){
        die("Must provide 2 indexed bam files.\n");
    }

    if ((!defined($bam_file_1) && !defined($bam_file_2)) && (!defined($model_id_1) || !defined($model_id_2))){
        die("Must provide 2 model IDs.\n");
    }

    if ((defined($bam_file_1) && defined($bam_file_2)) && !defined($reference_genome) && (!defined($model_id_1) && !defined($model_id_2))){
        die("Must provide reference genome of 2 bam files.\n");
    }

    if (defined($model_id_1) && defined($model_id_2)){
        # get model IDs
        my $model_1 = Genome::Model->get($model_id_1);
        my $model_2 = Genome::Model->get($model_id_2);

        # get the bam paths of the last succeeded build
        my $build_1 = $model_1->last_succeeded_build;
        my $build_2 = $model_2->last_succeeded_build;
        $bam_file_1 = $build_1->merged_alignment_result->bam_file;
        $bam_file_2 = $build_2->merged_alignment_result->bam_file;

        # get reference genome of the model
        my $ref_seq_build_id = $model_1->reference_sequence_build->build_id;
        my $ref_seq_build = Genome::Model::Build->get($ref_seq_build_id);
        $reference_genome = $ref_seq_build->full_consensus_path('fa');
    }

    # index reference genome and bam files
    my $index_file = "$reference_genome.fai";
    if (!(-e $index_file)){
        die "Index file for reference ($index_file) not found!\n";
    }

    my $index_bam_file_1 = "$bam_file_1.bai";
    if (!(-e $index_bam_file_1)){
        die "Index file for bam ($index_bam_file_1) not found!\n";
    }

    my $index_bam_file_2 = "$bam_file_2.bai";
    if (!(-e $index_bam_file_2)){
        die "Index file for bam ($index_bam_file_2) not found!\n";
    }

    # dbSNP dbsnp_138.hg19.sort, Three columns required: chr start end

    # create temp directory for munging
    my $tempdir = Genome::Sys->create_temp_directory();
    $tempdir or die "Unable to create temporary directory $!";

    # run bam-readcount
    my $bam_readcount_1 = Genome::Model::Tools::Sam::Readcount->execute(
        bam_file => $bam_file_1,
        reference_fasta => $reference_genome,
        region_list => $snp_file,
#        minimum_mapping_quality => $min_mapping_quality,
        output_file => "$tempdir/norm_rct",
    );

    my $bam_readcount_2 = Genome::Model::Tools::Sam::Readcount->execute(
        bam_file => $bam_file_2,
        reference_fasta => $reference_genome,
        region_list => $snp_file,
#        minimum_mapping_quality => $min_mapping_quality,
        output_file => "$tempdir/pre_rct",
    );

    # empty result
    if (-z "$tempdir/norm_rct" || -z "$tempdir/pre_rct"){
        die("Have zero read counts.\n");
    }

    # parse in perl
    # for i in *rct; do echo $i; cat $i | tr ":" "\t" | cut -f1-3,19,20,33,34,47,48,61,62,75,76 > $i.txt; done
    my $parse_file_1 = Genome::Sys->open_file_for_writing("$tempdir/norm_parse");
    my $open_readcount_1 = Genome::Sys->open_file_for_reading("$tempdir/norm_rct");
    while (my $line = $open_readcount_1->getline){
        $line =~ tr/:/\t/;
        my @field = split(/\t/,$line);
        $parse_file_1->print(join("\t", ($field[0],$field[1],$field[2],$field[18],$field[19],$field[32],$field[33],$field[46],$field[47],$field[60],$field[61],$field[74],$field[75])), "\n");
    }
    $open_readcount_1->close;
    $parse_file_1->close;

    my $parse_file_2 = Genome::Sys->open_file_for_writing("$tempdir/pre_parse");
    my $open_readcount_2 = Genome::Sys->open_file_for_reading("$tempdir/pre_rct");
    while (my $line = $open_readcount_2->getline){
        $line =~ tr/:/\t/;
        my @field = split(/\t/,$line);
        $parse_file_2->print(join("\t", ($field[0],$field[1],$field[2],$field[18],$field[19],$field[32],$field[33],$field[46],$field[47],$field[60],$field[61],$field[74],$field[75])), "\n");
    }
    $open_readcount_2->close;
    $parse_file_2->close;

    # run R script and output # goes in temp dir
    my $dir_name = dirname(__FILE__);
    my $r_script_file = "\"" . $dir_name . "/Concordance.R\"";
    my $cmd_1 = "Rscript $r_script_file '$tempdir/norm_parse' '$tempdir/norm_r'";
    my $return_1 = Genome::Sys->shellcmd(
        cmd => "$cmd_1",
    );
    unless($return_1) {
        $self->error_message("Failed to execute: $cmd_1.\n");
        die $self->error_message;
    }

    my $cmd_2 = "Rscript $r_script_file '$tempdir/pre_parse' '$tempdir/pre_r'";
    my $return_2 = Genome::Sys->shellcmd(
        cmd => "$cmd_2",
    );
    unless($return_2) {
        $self->error_message("Failed to execute: $cmd_2.\n");
        die $self->error_message;
    }

    my $norm = "$tempdir/norm_r";
    my $pre = "$tempdir/pre_r";

    my (@norm_total, @norm_snp, @pre_total, @pre_snp);

    open (my $IN1,'<',$norm) or die "$!"; #open 1st file
    while (<$IN1>) {
        chomp;                            #remove newline from each row of the 1st file
        my @field = split(/\t/);          #read line from 1st file
        if ($field[13] ne "NA"){
            push(@norm_total,join('_',$field[0],$field[1]));
            push(@norm_snp,join('_',$field[0],$field[1],$field[13]));
        }
    }  
    close $IN1;

    open (my $IN2,'<',$pre) or die "$!";  #open 2nd file
    while (<$IN2>) {
        chomp;                            #remove newline from each row of the 2nd file
        my @field = split(/\t/);          #read line from 2nd file
        if ($field[13] ne "NA"){
            push(@pre_total,join('_',$field[0],$field[1]));
            push(@pre_snp,join('_',$field[0],$field[1],$field[13]));
        }
    }  
    close $IN2;

    # sort and uniq
    my %hashTemp;

    %hashTemp = map { $_ => 1 } @norm_total;
    my @norm_total_uniq = sort keys %hashTemp;

    %hashTemp = map { $_ => 1 } @norm_snp;
    my @norm_snp_uniq = sort keys %hashTemp;

    %hashTemp = map { $_ => 1 } @pre_total;
    my @pre_total_uniq = sort keys %hashTemp;

    %hashTemp = map { $_ => 1 } @pre_snp;
    my @pre_snp_uniq = sort keys %hashTemp;

    my @norm_pre_total = (@norm_total_uniq, @pre_total_uniq);
    my @norm_pre_snp = (@norm_snp_uniq, @pre_snp_uniq);

    %hashTemp = map { $_ => 1 } @norm_pre_total;
    my @norm_pre_total_uniq = sort keys %hashTemp;

    %hashTemp = map { $_ => 1 } @norm_pre_snp;
    my @norm_pre_snp_uniq = sort keys %hashTemp;

#    print "norm total = ", scalar(@norm_total), "\n";
#    print "norm total uniq = ", scalar(@norm_total_uniq), "\n";
#    print "norm snp = ", scalar(@norm_snp), "\n";
#    print "norm snp uniq = ", scalar(@norm_snp_uniq), "\n";

#    print "pre total = ", scalar(@pre_total), "\n";
#    print "pre total uniq = ", scalar(@pre_total_uniq), "\n";
#    print "pre snp = ", scalar(@pre_snp), "\n";
#    print "pre snp uniq = ", scalar(@pre_snp_uniq), "\n";

#    print "norm pre total = ", scalar(@norm_pre_total), "\n";
#    print "norm pre total uniq = ", scalar(@norm_pre_total_uniq), "\n";
#    print "norm pre snp = ", scalar(@norm_pre_snp), "\n";
#    print "norm pre snp uniq = ", scalar(@norm_pre_snp_uniq), "\n";

    my $a = scalar(@norm_snp_uniq)+scalar(@pre_snp_uniq)-scalar(@norm_pre_snp_uniq);
    my $b = scalar(@norm_total_uniq)+scalar(@pre_total_uniq)-scalar(@norm_pre_total_uniq);

#    print "Sample1 vs Sample2 [% matches] = ", "\t", $a, "\t", $b, "\t", $a/$b*100, "%", "\n";
#    print "Sample1 vs Sample2 [% mismatches] = ", "\t", $b-$a, "\t", $b, "\t", ($b-$a)/$b*100, "%", "\n";

    open (my $OUT,'>', $output_file) or die "$!";
    print $OUT "Sample1 = ", "\t", $bam_file_1, "\n";
    print $OUT "Sample2 = ", "\t", $bam_file_2, "\n";
#    print $OUT "Sample1 vs Sample2 [% matches] = ", "\t", $a, "\t", $b, "\t", $a/$b*100, "%", "\n";
    print $OUT "Sample1 vs Sample2 [% matches] = ", "\t";
    printf $OUT ("%10d\t", $a);
    printf $OUT ("%10d\t", $b);
    printf $OUT ("%7.3f\t", $a/$b*100);
    print $OUT "%\n";
#    print $OUT "Sample1 vs Sample2 [% mismatches] = ", "\t", $b-$a, "\t", $b, "\t", ($b-$a)/$b*100, "%", "\n";
    print $OUT "Sample1 vs Sample2 [% mismatches] = ", "\t";
    printf $OUT ("%10d\t", $b-$a);
    printf $OUT ("%10d\t", $b);
    printf $OUT ("%7.3f\t", ($b-$a)/$b*100);
    print $OUT "%\n";
    close $OUT;

    return 1;

#    unlink glob('readcount_temp_file_*');
#    unlink glob('parse_file_*');
#    unlink glob('r_output_file_*');

}

1;

