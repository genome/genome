#!/usr/bin/env genome-perl

# P R A G M A S ###############################################################
use strict;
use warnings;

# M O D U L E S ###############################################################
use above 'Genome';
use Path::Class;
use Test::More;

# M A I N #####################################################################
#This is a bare minimum test that just compiles Perl and the UR class.
use_ok('Genome::Model::Tools::Transcriptome::ErccMapUnaligned');

SKIP: {
    skip "this is a long-running test & requires a human to view its outputs",
         2;
    #my $model = '39f7c31b4e28432f8cbee18c8a391d06'; # GRCh37 only align (Mix 2)
    #my $model = 'f8e2cf3aa1fe4e6da97370d8d7dd4f5a'; # GRCh37+ERCC align (Mix 2)
    my $model = '46fccb0ce8334d06bb09a26dd9903f9f'; # Heather's 1st Model [NO ERCC]
    #my $model = '576d3c68ae574779b92052daead2fa2d'; # Heather's 2nd Model [Mix1]
    my $ercc_dir = Path::Class::Dir->new(
        '/gscmnt/gc2764/cad/',
        'jwalker/scratch/ERCC/remap'
    );

    my $build = get_build($model);
    run_analysis($build, $ercc_dir);
}

done_testing();

# S U B R O U T I N E S #######################################################
sub run_analysis {
    my ($model, $ercc_dir) = @_;
    my $tool = Genome::Model::Tools::Transcriptome::ErccMapUnaligned->create(
        build           => $model,
        ercc_spike_in_file =>
          $ercc_dir->file('ERCC_Controls_Analysis.txt')->stringify,
    );
    ok($tool, 'Got a G::M::T::Transcriptome::ErccMapUnaligned instance');
    $tool->execute or die "[err] Trouble running tool!\n";
}

sub get_build {
    my $model_id = shift;
    my $model = Genome::Model->get(id => $model_id);
    my $build = $model->last_succeeded_build;
    unless ($build) {
        ok($build, 'got a build: ' . $build->id);
        done_testing();
        exit(0);
    }
    return $build;
}

sub get_bam {
    my $model_id = shift;
    my $model = Genome::Model->get(id => $model_id);
    my $build = $model->last_succeeded_build;
    my $bam = $build->merged_alignment_result->bam_path;
    ok($bam, 'got a bam path');
    ok(-e $bam, 'bam file exists');
    return $bam;
}
