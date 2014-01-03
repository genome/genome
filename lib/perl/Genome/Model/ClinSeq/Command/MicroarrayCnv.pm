package Genome::Model::ClinSeq::Command::MicroarrayCnv;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::MicroarrayCnv {
    is => 'Command::V2',
    has_input => [
        microarray_model_tumor => {
            is => 'Genome::Model::GenotypeMicroarray',
            doc => 'Tumor microarray model',
            is_optional => 1,
        },
        microarray_model_normal => {
            is => 'Genome::Model::GenotypeMicroarray',
            doc => 'Normal microarray model',
            is_optional => 1,
        },
        clinseq_model => {
            is => 'Genome::Model::ClinSeq',
            doc => 'Clinseq model(the microarray models will be extracted from the somatic-exome or somatic-wgs models of this model)',
            is_optional => 1,
        },
        outdir => {
            is => 'FilesystemPath',
            doc => 'Directory where output files will be written',
        },
        test => {
            is => 'Boolean',
            doc => 'True for tests, just makes CNView plots for kinase gene symbol list',
            default_value => 0,
        	is_optional => 1,
        },
    ],
    has_optional_input => [
        copynumber_tumor => {
            is => 'FilesystemPath',
            doc => '.copynumber file for tumor'
        },
        copynumber_normal => {
            is => 'FilesystemPath',
            doc => '.copynumber file for tumor'
        },
        cancer_annotation_db => {
            is => 'Genome::Db::Tgi::CancerAnnotation',
            doc => 'cancer annotation db to be used.(not required if using a clinseq model as input)',
            example_values  => ['tgi/cancer-annotation/human/build37-20130401.1'],
        },
        annotation_build_id => {
            is => 'Genome::Model::Build::ImportedAnnotation', 
            doc => 'Supply an annotation build id (e.g., 124434505 for NCBI-human.ensembl/67_37l_v2)',
            default_value => '124434505',
        },
    ],
    doc => 'Create somatic CopyNumber plots using MicroArray build files with CnView',
};

sub help_synopsis {
    return <<EOS
        genome model clin-seq microarray-cnv --outdir=/gscuser/gscuser1/tmp/ --clinseq-model=2887519760
        genome model clin-seq microarray-cnv --microarray-model-tumor=2878860299 --microarray-model-normal=2878860274 --outdir=/gscuser/gscuser1/tmp/ --cancer-annotation-db='tgi/cancer-annotation/human/build37-20130401.1'
        genome model clin-seq microarray-cnv --copynumber-tumor=/gscuser/aramu/tumor.HCC1.original --copynumber-normal=/gscuser/aramu/normal.HCC1.original--outdir=/gscuser/gscuser1/tmp/ --cancer-annotation-db='tgi/cancer-annotation/human/build37-20130401.1'
EOS
}

sub help_detail {
    return <<EOS
Create somatic copynumber plots using the ".original" files in the microarray builds. This tool can either take in two microarray models as input or it can take in a Clin-Seq model (with an exome-sv or wgs-sv model) as input or it can take in two ".original" files as input. In the second case the tool extracts the ".original" file by looking at the microarray models used in the exome-sv or wgs-sv model of the clinseq model. The segmentation is done using 'gmt copy-number cbs' and the copy number plots are made using 'gmt copy-number cn-view'. 
EOS
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);

    unless (-e $self->outdir && -d $self->outdir) {
        push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['outdir'],
	                                          desc => "Outdir: " . $self->outdir . " not found or not a directory",
                                            );
    }
    return @errors;
}

sub get_annotation_db {
    my $self = shift;
    my $cancer_annotation_db;
    if($self->clinseq_model) {
        $cancer_annotation_db = $self->clinseq_model->cancer_annotation_db;
    }
    elsif($self->cancer_annotation_db) {
        $cancer_annotation_db = $self->cancer_annotation_db;
        $self->status_message("Using cancer annotation db: " . $self->cancer_annotation_db);
    }
    else {
        die $self->error_message("Please specify cancer annotation db");
    }
    return $cancer_annotation_db;
}

sub get_microarray_models {
    my $self = shift;
    if($self->microarray_model_tumor and $self->microarray_model_tumor) {
        return ($self->microarray_model_tumor, $self->microarray_model_normal);
    } elsif($self->clinseq_model) {
        my $base_model;
        if($self->clinseq_model->exome_model) {
            $base_model = $self->clinseq_model->exome_model;
        } elsif($self->clinseq_model->wgs_model) {
            $base_model = $self->clinseq_model->wgs_model;
        } else {
            die $self->error_message("Please specify a clinseq model with either an exome-sv [or] WGS-sv model");
        }
        my $microarray_model_tumor = $base_model->tumor_model->genotype_microarray;
        my $microarray_model_normal = $base_model->normal_model->genotype_microarray;
        return $microarray_model_tumor, $microarray_model_normal;
    } else {
        die $self->error_message("Please specify one clinseq model [or] two microarray models as input!");
    }
}

sub get_copynumber_files {
    my $self = shift;
    my $copynumber_tumor;
    my $copynumber_normal;
    if($self->copynumber_tumor and $self->copynumber_normal) {
        if(-e $self->copynumber_tumor and -e $self->copynumber_normal) {
            $copynumber_tumor = $self->copynumber_tumor;
            $copynumber_normal = $self->copynumber_normal;
        } else {
            die $self->error_message("Unable to find one of the copynumber files");
        }
    } else {
        my ($microarray_tumor, $microarray_normal) = $self->get_microarray_models();
        if(-e $microarray_tumor->last_succeeded_build->original_genotype_file_path) {
            $copynumber_tumor = $microarray_tumor->last_succeeded_build->original_genotype_file_path;
        } else {
            die $self->error_message("Unable to find copynumber file for " . $microarray_tumor->name);
        }
        if(-e $microarray_normal->last_succeeded_build->original_genotype_file_path) {
            $copynumber_normal = $microarray_normal->last_succeeded_build->original_genotype_file_path;
        } else {
            die $self->error_message("Unable to find copynumber file for " . $microarray_normal->name);
        }
    }
    Genome::Sys->copy_file($copynumber_tumor, $self->outdir."/tumor.copynumber.original");
    Genome::Sys->copy_file($copynumber_normal, $self->outdir."/normal.copynumber.original");
    return ($copynumber_tumor, $copynumber_normal);
}

sub run_cbs {
    my $self = shift;
    my ($tumor_copynumber, $normal_copynumber, $cbs_op) = @_;
    my $cnv_diff_file = $self->outdir . "/cnvs.diff";
    my $create_diff_cmd;
    #copynumber = 2^(log_r_ratio + 1)
    if(not $self->test) {
        $create_diff_cmd = 'paste ' . $tumor_copynumber . ' ' . $normal_copynumber . ' | awk \'!/NaN|chr/ { print $1"\t"$2"\t"2^($6 + 1) - 2^($17 + 1) }\' > ' . $cnv_diff_file;
    } else { #use only chr1 for test
        $create_diff_cmd = 'paste ' . $tumor_copynumber . ' ' . $normal_copynumber . ' | awk \'!/NaN|chr/ { if($1 == 6) print $1"\t"$2"\t"2^($6 + 1) - 2^($17 + 1) }\' > ' . $cnv_diff_file;
    }
    Genome::Sys->shellcmd(cmd => $create_diff_cmd);
    my $cbs = Genome::Model::Tools::CopyNumber::Cbs->create(array_file => $cnv_diff_file, output_file => $cbs_op);
    $cbs->execute();
}

sub run_cnview {
    my $self = shift;
    my $tumor_copynumber = shift;
    my $normal_copynumber = shift;
    my $cbs_op = shift;
    my $cancer_annotation_db = shift;
    
    #Create cnvhq file
    my $cnv_file = $self->outdir . "/cnvs.hq";
    #copynumber = 2^(log_r_ratio + 1)
    my $create_cnvhq_cmd; 
    if(not $self -> test) { 
        $create_cnvhq_cmd = 'paste ' . $tumor_copynumber . ' ' .  $normal_copynumber . ' | awk \' BEGIN { print "CHR\tPOS\tTUMOR\tNORMAL\tDIFF"; } !/NaN|chr/ { print $1"\t"$2"\t"2^($6 + 1)"\t"2^($17 + 1)"\t"2^($6 + 1)-2^($17 + 1); }\' > ' . $cnv_file;
    } else {#use just chr1 for test
        $create_cnvhq_cmd = 'paste ' . $tumor_copynumber . ' ' .  $normal_copynumber . ' | awk \' BEGIN { print "CHR\tPOS\tTUMOR\tNORMAL\tDIFF"; } !/NaN|chr/ { if($1 == 6) print $1"\t"$2"\t"2^($6 + 1)"\t"2^($17 + 1)"\t"2^($6 + 1)-2^($17 + 1); }\' > ' . $cnv_file;
    }
    
    Genome::Sys->shellcmd(cmd => $create_cnvhq_cmd);
    
    #Create cnvhmm file
    my $cnv_hmm_file = $cbs_op . ".cnvhmm"; 
    #print only cnv segments with atleast five snp markers.
    my $make_hmmfile_cmd = 'awk \'{ size = $3-$2; nmarkers=size; event = "NA"; if($5>0) { event = "Gain" } else if($5<0) { event = "Loss" } cn1 = $5 +2; cn1 = int(cn1 + 0.5); cn2 = 2; if((event == "Gain" || event == "Loss") && cn1 !=2 && $4 >=5) print $1"\t"$2"\t"$3"\t"size"\t"nmarkers"\t"cn1"\t"cn1"\t"cn2"\t"cn2"\tNA\t"event; } \' ' . $cbs_op . ' > ' .  $cnv_hmm_file;
    Genome::Sys->shellcmd(cmd => $make_hmmfile_cmd);
    
    #For each list of gene symbols, run the CNView analysis
    my @cnv_symbols;
    if($self->test) {
        @cnv_symbols = qw(Kinase_dGene);
    } else {
        @cnv_symbols = qw (Kinase_dGene CancerGeneCensusPlus_Sanger AntineoplasticTargets_DrugBank All);
    }
    my $gene_symbol_dir = $cancer_annotation_db->data_directory . "/GeneSymbolLists/";
    foreach my $symbol(@cnv_symbols){
        my $symbol_outdir =  $self->outdir;
        if ($symbol eq "All") {
            my $cnview_cmd = Genome::Model::Tools::CopyNumber::CnView->create(annotation_build => $self->annotation_build_id, cnv_file => $cnv_file, segments_file => $cnv_hmm_file, output_dir => $symbol_outdir, name => $symbol, cancer_annotation_db => $cancer_annotation_db, window_size => 0, verbose => 1);
            $cnview_cmd->execute();
        } else {
            my $gene_targets_file = "$gene_symbol_dir/$symbol" . ".txt";
            my $cnview_cmd = Genome::Model::Tools::CopyNumber::CnView->create(annotation_build =>  $self->annotation_build_id, cnv_file => $cnv_file, segments_file => $cnv_hmm_file, output_dir => $symbol_outdir, gene_targets_file => $gene_targets_file, name => $symbol, cancer_annotation_db => $cancer_annotation_db, window_size => 0, verbose => 1);
            $cnview_cmd->execute();
        }
    }
}

sub execute {
    my $self = shift;
    my $copynumber_tumor;
    my $copynumber_normal;
    my $cbs_op = $self->outdir . "/cnvs.diff.cbs";
    my $cancer_annotation_db = $self->get_annotation_db();
    ($copynumber_tumor, $copynumber_normal) = $self->get_copynumber_files();
    $self->run_cbs($copynumber_tumor, $copynumber_normal, $cbs_op);
    $self->run_cnview($copynumber_tumor, $copynumber_normal, $cbs_op, $cancer_annotation_db);
    return 1;
}

1;
