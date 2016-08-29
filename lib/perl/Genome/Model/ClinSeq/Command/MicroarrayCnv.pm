package Genome::Model::ClinSeq::Command::MicroarrayCnv;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::MicroarrayCnv {
    is        => 'Command::V2',
    has_param => [lsf_resource => {default => "-R 'span[hosts=1] rusage[mem=8000]' -M 8000000",}],
    has_input => [
        microarray_model_tumor => {
            is          => 'Genome::Model::GenotypeMicroarray',
            doc         => 'Tumor microarray model',
            is_optional => 1,
        },
        microarray_model_normal => {
            is          => 'Genome::Model::GenotypeMicroarray',
            doc         => 'Normal microarray model',
            is_optional => 1,
        },
        microarray_model_single => {
            is          => 'Genome::Model::GenotypeMicroarray',
            doc         => 'Microarray model of a single sample, not somatic.',
            is_optional => 1,
        },
        clinseq_build => {
            is => 'Genome::Model::Build::ClinSeq',
            doc =>
                'Clinseq build (the microarray models will be extracted from the somatic-exome or somatic-wgs builds of this build)',
            is_optional => 1,
        },
        outdir => {
            is  => 'FilesystemPath',
            doc => 'Directory where output files will be written',
        },
        min_cnv_diff => {
            is            => 'Float',
            doc           => 'Cutoff for the minimum cnv difference between tumor and normal [the absolute value]',
            default_value => 0.2,
            is_optional   => 1,
        },
        test => {
            is            => 'Boolean',
            doc           => 'True for tests, just makes CNView plots for kinase gene symbol list',
            default_value => 0,
            is_optional   => 1,
        },
    ],
    has_optional_input => [
        copynumber_tumor => {
            is  => 'FilesystemPath',
            doc => '.copynumber file for tumor'
        },
        copynumber_normal => {
            is  => 'FilesystemPath',
            doc => '.copynumber file for tumor'
        },
        cancer_annotation_db => {
            is             => 'Genome::Db::Tgi::CancerAnnotation',
            doc            => 'cancer annotation db to be used.(not required if using a clinseq model as input)',
            example_values => ['tgi/cancer-annotation/human/build37-20130401.1'],
        },
        annotation_build_id => {
            is            => 'Genome::Model::Build::ImportedAnnotation',
            doc           => 'Supply an annotation build id (e.g., 124434505 for NCBI-human.ensembl/67_37l_v2)',
            default_value => '124434505',
        },
    ],
    doc => 'create somatic copy number plots using genotype-microarray build files with CnView',
};

sub help_synopsis {
    return <<EOS
        genome model clin-seq microarray-cnv --outdir=/gscuser/gscuser1/tmp/ --clinseq-model=2887519760
        genome model clin-seq microarray-cnv --microarray-model-tumor=2878860299 --microarray-model-normal=2878860274 --outdir=/gscuser/gscuser1/tmp/ --cancer-annotation-db='tgi/cancer-annotation/human/build37-20130401.1'
        genome model clin-seq microarray-cnv --copynumber-tumor=/gscuser/aramu/tumor.HCC1.original --copynumber-normal=/gscuser/aramu/normal.HCC1.original--outdir=/gscuser/gscuser1/tmp/ --cancer-annotation-db='tgi/cancer-annotation/human/build37-20130401.1'
        genome model clin-seq microarray-cnv --microarray-model-single=cf3ab49f066945dd8956856c31aa8611 --outdir /tmp/ma_cnv/ --cancer-annotation-db='tgi/cancer-annotation/human/build37-20130401.1'
EOS
}

sub help_detail {
    return <<EOS
Create somatic copynumber plots using the ".original" files in the microarray builds. This tool can either take in two microarray models as input or it can take in a Clin-Seq model (with an exome-sv or wgs-sv model) as input or it can take in two ".original" files as input. In the second case the tool extracts the ".original" file by looking at the microarray models used in the exome-sv or wgs-sv model of the clinseq model. The segmentation is done using 'gmt copy-number cbs' and the copy number plots are made using 'gmt copy-number cn-view'.
EOS
}

sub __errors__ {
    my $self   = shift;
    my @errors = $self->SUPER::__errors__(@_);
    unless (-e $self->outdir && -d $self->outdir) {
        push @errors,
            UR::Object::Tag->create(
            type       => 'error',
            properties => ['outdir'],
            desc       => "Outdir: " . $self->outdir . " not found or not a directory",
            );
    }
    return @errors;
}

sub get_annotation_db {
    my $self = shift;
    my $cancer_annotation_db;
    if ($self->clinseq_build) {
        $cancer_annotation_db = $self->clinseq_build->cancer_annotation_db;
    }
    elsif ($self->cancer_annotation_db) {
        $cancer_annotation_db = $self->cancer_annotation_db;
        $self->status_message("Using cancer annotation db: " . $self->cancer_annotation_db->id);
    }
    else {
        die $self->error_message("Please specify cancer annotation db");
    }
    return $cancer_annotation_db;
}

sub get_microarray_builds {
    my $self = shift;
    if ($self->microarray_model_tumor and $self->microarray_model_normal) {
        return ($self->microarray_model_tumor->last_succeeded_build, $self->microarray_model_normal->last_succeeded_build);
    }
    elsif (my $clinseq_build = $self->clinseq_build) {
        if ($clinseq_build->has_microarray_build) {
            return ($clinseq_build->tumor_microarray_build, $clinseq_build->normal_microarray_build);
        }
        else {
            die $self->error_message(
                "Please specify a clinseq model with either an exome-sv [or] wgs-sv model with microarray builds");
        }
    }
    else {
        die $self->error_message("Please specify one clinseq model [or] two microarray models as input!");
    }
}

sub get_copynumber_files {
    my $self = shift;
    my $copynumber_tumor;
    my $copynumber_normal;
    if ($self->copynumber_tumor and $self->copynumber_normal) {
        if (-e $self->copynumber_tumor and -e $self->copynumber_normal) {
            $copynumber_tumor  = $self->copynumber_tumor;
            $copynumber_normal = $self->copynumber_normal;
        }
        else {
            die $self->error_message("Unable to find one of the copynumber files");
        }
    }
    else {
        my ($microarray_tumor, $microarray_normal) = $self->get_microarray_builds;
        if (-e $microarray_tumor->original_genotype_file_path) {
            $copynumber_tumor = $microarray_tumor->original_genotype_file_path;
        }
        else {
            die $self->error_message("Unable to find copynumber file for " . $microarray_tumor->model->name);
        }
        if (-e $microarray_normal->original_genotype_file_path) {
            $copynumber_normal = $microarray_normal->original_genotype_file_path;
        }
        else {
            die $self->error_message("Unable to find copynumber file for " . $microarray_normal->model->name);
        }
    }
    my $copynumber_tumor_copy  = $self->outdir . "/tumor.copynumber.original";
    my $copynumber_normal_copy = $self->outdir . "/normal.copynumber.original";
    Genome::Sys->copy_file($copynumber_tumor,  $copynumber_tumor_copy);
    Genome::Sys->copy_file($copynumber_normal, $copynumber_normal_copy);
    return ($copynumber_tumor_copy, $copynumber_normal_copy);
}

sub intersect_files {
    my $self = shift;
    my $f1   = shift;
    my $f2   = shift;
    my $f1_o = $f1 . ".intersected";
    my $f2_o = $f2 . ".intersected";

    open(my $f1_fh,  "<" . $f1);
    open(my $f1_ofh, ">" . $f1_o);
    open(my $f2_fh,  "<" . $f2);
    open(my $f2_ofh, ">" . $f2_o);

    my %f1;
    my %chr_pos1;
    my $is_header = 1;

    while (<$f1_fh>) {
        my $line = $_;
        if ($is_header) {
            $f1{"header"} = $line;
            $is_header = 0;
            next;
        }
        my @fields = split("\t", $line);
        my $key = $fields[0] . ":" . $fields[1];
        $chr_pos1{$key} = 1;
        $f1{$key}       = $line;
    }

    $is_header = 1;
    while (<$f2_fh>) {
        my $line = $_;
        if ($is_header) {
            print $f2_ofh $line;
            print $f1_ofh $f1{"header"};
            $is_header = 0;
            next;
        }
        my @fields = split("\t", $line);
        my $key = $fields[0] . ":" . $fields[1];
        if ($chr_pos1{$key}) {
            print $f2_ofh $line;
            print $f1_ofh $f1{$key};
        }
    }
    return ($f1_o, $f2_o);
}

sub create_single_cnv_diff_hq_file {
    my $self    = shift;
    my $cn_f    = shift;
    my $diff_f  = shift;
    my $cnvhq_f = shift;
    #at the time of writing the headers for these files look like
    #"chromosome position alleles reference id sample_name log_r_ratio gc_score cnv_value cnv_confidence b_allele_freq allele1 allele2"
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        separator => "\t",
        input     => $cn_f,
    );
    my @headers_diff = qw/chr pos cnv_diff/;
    my $writer_diff  = Genome::Utility::IO::SeparatedValueWriter->create(
        output        => $diff_f,
        separator     => "\t",
        headers       => \@headers_diff,
        print_headers => 0,
    );
    my @headers_cnvhq = qw/CHR POS TUMOR NORMAL DIFF/;
    my $writer_cnvhq  = Genome::Utility::IO::SeparatedValueWriter->create(
        output        => $cnvhq_f,
        separator     => "\t",
        headers       => \@headers_cnvhq,
        print_headers => 1,
    );
    unless ($reader) {
        die $self->error_message('Unable to create SeparatedValueReader for ' . $cn_f);
    }
    while (my $data = $reader->next) {
        if (   $data->{chromosome} eq "chr"
            or $data->{log_r_ratio} eq "NaN"
            or ($self->test and $data->{chromosome} ne "6"))
        {
            next;
        }
        my $diff_data;
        my $cnvhq_data;
        $diff_data->{chr}  = $data->{chromosome};
        $cnvhq_data->{CHR} = $data->{chromosome};
        $diff_data->{pos}  = $data->{position};
        $cnvhq_data->{POS} = $data->{position};
        my $cnv_ratio = $data->{log_r_ratio};
        #copynumber ~ 2^(log_r_ratio + 1)
        $cnvhq_data->{TUMOR}   = $data->{log_r_ratio};
        $cnvhq_data->{NORMAL}  = 0;
        $diff_data->{cnv_diff} = $cnv_ratio;                                   #segment the LRR
        $cnvhq_data->{DIFF}    = $cnvhq_data->{TUMOR} - $cnvhq_data->{NORMAL};
        $writer_diff->write_one($diff_data);
        $writer_cnvhq->write_one($cnvhq_data);
    }
}

sub create_somatic_cnv_diff_hq_file {
    my $self        = shift;
    my $tumor_cn_f  = shift;
    my $normal_cn_f = shift;
    my $diff_f      = shift;
    my $cnvhq_f     = shift;
    #at the time of writing the headers for these files look like "chromosome position alleles reference id sample_name log_r_ratio gc_score cnv_value cnv_confidence allele1 allele2"
    my $reader_t = Genome::Utility::IO::SeparatedValueReader->create(
        separator => "\t",
        input     => $tumor_cn_f,
    );
    my $reader_n = Genome::Utility::IO::SeparatedValueReader->create(
        separator => "\t",
        input     => $normal_cn_f,
    );
    my @headers_diff = qw/chr pos cnv_diff/;
    my $writer_diff  = Genome::Utility::IO::SeparatedValueWriter->create(
        output        => $diff_f,
        separator     => "\t",
        headers       => \@headers_diff,
        print_headers => 0,
    );
    my @headers_cnvhq = qw/CHR POS TUMOR NORMAL DIFF/;
    my $writer_cnvhq  = Genome::Utility::IO::SeparatedValueWriter->create(
        output        => $cnvhq_f,
        separator     => "\t",
        headers       => \@headers_cnvhq,
        print_headers => 1,
    );
    unless ($reader_t and $reader_n) {
        die $self->error_message('Unable to create SeparatedValueReader for ' . $tumor_cn_f . ' and ' . $normal_cn_f);
    }
    while (my $data_t = $reader_t->next and my $data_n = $reader_n->next) {
        if (   $data_n->{chromosome} ne $data_t->{chromosome}
            or $data_n->{position} ne $data_t->{position})
        {
            $self->warning_message(
                'Mismatch in chr:pos between ' . $tumor_cn_f . ' and ' . $normal_cn_f . '. Taking intersection.');
            ($tumor_cn_f, $normal_cn_f) = $self->intersect_files($tumor_cn_f, $normal_cn_f);
            Genome::Sys->shellcmd(cmd => "rm -f $diff_f $cnvhq_f");
            $self->create_somatic_cnv_diff_hq_file($tumor_cn_f, $normal_cn_f, $diff_f, $cnvhq_f);
            return;
        }
        if ($data_n->{chromosome} eq "chr" or $data_n->{log_r_ratio} eq "NaN" or $data_t->{log_r_ratio} eq "NaN") {
            next;
        }
        my $diff_data;
        my $cnvhq_data;
        $diff_data->{chr}  = $data_n->{chromosome};
        $cnvhq_data->{CHR} = $data_n->{chromosome};
        $diff_data->{pos}  = $data_n->{position};
        $cnvhq_data->{POS} = $data_n->{position};
        #copynumber ~ 2^(log_r_ratio + 1)
        my $log_cn_normal = $data_n->{log_r_ratio};
        my $log_cn_tumor  = $data_t->{log_r_ratio};
        my $cnv_ratio     = $log_cn_tumor - $log_cn_normal;
        $cnvhq_data->{TUMOR}   = 2.0**($data_t->{log_r_ratio} + 1.0);
        $cnvhq_data->{NORMAL}  = 2.0**($data_n->{log_r_ratio} + 1.0);
        $diff_data->{cnv_diff} = $cnv_ratio;                                   #segment the ratio of LRR
        $cnvhq_data->{DIFF}    = $cnvhq_data->{TUMOR} - $cnvhq_data->{NORMAL};

        if (not $self->test) {
            $writer_diff->write_one($diff_data);
            $writer_cnvhq->write_one($cnvhq_data);
        }
        elsif ($data_n->{chromosome} eq 6) {                                   #for test, use only chr6
            $writer_diff->write_one($diff_data);
            $writer_cnvhq->write_one($cnvhq_data);
        }
    }
}

#Get the .original file from a microarray build
sub get_copynumber_file {
    my $self             = shift;
    my $microarray_model = shift;
    my $build            = $microarray_model->last_succeeded_build;
    unless ($build) {
        die $self->error_message("No succeeded build found for" . $microarray_model->id);
    }
    my $dd      = $build->data_directory;
    my $cn_file = $build->original_genotype_file_path;
    if (-e $cn_file) {
        return $cn_file;
    }
    else {
        die $self->error_message("*.original file not found in $dd");
    }
}

sub run_cbs {
    my $self = shift;
    my ($tumor_copynumber, $normal_copynumber, $cbs_op) = @_;
    my $cnv_diff_file = $self->outdir . "/cnvs.diff";
    my $cnv_hq_file   = $self->outdir . "/cnvs.hq";
    my $cbs;
    if ($self->microarray_model_single) {
        my $single_copynumber = $self->get_copynumber_file($self->microarray_model_single);
        $self->status_message("Using microarray file $single_copynumber");
        $cbs = Genome::Model::Tools::CopyNumber::Cbs->create(
            microarray_original_file => $single_copynumber,
            output_file              => $cbs_op
        );
        $self->create_single_cnv_diff_hq_file($single_copynumber, $cnv_diff_file, $cnv_hq_file);
    }
    else {
        $self->create_somatic_cnv_diff_hq_file($tumor_copynumber, $normal_copynumber, $cnv_diff_file, $cnv_hq_file);
        $cbs = Genome::Model::Tools::CopyNumber::Cbs->create(
            array_file  => $cnv_diff_file,
            output_file => $cbs_op
        );
    }
    $cbs->execute();
}

sub run_cnview {
    my $self                 = shift;
    my $tumor_copynumber     = shift;
    my $normal_copynumber    = shift;
    my $cbs_op               = shift;
    my $cancer_annotation_db = shift;
    my $min_cnv_diff         = $self->min_cnv_diff;

    if ($min_cnv_diff < 0) {
        $min_cnv_diff = $min_cnv_diff * -1;  #use the absolute value
    }

    my $cnv_hq_file = $self->outdir . "/cnvs.hq";
    #Create cnvhmm file
    my $cnv_hmm_file = $cbs_op . ".cnvhmm";
    #print only cnv segments with atleast five snp markers.
    my ($cn_diff, $cn1, $cn2);
    if ($self->microarray_model_single) {
        $cn_diff = '$5';
        $cn1     = '$5';
        $cn2     = '0';
    }
    else {
        $cn_diff = '(2 ^ $5 * 2) - 2';
        $cn1     = 'cn_diff + 2';
        $cn2     = '2';
    }
    my $make_hmmfile_cmd =
          'awk \'{ size = $3-$2; nmarkers=$4; cn_diff = ' 
        . $cn_diff
        . '; event = "NA"; if(cn_diff>0) { event = "Gain" } else if(cn_diff<0) { event = "Loss" } cn1 = '
        . $cn1
        . '; cn2 = '
        . $cn2
        . '; if((event == "Gain" || event == "Loss") && (cn_diff > '
        . $min_cnv_diff
        . ' || cn_diff < -1 * '
        . $min_cnv_diff
        . ' ) && $4 >=5) print $1"\t"$2"\t"$3"\t"size"\t"nmarkers"\t"cn1"\t"cn1"\t"cn2"\t"cn2"\tNA\t"event; } \' '
        . $cbs_op . ' > '
        . $cnv_hmm_file;
    Genome::Sys->shellcmd(cmd => $make_hmmfile_cmd);

    #For each list of gene symbols, run the CNView analysis
    my @cnv_symbols;
    if ($self->test) {
        @cnv_symbols = qw(Kinase_dGene);
    }
    else {
        @cnv_symbols = qw (Kinase_dGene CancerGeneCensusPlus_Sanger AntineoplasticTargets_DrugBank All);
    }
    my $gene_symbol_dir = $cancer_annotation_db->data_directory . "/GeneSymbolLists/";

    #Run only chr6 for test-cases
    my $plot_chr = "ALL";
    if ($self->test) {
        $plot_chr = "6";
    }

    foreach my $symbol (@cnv_symbols) {
        my $symbol_outdir = $self->outdir;
        my %params        = (
            annotation_build     => $self->annotation_build_id,
            cnv_file             => $cnv_hq_file,
            segments_file        => $cnv_hmm_file,
            output_dir           => $symbol_outdir,
            name                 => $symbol,
            cancer_annotation_db => $cancer_annotation_db,
            window_size          => 0,
            verbose              => 1,
            chr                  => $plot_chr
        );
        if ($symbol ne 'All') {
            $params{gene_targets_file} = "$gene_symbol_dir/$symbol" . ".txt";
        }
        my $cnview_cmd = Genome::Model::Tools::CopyNumber::CnView->create(%params);
        $cnview_cmd->execute;
    }
}

sub execute {
    my $self                 = shift;
    my $copynumber_tumor     = "NA";
    my $copynumber_normal    = "NA";
    my $cbs_op               = $self->outdir . "/cnvs.diff.cbs";
    my $cancer_annotation_db = $self->get_annotation_db();
    unless ($self->microarray_model_single) {
        ($copynumber_tumor, $copynumber_normal) = $self->get_copynumber_files();
    }
    $self->run_cbs($copynumber_tumor, $copynumber_normal, $cbs_op);
    $self->run_cnview($copynumber_tumor, $copynumber_normal, $cbs_op, $cancer_annotation_db);
    return 1;
}

1;
