package Genome::Model::RnaSeq::Command::ExpressionCvFilter;

use strict;
use warnings;

use Genome;
use Statistics::Descriptive;

class Genome::Model::RnaSeq::Command::ExpressionCvFilter {
    is => 'Command::V2',
    has => [
        models => {
            is => 'Genome::Model::RnaSeq',
            is_many => 1,
            shell_args_position => 1,
            doc => 'RNAseq models to generate expression matrix.',
        },
        fpkm_tsv_file => {
            doc => 'The output tsv file of genes passing the CV filter with FPKMvalues per sample.',
        },
        cv_tsv_file => {
            doc => 'An output tsv file with per gene FPKM metrics, like mean, standard_deviation, coefficient_of_variation(CV)',
        },
    ],
    has_optional => [
        standard_deviations => {
            doc => 'The number or ratio of standard deviations from the mean to keep.',
            default_value => '1',
        },
    ],

};

sub help_synopsis {
    return <<"EOS"
    genome model rna-seq expression-cv-filter
EOS
}

sub help_brief {
    return "Filter RNAseq FPKM values by CV.";
}

sub help_detail {
    return <<EOS
Filter FPKM values for genes across all samples by the coefficient of variation (CV) for RNAseq models.
EOS
}


sub execute {
    my $self = shift;
    my @models = $self->models;
    my @non_rna_models = grep { !$_->isa('Genome::Model::RnaSeq') } @models;
    if (@non_rna_models) {
        die('Found a non-RNAseq model: '. Data::Dumper::Dumper(@non_rna_models));
    }
    my @builds;
    my $annotation_build;
    my $reference_build;
    my %subjects;
    for my $model (@models) {
        if ( defined($subjects{$model->name}) ) {
            die('Multiple models for subject: '. $model->name);
        } else {
            $subjects{$model->name} = 1;
        }
        my $build = $model->last_succeeded_build;
        unless ($build) {
            $build = $model->latest_build;
            unless ($build) {
                die('Failed to find build for model: '. $model->id);
            }
        }
        push @builds, $build;
        my $model_reference_sequence_build = $model->reference_sequence_build;
        if ($reference_build) {
            unless ($reference_build->id eq $model_reference_sequence_build->id) {
                die('Mis-match reference sequence builds!');
            }
        } else {
            $reference_build = $model_reference_sequence_build;
        }
        my $model_annotation_build = $model->annotation_build;
        if ($annotation_build) {
            unless ($annotation_build->id eq $model_annotation_build->id) {
                die('Mis-match annotation builds!');
            }
        } else {
            $annotation_build = $model_annotation_build;
        }
    }
    my @subjects = sort keys %subjects;
    my @headers = ('gene_id',@subjects);
    my $tsv_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->fpkm_tsv_file,
        separator => "\t",
        headers => \@headers,
    );
    my $cv_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->cv_tsv_file,
        separator => "\t",
        headers => ['gene_id','mean','standard_deviation','coefficient_of_variation'],
    );
    my $gtf_path = $annotation_build->annotation_file('gtf',$reference_build->id);
    my $gff_reader = Genome::Utility::IO::GffReader->create(
        input => $gtf_path,
    );
    unless ($gff_reader) {
        die('Failed to read GTF file: '. $gtf_path);
    }
    my %genes;
    while (my $data = $gff_reader->next_with_attributes_hash_ref) {
        my $attributes = delete($data->{attributes_hash_ref});
        $genes{$attributes->{gene_id}}{gene_id} = $attributes->{gene_id};
    }
    $self->status_message('There are '. scalar(keys %genes) .' genes in annotation file: '. $gtf_path);
    for my $build (@builds) {
        my $gene_fpkm_tracking = $build->data_directory .'/expression/genes.fpkm_tracking';
        unless (-e $gene_fpkm_tracking) {
            die ('Failed to find gene FPKM file: '. $gene_fpkm_tracking);
        }
        my $gene_fpkm_reader = Genome::Utility::IO::SeparatedValueReader->create(
            input => $gene_fpkm_tracking,
            separator => "\t",
        );
        my $match = 0;
        while (my $fpkm_data = $gene_fpkm_reader->next) {
            if ( defined($genes{$fpkm_data->{gene_id}}) ) {
                if ( defined($genes{$fpkm_data->{gene_id}}{$build->model->name}) ) {
                    if ($genes{$fpkm_data->{gene_id}}{$build->model->name} < $fpkm_data->{FPKM}) {
                        $genes{$fpkm_data->{gene_id}}{$build->model->name} = $fpkm_data->{FPKM};
                    }
                } else {
                    $genes{$fpkm_data->{gene_id}}{$build->model->name} = $fpkm_data->{FPKM};
                    $match++;
                }
            }
        }
        $self->status_message('There are '. $match .' matching genes in FPKM file: '. $gene_fpkm_tracking);
    }
    my $cv_stats = Statistics::Descriptive::Sparse->new();
    for my $gene (sort keys %genes) {
        my %data = %{$genes{$gene}};
        # Depending on the mode cufflinks was run, there may not be an entry for all genes in every FPKM file,  stick to reference only mode
        if (scalar(keys %data) == (scalar(@headers))) {
            my @values = map { $data{$_} } @subjects;
            my $stat = Statistics::Descriptive::Sparse->new();
            $stat->add_data(@values);
            # Coefficient of Variation
            my $cv = 0;
            if ($stat->mean) {
                $cv = $stat->standard_deviation / $stat->mean;
                $cv_stats->add_data($cv);
            }
            my %cv_data = (
                gene_id => $gene,
                mean => $stat->mean,
                standard_deviation => $stat->standard_deviation,
                coefficient_of_variation => $cv,
            );
            $cv_writer->write_one(\%cv_data);
        } else {
            #is there a minimum number of samples(90%) that is required....
        }
    }
    my $min_cv_filter = $cv_stats->mean - ($cv_stats->standard_deviation * $self->standard_deviations);
    my $max_cv_filter = $cv_stats->mean + ($cv_stats->standard_deviation * $self->standard_deviations);
    $self->status_message('The mean coefficient of variation is: '. $cv_stats->mean);
    $self->status_message('The standard deviation from the mean coefficient of variation is: '. $cv_stats->standard_deviation);
    $self->status_message('The minimum coefficient of variation is: '. $min_cv_filter);
    $self->status_message('The maximum coefficient of variation is: '. $max_cv_filter);
    for my $gene (sort keys %genes) {
        my %data = %{$genes{$gene}};
        # Depending on the mode cufflinks was run, there may not be an entry for all genes in every FPKM file,  stick to reference only mode
        if (scalar(keys %data) == (scalar(@headers))) {
            my @values = map { $data{$_} } @subjects;
            my $stat = Statistics::Descriptive::Sparse->new();
            $stat->add_data(@values);
            # Coefficient of Variation
            if ($stat->mean) {
                my $cv = $stat->standard_deviation / $stat->mean;
                if ( ($cv >= $min_cv_filter) && ($cv <= $max_cv_filter) ){
                    $tsv_writer->write_one(\%data);
                }
            }
        }
    }
    return 1;
}
