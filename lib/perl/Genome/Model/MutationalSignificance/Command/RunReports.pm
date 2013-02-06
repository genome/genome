package Genome::Model::MutationalSignificance::Command::RunReports;

use strict;
use warnings;
use Genome;
use IO::File;

class Genome::Model::MutationalSignificance::Command::RunReports {
    is => 'Command::V2',
    has_input => [
        output_dir => {
            is => 'Text', 
            doc => 'Output directory path',
            is_input => 1,
        },
        maf_file => {
            is => 'Text', 
            doc => 'final maf output file',
            is_input => 1,
        },
        model_list => {
            is => 'Text',
            doc => '', #TODO: write me
            is_input => 1,
        },
        annotation_build_id => {
            type => 'Text',
            doc => 'The id of the annotation build to use',
            is_optional => 1,
        },
    ],
};

sub execute {
    my $self = shift;
    # $self->create_mutation_diagrams;
    $self->create_mutation_spectrum_plots;
    return 1;
}

sub create_mutation_diagrams {
    my $self = shift;
    my $output_dir = join('/', $self->output_dir, 'mutation_diagrams');
    Genome::Sys->create_directory($output_dir);
    my $annotation_file = $self->create_mutation_diagram_annotation_file;
    my $cmd = Genome::Model::Tools::Graph::MutationDiagram->create(
        output_directory => $output_dir,
        annotation_format => 'tgi',
        annotation => $annotation_file,
        annotation_build_id => $self->annotation_build_id,
    );
    my $rv = eval{$cmd->execute()};
    if($@){
        my $error = $@;
        $self->error_message('Error running ' . $cmd->command_name . ': ' . $error);
        return;
    }

    return 1;
}

sub create_mutation_diagram_annotation_file {
    my $self = shift;
    #create a file from $self->final_maf using (1-based) columns 33-53
    my $output_dir = $self->output_dir;
    my $annotation_file = join('/', $output_dir, 'final.annotated');
    my $input_fh = Genome::Sys->open_file_for_reading($self->maf_file); #Use a separated value reader here?
    my $output_fh = Genome::Sys->open_file_for_writing($annotation_file);
    my $header = 1;
    for my $line (<$input_fh>){
        $header-- and next if $header;
        chomp $line;
        my @columns = split("\t", $line);
        $output_fh->print(join("\t", @columns[32..52]), "\n");
    }
    $input_fh->close;
    $output_fh->close;
    return $annotation_file;
}

sub create_mutation_spectrum_plots {
    my $self = shift;
    my $output_dir = join('/', $self->output_dir, 'mutation_spectrum_plots');
    Genome::Sys->create_directory($output_dir);
    my $output_file = join('/', $output_dir, 'mutation_spectrum_plot.pdf');
    my $mut_spec_file = join('/', $output_dir, 'mutation_spectrum_file.tsv');
    my $input_file = $self->create_mutation_spectrum_plots_input_file($output_dir);
    my $cmd = Genome::Model::Tools::Analysis::SummarizeMutationSpectrum->create(
        input_file => $input_file,
        output_file => $output_file,
        mut_spec_file => $mut_spec_file,
    );
    my $rv = eval{$cmd->execute()};
    if($@){
        my $error = $@;
        $self->error_message('Error running ' . $cmd->command_name . ': ' . $error);
        return;
    }
    return 1;
}

sub create_mutation_spectrum_plots_input_file {
    my $self = shift;
    my $output_dir = shift;
    my $model_list = $self->model_list;
    my %bed_label_mapping;
    my @model_ids = split(',', $model_list);
    my $bed_dir = join('/', $output_dir, 'input_bed_files');
    Genome::Sys->create_directory($bed_dir);
    for my $model_id (@model_ids){
        my $model = Genome::Model->get($model_id); 
        my ($input_file, $automatic_label) = $self->make_input_file_from_model($model, $bed_dir);
        $bed_label_mapping{$input_file} = $automatic_label;
    }
    my $input_file = join("/", $output_dir, 'labeled_bed_list.tsv');
    my $fh = Genome::Sys->open_file_for_writing($input_file);
    for my $bed_file (keys %bed_label_mapping){
        $fh->print(join("\t", $bed_label_mapping{$bed_file}, $bed_file), "\n");
    }
    $fh->close;
    return $input_file;
}

sub make_input_file_from_model {
    my $self = shift;
    my ($somatic_model, $output_dir) = @_; #requires a somatic variation model object

    my $model_id = $somatic_model->id;
    if(!defined($somatic_model)){
        $self->error_message("Invalid somatic model for $model_id!  Aborting....");
        return;
    }
    my $build=$somatic_model->last_succeeded_build;
    if(!defined($build)){
        $self->error_message("No successful build for model $model_id found! Aborting...");
        return;
    }

    my $sample_label = $build->tumor_build->model->subject->source_common_name;
    $sample_label = $somatic_model->id if(!$sample_label); #label defaults to model ID if common name cannot be found.
    my $type = $build->tumor_build->model->subject->common_name;
    $sample_label = uc("${sample_label}_${type}");

    #find the tier 1,2,3 SNV bed file
    my $dir = $build->data_directory . "/effects";
    my $tier1 = "$dir/snvs.hq.novel.tier1.v2.bed";
    my $tier2 = "$dir/snvs.hq.novel.tier2.v2.bed";
    my $tier3 = "$dir/snvs.hq.novel.tier3.v2.bed";

    #these files should all exist for each build but just in case.....
    $self->warning_message("$tier1 does not exist for $model_id.") if(!-e $tier1);
    $self->warning_message("$tier2 does not exist for $model_id.") if(!-e $tier2);
    $self->warning_message("$tier3 does not exist for $model_id.") if(!-e $tier3);

    my $bed_file = join('/', $output_dir, join('-', $build->id, "$sample_label.bed"));
    Genome::Sys->shellcmd('cmd' => "cat $tier1 $tier2 $tier3 > $bed_file");

    return ($bed_file,$sample_label);
}

1;
