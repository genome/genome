package Genome::Model::Tools::Relationship::PolymuttPipeline;

use strict;
use warnings;
use File::Basename;
use Genome;
use Genome::Info::IUB;
use Workflow;
use Workflow::Simple;

use List::MoreUtils qw/ uniq /;
class Genome::Model::Tools::Relationship::PolymuttPipeline {
    is => 'Command',
    has => [
    output_dir => {
        is => 'Text',
        is_optional=>0,
        doc => 'the directory where you want results stored',
    },
    ped_file=> {
        is=>'Text',
        is_optional=>0,
        doc=>'ped file for individuals you have supplied',
    },
    model_group=> {
        is=>'Number',
        is_optional=>0,
    },
    _ref_fasta=> {
        is_optional=>1,
    },
    ],
};

sub help_synopsis {
    return <<EOS

EOS
}


sub execute {
    my $self = shift;
    $DB::single=1;
    unless(-e $self->ped_file) {
        my $ped_file = $self->ped_file;
        $self->error_message("supplied ped file does not exist: $ped_file");
        return;
    }
    unless(-d $self->output_dir) {
        my $dir = Genome::Sys->create_directory($self->output_dir);
        unless(-d $dir) {
            $self->error_message("Unable to create output directory");
            return;
        }
    }
    my $mg = Genome::ModelGroup->get($self->model_group);
    my @models = sort $mg->models;
    unless(@models) { 
        $self->error_message("Unable to get models from model group");
    }
    #self->verify_models_have_succeeded_builds_and_bams FIXME 
    $self->_ref_fasta($models[0]->reference_sequence_build->full_consensus_path("fa"));
    my @glfs = $self->generate_glfs(@models);
#    my @glfs = glob($self->output_dir . "/*glf");  
    my $dat_file = $self->generate_dat();
    my $glf_index = $self->generate_glfindex(@glfs);
    $self->run_polymutt($dat_file, $glf_index);

}
sub run_polymutt {
    my($self, $dat_file, $glf_index) = @_;
    my $ped_file = $self->ped_file;
    my %inputs;
    $inputs{dat_file}=$dat_file;
    $inputs{ped_file}=$ped_file;
    $inputs{denovo}=1;
    $inputs{glf_index}=$glf_index;
    chomp(my $family_id = `head -n 1 $ped_file | cut -f 1`);
    $inputs{output_denovo} = $self->output_dir . "/$family_id.denovo.vcf";
    $inputs{output_standard} = $self->output_dir . "/$family_id.standard.vcf";
    my $workflow = Workflow::Model->create(
        name=> "Run polymutt standard and denov",
        input_properties => [
        'dat_file',
        'glf_index',
        'ped_file',
        'output_denovo',
        'output_standard',
        'denovo',
        ],
        output_properties => [
        'output',
        ],
    );
    my $denovo_op = $workflow->add_operation(
        name=>"denovo polymutt",
        operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Relationship::RunPolymutt"),
    );
    my $standard_op = $workflow->add_operation(
        name=>"standard polymutt",
        operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Relationship::RunPolymutt"),
    );
    for my $op ($denovo_op, $standard_op) {
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"dat_file",
            right_operation=>$op,
            right_property=>"dat_file",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"glf_index",
            right_operation=>$op,
            right_property=>"glf_index",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"ped_file",
            right_operation=>$op,
            right_property=>"ped_file",
        );
        $workflow->add_link(
            left_operation=>$op,
            left_property=>"output_vcf",
            right_operation=>$workflow->get_output_connector,
            right_property=>"output",
        );
    }
    $workflow->add_link(
        left_operation=>$workflow->get_input_connector,
        left_property=>"output_denovo",
        right_operation=>$denovo_op,
        right_property=>"output_vcf",
    );
    $workflow->add_link(
        left_operation=>$workflow->get_input_connector,
        left_property=>"output_standard",
        right_operation=>$standard_op,
        right_property=>"output_vcf",
    );
    $workflow->add_link(
        left_operation=>$workflow->get_input_connector,
        left_property=>"denovo",
        right_operation=>$denovo_op,
        right_property=>"denovo",
    );
    my @errors = $workflow->validate;
    $workflow->log_dir($self->output_dir);
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }
    $self->debug_message("Now launching 2 polymutt jobs");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %inputs);
    unless($result) {
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        $self->error_message("parallel polymutt did not return correctly.");
        die;
    }

}




sub generate_dat {
    my $self = shift;
    my $out_file_name = $self->output_dir . "/" . "polymutt.dat";
    my $dat_fh = IO::File->new($out_file_name, ">");
    $dat_fh->print("T\tGLF_Index\n");
    $dat_fh->close;
    return $out_file_name;
}

sub generate_glfindex {
    my $self=shift;
    my @glfs = @_;
    my $glf_name = $self->output_dir ."/" . "polymutt.glfindex";
    my %sample_index = $self->parse_ped();
    my $glf_fh = IO::File->new($glf_name, ">");
    for my $glf (sort @glfs) {  #again we just assume the incoming ped file is alpha sorted, so we alpha sort
        my ($file, $path, $suffix) = fileparse($glf, ".glf");
        my $index = $sample_index{$file}; 
        unless($index) {
            $self->error_message("I am unable to match up the generated glfs with the supplied ped file. Please assist.  I was searching for a sample name match for this file: $glf\n");
            die;
        }
        $glf_fh->print("$index\t$glf\n");
    }
    $glf_fh->close;
    return $glf_name;
}

sub parse_ped {
    my $self = shift;
    my $ped_fh = Genome::Sys->open_file_for_reading($self->ped_file);
    my %sample_index;
    while(my $ped_line = $ped_fh->getline) {
        chomp($ped_line);
        my ($family_id, $individual, $mother, $father, $sex, $glf_index) = split "\t", $ped_line;
        $sample_index{"$individual"}=$glf_index; 
    }
    $DB::single=1;
    return %sample_index;
}

sub generate_glfs {
    my $self = shift;
    my @models = @_;
    my %inputs;
    my (@outputs, @inputs);
    $inputs{ref_fasta} = $self->_ref_fasta;
    for (my $i =0; $i < scalar(@models); $i++) {
        my $output_name = $self->output_dir . "/" . $models[$i]->subject->name . ".glf";
        push @outputs, $output_name;
        my $build = $models[$i]->last_succeeded_build;
        $inputs{"bam_$i"}=$build->whole_rmdup_bam_file;
        $inputs{"output_glf_$i"}=$output_name;
        push @inputs, ("bam_$i", "output_glf_$i");
    }
    my $workflow = Workflow::Model->create(
        name=> "polymutt parallel glf file creation",
        input_properties => [
        'ref_fasta',
        @inputs,
        ],
        output_properties => [
        'output',
        ],
    );
    for(my $i=0; $i< scalar(@models); $i++) {
        my $hybridview_op = $workflow->add_operation(
            name=>"glf creation $i",
            operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Samtools::HybridView"),
        );

        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"ref_fasta",
            right_operation=>$hybridview_op,
            right_property=>"ref_fasta",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"bam_$i",
            right_operation=>$hybridview_op,
            right_property=>"bam",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"output_glf_$i",
            right_operation=>$hybridview_op,
            right_property=>"output_glf",
        );
        $workflow->add_link(
            left_operation=>$hybridview_op,
            left_property=>"output_glf",
            right_operation=>$workflow->get_output_connector,
            right_property=>"output",
        );
    }
    my @errors = $workflow->validate;
    $workflow->log_dir($self->output_dir);
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }
    $self->debug_message("Now launching glf generation jobs");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %inputs);
    unless($result) {
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("parallel glf generation workflow did not return correctly.");
    }
    return @outputs;
}




1;
