package Genome::ModelGroup::Command::CreateSimulationPopulation;

use strict;
use warnings;

use Genome;
use Workflow;
use Workflow::Simple;
use File::Basename;
class Genome::ModelGroup::Command::CreateSimulationPopulation {
    is => 'Genome::Command::Base',
    has_input => [
    output_directory => {
        is => 'Text',
        is_optional=>0,
        doc => 'the directory where you want results stored',
    },
    coverage_target => {
        is => 'Text',
        is_optional => 1,
        doc => 'How much coverage should each bam end up with',
    },
    region => {
        is => 'Text', 
        is_optional => 1,
        doc=>"chr:start-stop format",
    },
    ref_fasta => {
        example_values =>["/gscmnt/gc4096/info/model_data/2741951221/build101947881/all_sequences.fa"],
        doc=>"build36 default",
    },
    num_cases => {
        is => 'Text',
        is_optional=>1,
        default=>'2'
    },
    num_controls => {
        is => 'Number', 
        is_optional=>1,
        default=>2,
    },
    disease_snps => {
        is =>'Text',
        is_optional=>0,
        doc => "see generate haplotypes doc on how to specify this (basically, <pos>, <which allele is the disease allele>, <het rel risk>, <hom rel risk>",
    },
    exome_regions=> {
        is =>'Text',
        is_optional=>0,
        doc => "a bed file to use to simulate exome sequencing.  supply a sureselect/nimblegen feature list for build36, it will emulate wingspan for you",
    },
    model_group_name => {
        is =>'Text',
        is_optional=>0,
        doc =>"what to name the resulting model group.  Don't make this too long as its integrated into the lib name which has a 64 character limit, and some stuff gets appended on",
    }
    ],
    has_transient_optional=> [
    _cases_bed_dir => {
        is=>'Text',
        calculate_from => ['output_directory'],
        calculate => q| return $output_directory . "/cases/case"|,
    },
    _controls_bed_dir => {
        is=>'Text',
        calculate_from => ['output_directory'],
        calculate => q| return $output_directory . "/controls/control";|,
    },

    ],
    doc => '', #TODO
};

sub help_synopsis {
    return <<EOS
genome model-group create-cross-sample-vcf --model-group=1745 --output-dir=/foo/bar/

EOS
}


sub execute {
    my $self=shift;
    my $output_directory = $self->output_directory;
     $DB::single=1;
    #Check for output directory
    unless(-d $output_directory) {
        $self->error_message("Unable to find output directory: " . $output_directory);
        return;
    }   

    my %inputs;
    $inputs{output_directory}=$output_directory;
    $inputs{disease_snps}=$self->disease_snps;
    $inputs{coverage_target}=$self->coverage_target;
    $inputs{number_of_cases}=$self->num_cases;
    $inputs{number_of_controls}=$self->num_controls;
    my $case_bed_dir = $output_directory . "/cases/";
    my $control_bed_dir = $output_directory . "/controls/";
    unless(-e $case_bed_dir) {
        mkdir($case_bed_dir);
    }
    unless(-e $control_bed_dir) {
        mkdir($control_bed_dir);
    }
        $case_bed_dir .= "case";
        $control_bed_dir .= "control";
    $inputs{cases_bed_dir}= $case_bed_dir;
    $inputs{controls_bed_dir}=$control_bed_dir;
    my $workflow = $self->_generate_workflow($output_directory, \%inputs);
    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }


    $self->status_message("Now launching the hapgen2 to generate haplotypes, then converting them to bed");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %inputs);
    unless($result){
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("hapgen2 workflow did not return correctly.");
    }
    $self->generate_fastas($case_bed_dir);
    $self->generate_fastas($control_bed_dir);
    my @case_bams = $self->generate_bams($case_bed_dir, $self->coverage_target);
    my @control_bams = $self->generate_bams($control_bed_dir, $self->coverage_target);


    #generate samples,import bams, make models
    my $mg_name = $self->model_group_name;
    $mg_name =~ s/ /_/g;
    my @finished_models;
    for my $bam (@case_bams, @control_bams) {
        my ($bam_name, $path, $bam_suffix) = fileparse($bam, ".bam");
        my $sample_name = "TEST-$mg_name-$bam_name";
        my $sample = Genome::Sample->get_or_create(name=>"$sample_name"); 
        my $import_cmd = Genome::InstrumentData::Command::Import::Bam->create(
            target_region=> 'none',
            original_data_path=> $bam,
            sample => $sample_name,
            create_library=>1,
        );
        eval {
            $import_cmd->execute();
        };
        if($@) {
            $self->status_message("Import Command Failure:$@");
            $self->cleanup_cruft($sample_name, $import_cmd->import_allocation_id);
        }
        else {
            $self->status_message("Import Command Success!");
        }
        my $model = Genome::Model::ReferenceAlignment->create(
            reference_sequence_build_id=>101947881,
            dbsnp_build_id=>106227442,
            processing_profile_id=>2638105,
            name=>"$mg_name-$bam_name-germline-pipeline",
            subject=>$sample,
        );
        if(!$model) {
            $self->status_message("Model Creation Command Failure: $@");
            $self->cleanup_cruft($sample_name, $import_cmd->import_allocation_id);
        }
        $DB::single=1;
        my $id_cmd = Genome::Model::Command::InstrumentData::Assign::Expression->create(
            model=>$model,
            instrument_data=>[$import_cmd->_inst_data],
        );
        eval{ 
            $id_cmd->execute();
        };     
        if($@) {
            $self->status_message("InstrumentData Import Command Failure: $@");
            $self->cleanup_cruft($sample_name, $import_cmd->import_allocation_id);
        }
        push @finished_models, $model;


    }
    $self->status_message("Adding these models: " . join " ", @finished_models . "\n");
    my $mg = Genome::ModelGroup->create(
        name=> $mg_name,
    );
    $mg->assign_models(@finished_models);

    #start modelgroup

    return 1;


}

sub cleanup_cruft {
    my $self = shift;
    my $sample_name = shift;
    my $allocation_id = shift;
    $self->status_message("###Cruft cleanup###:");
    $self->status_message(qq|perl -e "use Genome; Genome::Sample->create(name=>'$sample_name'); UR::Context->commit()"|);
    if($allocation_id) {
    $self->status_message(qq|genome model deallocate $allocation_id|);
}
else {
    $self->status_message("No allocation id passed, no cleanup necessary?");
 }
 $self->status_message("###################");
 return;

}


sub generate_bams {
    my $self = shift;
    my $fasta_dir = shift;
    my $cov_level = shift;
    my @files = glob("$fasta_dir*.fasta");
    $self->error_message("Found: " . scalar(@files) . " files in $fasta_dir to make bams from");
    my %inputs;
    $inputs{diploid_fasta} = \@files;
    $inputs{coverage_level}=$cov_level;
    my $op = Workflow::Operation->create( 
        name=>'Generate Bams Parallel',
        operation_type=> Workflow::OperationType::Command->get('Genome::Model::Tools::Simulation::GenerateReads'),
    );
    $op->parallel_by('diploid_fasta');
    my $output= Workflow::Simple::run_workflow_lsf($op, %inputs);
    return glob("$fasta_dir*.bam");
}





sub generate_fastas {
    my $self = shift;
    my $bed_dir=shift;
    my @files = glob("$bed_dir*.bed");
    my %inputs ; 
    $inputs{ref_fasta} = $self->ref_fasta;
    $inputs{region} = $self->region; #FIXME: pass this through when we want to support more or less than all of 22
    $inputs{limit_regions}=$self->exome_regions;
    $inputs{mutation_bed}=\@files;
    my $op = Workflow::Operation->create(
        name=> 'AlterRefSeq Parallel',
        operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Simulation::AlterReferenceSequence'),
    );
    $op->parallel_by('mutation_bed');
    my $output = Workflow::Simple::run_workflow_lsf( $op, %inputs);
}


sub _generate_workflow { 
    my $self = shift;
    my $output_dir = shift;
    my $input_ref = shift;
    my %inputs = %$input_ref;

    my $workflow = Workflow::Model->create(
        name => 'Human Genetics ModelGroup Simulator',
        input_properties => [
        'output_directory',
        'disease_snps',
        'coverage_target',
        'number_of_cases',
        'number_of_controls',
        'cases_bed_dir',
        'controls_bed_dir'
        ],
        output_properties => [
        'output',
        ],
    );
    $workflow->log_dir($output_dir);
    my $gen_haplotypes_operation = $self->_add_haplotype_generation(\$workflow);

    my $gen_case_beds = $self->hapgen_to_bed("cases", \$workflow, $gen_haplotypes_operation);
    my $gen_control_beds = $self->hapgen_to_bed("controls", \$workflow, $gen_haplotypes_operation);

    $workflow->add_link(
        left_operation => $gen_case_beds,
        left_property => "output_prefix",
        right_operation => $workflow->get_output_connector,
        right_property => "output",
    );
    $workflow->add_link(
        left_operation => $gen_control_beds,
        left_property => "output_prefix",
        right_operation => $workflow->get_output_connector,
        right_property => "output",
    );

    return $workflow;
}

sub hapgen_to_bed { 
    my $self = shift;
    my $type = shift;
    my $workflow = shift;
    $workflow = $$workflow;
    my $prev_operation= shift;
    my $hapgen2bed_op = $workflow->add_operation(
        name=>"Hapgen To Bed $type",
        operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Simulation::HapgenToBed"),
    );
    $workflow->add_link(
        left_operation=>$prev_operation,
        left_property=>"_$type\_file",
        right_operation=> $hapgen2bed_op,
        right_property=>"input_haps_file",
    );
    $workflow->add_link(
        left_operation=>$prev_operation,
        left_property=>"_legend\_file",
        right_operation=> $hapgen2bed_op,
        right_property=>"input_legend_file",
    );
    $workflow->add_link(
        left_operation=>$workflow->get_input_connector,
        left_property=>"$type\_bed_dir",
        right_operation=>$hapgen2bed_op,
        right_property=>'output_prefix',
    );
    return $hapgen2bed_op;
}



sub _add_haplotype_generation {
    my $self = shift;
    my $workflow = shift;
    $workflow = $$workflow;

    my $generate_haplotypes_operation = $workflow->add_operation( 
        name=> "Generate Haplotypes", 
        operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Simulation::GenerateHaplotypes")
    );
    #link disease-snps, number-of-cases, number-of-controls, output_directory
    for my $property ('disease_snps', 'number_of_cases', 'number_of_controls', 'output_directory') {
        $workflow->add_link(
            left_operation=> $workflow->get_input_connector,
            left_property=> $property,
            right_operation=> $generate_haplotypes_operation,
            right_property=> $property,
        );
    }
    return $generate_haplotypes_operation;
}







1;
