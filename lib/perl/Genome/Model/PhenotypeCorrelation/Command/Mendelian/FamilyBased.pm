package Genome::Model::PhenotypeCorrelation::Command::Mendelian::FamilyBased;

use strict;
use warnings;
use File::Basename;
use Carp 'confess';
use Workflow::Simple;
use Cwd;
class Genome::Model::PhenotypeCorrelation::Command::Mendelian::FamilyBased {
    is  => 'Genome::Model::PhenotypeCorrelation::Command::Base',
    has_optional_input => [
    multisample_vcf => {
        is=>'Text',
        is_optional=>0,
        is_output=>1,
    },  
    ped_file => {
        is=>'Text',
        is_optional=>0,
    },  
    bgzip => {
        is_optional=>1,
        default=>1,
        doc=>'set this to 0 if you prefer uncompressed',
    },  
    output_directory => {
        is => "String",
        doc => 'Where to place output files',
    },  
    dbsnp_vcf => {
        is_optional=>1,
        example_values =>["/gscmnt/gc8001/info/build_merged_alignments/detect-variants--blade10-2-5.gsc.wustl.edu-acoffman-19352-120487229/snvs.hq.vcf"]
    },  
    thousand_genomes_vcf=> {
        is_optional=>1,
        example_values =>["/gscmnt/gc6132/info/medseq/1000_genomes/downloads/2012-03-27/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf"],
          
    }
]

};

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}

sub execute {
    my $self = shift;
    my $input_vcf = $self->multisample_vcf;

#    my $annotated_vcf = $self->annotate_vcf($input_vcf, $self->dbsnp_vcf, "GMAF:dbSNPBuildID=dbSNPBuildID,per-alt");
    #    my $really_annotated_vcf  = $self->annotate_vcf($annotated_vcf, $self->thousand_genomes_vcf, "ASN_AF:AMR_AF:AFR_AF:EUR_AF,per-alt");
    #currently joinx chokes on 1kg vcf so can't use it
    #   my $really_annotated_vcf = $annotated_vcf;
    $self->do_mendelian_things_in_parallel($input_vcf, $self->ped_file, $self->output_directory);
    return 1;
}

sub do_mendelian_things_in_parallel {
    my ($self, $input_vcf, $ped_file, $output_dir) = @_;
    my %inputs;
    my ($filename, $path, $suffix) = fileparse($input_vcf, (".vcf", ".vcf.gz"));
    $inputs{input_vcf} = $input_vcf;
    $inputs{dominant_output}=$output_dir . "/" . $filename . ".dominant" . $suffix;
    $inputs{recessive_output}=$output_dir . "/" . $filename . ".recessive" . $suffix;
    $inputs{compound_het_output}=$output_dir . "/" . $filename . ".compound_het" . $suffix;
    $inputs{ped_file}=$ped_file;

    my $workflow = Workflow::Model->create(
        name=>"Run Mendelian Analysis tools",
        input_properties => [ keys %inputs ],
        output_properties => [ 'output' ],
    );
    my $dom_op = $workflow->add_operation(
        name=>"Het-Dom Analysis Module",
        operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Relationship::HeterozygousDominant"),
    );
    $workflow->add_link(
        left_operation=>$workflow->get_input_connector,
        left_property=>"dominant_output",
        right_operation=>$dom_op,
        right_property=>"vcf_output",
    );
    my $recess_op = $workflow->add_operation(
        name=>"Homo-Recess Analysis Module",
        operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Relationship::HomozygousRecessive"),
    );
    $workflow->add_link(
        left_operation=>$workflow->get_input_connector,
        left_property=>"recessive_output",
        right_operation=>$recess_op,
        right_property=>"vcf_output",
    );
    my $comphet_op = $workflow->add_operation(
        name=>"Comp-Het Analysis Module",
        operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Relationship::CompoundHet"),
    );
    $workflow->add_link(
        left_operation=>$workflow->get_input_connector,
        left_property=>'compound_het_output',
        right_operation=>$comphet_op,
        right_property=>"vcf_output",
    );
    for my $op ($comphet_op, $recess_op, $dom_op) {
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>'input_vcf',
            right_operation=>$op,
            right_property=>'input_vcf',
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>'ped_file',
            right_operation=>$op,
            right_property=>'ped_file'
        );
        $workflow->add_link(
            left_operation=>$op,
            left_property=>'vcf_output',
            right_operation=>$workflow->get_output_connector,
            right_property=>'output'
        );
    };
    my @errors = $workflow->validate;
    $workflow->log_dir($self->output_directory);
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }
    $self->status_message("Now launching glf generation jobs");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %inputs);
    unless($result) {
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("parallel glf generation workflow did not return correctly.");
    }
}



1;
