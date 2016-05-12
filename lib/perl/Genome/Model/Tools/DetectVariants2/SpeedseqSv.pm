package Genome::Model::Tools::DetectVariants2::SpeedseqSv;

use warnings;
use strict;

use Genome;

use File::Basename qw(fileparse);
use File::Spec;

class Genome::Model::Tools::DetectVariants2::SpeedseqSv {
    is => 'Genome::Model::Tools::DetectVariants2::Detector',
    has_param => [
        lsf_resource => {
            default_value => Genome::Config::get('lsf_resource_dv2_speedseq_sv'),
        },
    ],
};

# For the Parameters what I need to do is make two hashes.
# The first one will have the name of the command in the SV class and the letter calling it.  
# The second Hash will have the class string name and the value wil be either a boolean or a file location.

sub _detect_variants {
    my $self = shift;

    my @fullBam = ($self->aligned_reads_input,$self->control_aligned_reads_input);

    my $aligned_bams = join(',',@fullBam);

    my %tool_param_to_property_name = $self->get_tool_parameter_hash();

    while (my ($param_name, $property_name) = each(%tool_param_to_property_name)){
        $self->debug_message($param_name .' => '. $property_name);
    }
    my %final_cmd = ();

    my %dv2_params = $self->split_params_to_letter();

    while (my ($param_name, $property_name) = each(%tool_param_to_property_name)) {
        $final_cmd{$property_name} = $dv2_params{$param_name} if exists $dv2_params{$param_name};
    }

    %final_cmd = (
        %final_cmd,
        output_prefix => $self->_sv_staging_output,
        full_bam_file => $aligned_bams,
        version => $self->version,
        temp_directory => Genome::Sys->create_temp_directory(),
        reference_fasta => $self->reference_sequence_input,
        $self->find_file("splitters","split_read_bam_file",@fullBam),
        $self->find_file("discordants","discordant_read_bam_file",@fullBam),
    );

    my $set = Genome::Model::Tools::Speedseq::Sv->create(%final_cmd);
    $set->execute();

    # Make a symlink to the actual VCF as the DV2 Dispatcher expects an output of svs.hq rather than BED or VCF output
    my ($vcf_file) = glob($self->_sv_staging_output .'*.vcf.gz');
    if ($vcf_file) {
        symlink($vcf_file,$self->_sv_staging_output);
    }
}

sub find_file {
    my $self = shift;
    my $file = shift;
    my $value = shift;
    my @bam_dir = @_;
    my @final = ();

    foreach (@bam_dir){
        my ($editor, $dir, $suffix) = fileparse($_, '.bam');
        my $newFile = "$dir$editor.$file$suffix";
        if (!-s $newFile) {die $self->error_message("File couldn't be found: $newFile Bam Files Must be aligned with Speedseq.")};
        push (@final, $newFile);
    }
    my $combined_splits = join (',',@final);
    return (
        $value => $combined_splits,
    );
}

sub split_params_to_letter {
    my $self = shift;

    my $params = $self->params;
    
    my %params_hash = ();
    my @params = split(',',$params);
    foreach (@params){
        if ($_ =~ /:/){
            my ($letter, $value) = split(':',$_, 2);
            $letter =~ s/-//gi;
            $params_hash{$letter} = $value;
        }
        else {
            my $letter = substr($_,1,1);
            $params_hash{$letter} = 1;
        }
    }
    return %params_hash;
}

sub get_tool_parameter_hash {
    my $self = shift;
    my @meta_array = Genome::Model::Tools::Speedseq::Sv-> _tool_param_metas();

    my %library = ();

    foreach my $meta (@meta_array){
        $library{$meta->tool_param_name} = $meta->property_name; 
    }
    return %library;
}

sub has_version {
    my $self    = shift;
    my $version = shift;

    my $speedseq_path = Genome::Model::Tools::Speedseq::Base->path_for_version($version);
    if (-e $speedseq_path) {
        return 1;
    }
    else {
        return 0;
    }
}


1;
