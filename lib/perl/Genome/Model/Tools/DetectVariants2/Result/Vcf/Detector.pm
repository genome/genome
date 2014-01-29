package Genome::Model::Tools::DetectVariants2::Result::Vcf::Detector;

use strict;
use warnings;

use Genome;
use File::Copy;
use Sys::Hostname;
use File::stat qw();

class Genome::Model::Tools::DetectVariants2::Result::Vcf::Detector {
    is  => ['Genome::Model::Tools::DetectVariants2::Result::Vcf'],
};

sub _generate_vcf {
    my $self = shift;
    my $path = $self->input_directory;
    my $retval = 1;
    my %files;

    my $detector_class;
    if($self->input->can("detector_name")){
        $detector_class = $self->input->detector_name;
    } else {
        die $self->error_message("Could not call detector_name on input!");
    }
    my $convertable = 0;
    for my $variant_type ("snvs","indels"){
        my $variant_file = $path."/".$variant_type.".hq";
        unless(-e $variant_file){
            next;
        }

        my $vcf_module = $self->conversion_class_name($detector_class,$variant_type);
        if($vcf_module){
            $convertable++;
            $self->debug_message("Generating Vcf");
            $self->debug_message("executing $vcf_module on file $variant_file");
            $retval &&= $self->_run_vcf_converter($vcf_module, $variant_file,$variant_type);
        } else {
            $self->debug_message("Couldn't find a working vcf converter for $detector_class $variant_type");
            next;
        }

        if($variant_type eq 'indels'){
            $self->_run_vcf_indel_normalizer();
        }
    }
    unless($convertable>0){
        die $self->error_message("Was not able to convert any outputs to VCF");
    }

    return $retval;
}

sub _run_vcf_converter {
    my $self = shift;
    my $converter = shift;
    my $input_file = shift;
    my $type = shift;
    my $detector_result = $self->input;
    my $dirname = $self->output_dir;
    unless($dirname){
        die $self->error_message("Could not get dirname!");
    }
    my $output_file = $dirname . '/'.$type.'.vcf.gz';
    my $aligned_reads_sample = $self->aligned_reads_sample;
    my $control_aligned_reads_sample = $self->control_aligned_reads_sample;
    my $reference_sequence_build = Genome::Model::Build->get($detector_result->reference_build_id);
    my %params = (
        input_file => $input_file,
        output_file => $output_file,
        aligned_reads_sample => $aligned_reads_sample,
        sequencing_center => 'WUSTL',
        reference_sequence_build => $reference_sequence_build,
    );
    $params{control_aligned_reads_sample} = $control_aligned_reads_sample if $control_aligned_reads_sample;

    my $command = $converter->create(%params);
    unless($command->execute) {
        $self->error_message('Failed to convert ' . $input_file . ' to the standard format.');
        return;
    }
    return 1;
}

sub _run_vcf_indel_normalizer {
    my $self = shift;
    my $input_file = join('/', $self->output_dir, 'indels.vcf.gz'); #TODO: this is probably not the correct way to do this
    unless (-e $input_file){
        die $self->error_message("No file $input_file to normalize");
    }
    my $dir_name = $self->output_dir;
    my $output_file = $dir_name . '/normalized_indel.vcf.gz';
    my $detector_result = $self->input;
    my $reference_sequence_build = Genome::Model::Build->get($detector_result->reference_build_id);
    my $reference_fasta = $reference_sequence_build->full_consensus_path('fa');
    my %params = (
        input_file => $input_file,
        output_file => $output_file,
        reference => $reference_fasta,
        use_bgzip => 1,
        use_version => 1.7, #1.6 is trouble
    );

    my $command = Genome::Model::Tools::Joinx::VcfNormalizeIndels->create(%params);
    unless ($command->execute) {
        die $self->error_message("Failed to normalize $input_file.");
    }

    $self->_debug_rename($output_file, $input_file);

    return 1;
}

sub _dev {
    my $path = shift;
    my $stat = File::stat::stat($path);
    return $stat->dev();
}

sub _debug_rename {
    my $self = shift;
    my ($current, $destination) = @_;

    if (_dev($current) ne _dev($destination)) {
        $self->debug_message('Cross-device renames are not expected but has been detected:');
        $self->debug_message('current: ' . $current);
        $self->debug_message('destination: ' . $destination);
    }

    my $m = 'Failed to replace vcf file with normalized vcf file using %s: %s';

    my $rename_success = Genome::Sys::retry(
        delay    => 30,
        tries  => 3,
        callback => sub {
            my $rv = rename($current, $destination);
            my $error = $!;
            if ($rv) {
                return $rv
            } else {
                $self->warning_message(sprintf($m, 'rename', $error));
                return;
            }
        },
    );

    unless ($rename_success) {
        # we'll try one more time with move() since we want builds to succeed
        unless (File::Copy::move($current, $destination)) {
            my $move_error = $!;

            # if move() fails we may have a partial destination file
            if (-f $destination) {
                unless (unlink $destination) {
                    $self->warning_message('Failed to cleanup partial destination file after failed move(): ' . $!);
                }
            }

            die $self->error_message(sprintf($m, 'File::Copy::move', $move_error));
        }
    }
}

sub _needs_symlinks_followed_when_syncing {
    return 0;
}

sub _working_dir_prefix {
    return "detector_vcf_results";
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub allocation_subdir_prefix {
    return "detector_vcf_results";
}

sub _combine_variants {
    die "overload this function to do work";
}

sub estimated_kb_usage {
    return 10_000_000;
}

sub _staging_disk_usage {
    return 10_000_000;
}

sub _add_as_user_of_inputs {
    my $self = shift;

    my $input = $self->input;

    return (
        $input->add_user(user => $self, label => 'uses')
    );
}

1;
