package Genome::Model::Tools::DetectVariants2::Filter::Loh;

use warnings;
use strict;

use Genome;
use Genome::Info::IUB;

class Genome::Model::Tools::DetectVariants2::Filter::Loh{
    is => 'Genome::Model::Tools::DetectVariants2::Filter',
    doc => 'Separate LOH calls from non-LOH calls. Requires bed file input.',
    has_optional => [
        normal_snv_file  => {
            type => 'String',
            doc => 'Snv file for the LoH filter to use as a control. This will be generated if not provided.',
        },
    ],
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants2 filter loh --input-directory /path/to/sniper/outputdir --output-directory /some/output/path --control-aligned-reads-input normal.bam
EOS
}

sub help_detail {                           
    return <<EOS 
This filters out SNVs that are likely to be the result of Loss of Heterozygosity (LOH) events. The Somatic Pipeline will pass these through on its own as they are tumor variants that differ from the normal. This script defines a variant as LOH if it is homozygous, there is a heterozygous SNP at the same position in the normal and the tumor allele is one of the two alleles in the normal.

The input directory must contain snvs.hq.bed and this file must be in bed format. It will output only bed files.
EOS
}

sub _variant_type { 'snvs' };

# FIXME this should be reframed, and likely not stay a filter. For now we need to run samtools to detect snvs on the normal sample in order to replicate old behavior.
sub _filter_variants {
    my $self = shift;

    my $hq_fh = Genome::Sys->open_file_for_writing($self->_temp_staging_directory . "/snvs.hq.bed");
    my $lq_fh = Genome::Sys->open_file_for_writing($self->_temp_staging_directory . "/snvs.lq.bed");

    my $control_variant_file;
    if ($self->normal_snv_file) {
        $control_variant_file = $self->normal_snv_file;
        unless (-s $control_variant_file) {
            die $self->error_message("Normal snv file $control_variant_file does not exist or has zero size");
        }
    } else {
        $control_variant_file = $self->_generate_control_file;
    }
    my $normal_snp_fh = Genome::Sys->open_file_for_reading($control_variant_file);
    my $input_fh = Genome::Sys->open_file_for_reading($self->input_directory . "/snvs.hq.bed");

    #MAKE A HASH OF NORMAL SNPS!!!!!!!!!!!!!
    #Assuming that we will generally be doing this on small enough files (I hope). I suck. -- preserved in time from dlarson
    my %normal_variants;
    while(my $line = $normal_snp_fh->getline) {
        chomp $line;
        my ($chr, $start, $pos2, $ref_and_var) = split /\t/, $line;
        my ($ref, $var_iub) = split("/", $ref_and_var);
        #first find all heterozygous sites in normal
        next if($var_iub =~ /[ACTG]/);
        my @alleles = Genome::Info::IUB->iub_to_alleles($var_iub);
        $normal_variants{$chr}{$start} = join '',@alleles;
    }
    $normal_snp_fh->close;
    
    # Go through input variants. If a variant was called in both the input set and the control set (normal samtools calls):
    # If that variant was heterozygous in the control call and became homozygous in the input set, it is considered a loss of heterozygocity event, and goes in the LQ file
    # Otherwise it is not filtered out, and remains in the HQ output
    while(my $line = $input_fh->getline) {
        chomp $line;

        my ($chr, $start, $stop, $ref_and_iub) = split /\t/, $line;
        my ($ref, $var_iub) = split("/", $ref_and_iub);
        
        #now compare to homozygous sites in the tumor
        if ($var_iub =~ /[ACTG]/ && exists($normal_variants{$chr}{$start})) {
            if(index($normal_variants{$chr}{$start},$var_iub) > -1) {
                #then they share this allele and it is LOH
                $lq_fh->print("$line\n");
            }
            else {
                $hq_fh->print("$line\n");
            }
        }
        else {
            $hq_fh->print("$line\n");
        }
    }
    $input_fh->close;

    return 1;
}

# This is hacky. The Loh filter relies on a set of "control" snvs produced from the normal sample.
sub _generate_control_file {
    my $self = shift;

    # Temporarily hardcode, but this could be flexible later
    my $detector_version = "r599";
    my $detector_params = "";
    my $detector_name = "samtools";
    my $normal_detector = "Genome::Model::Tools::DetectVariants2::" . ucfirst($detector_name);
    my $detector_output_dir = $self->_temp_scratch_directory . "/normal_detector";

    # Detect snvs on the normal sample
    my $detector_command = $normal_detector->create(
        version => "$detector_version",
        params => "$detector_params",
        output_directory => $detector_output_dir,
        aligned_reads_input => $self->control_aligned_reads_input,
        reference_build_id => $self->reference_build_id,
        aligned_reads_sample => 'TEST',
    );

    unless($detector_command->execute) {
        die $self->error_message("Could not execute samtools to detect snvs on the normal sample");
    }

    my $filter_output = $self->_temp_scratch_directory . "/snpfilter";
    my $filter_command = Genome::Model::Tools::DetectVariants2::Filter::SnpFilter->create(
        previous_result_id => $detector_command->_result_id,
        output_directory => $filter_output,
    );

    unless($filter_command->execute) {
        die $self->error_message("Could not execute samtools to detect snvs on the normal sample");
    }

    my $normal_snv_file = "$filter_output/snvs.hq.bed";
    unless (-e $normal_snv_file) {
        die $self->error_message("Normal snv file $normal_snv_file does not exist");
    }

    my $copy_destination = $self->_temp_staging_directory . "/samtools.normal.snvs.hq.bed";
    Genome::Sys->copy_file($normal_snv_file, $copy_destination);

    return $copy_destination;
}

sub _check_file_counts {
    return 1;
}

sub _create_detector_file {
    return 1;
}

1;
