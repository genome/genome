package Genome::Model::Tools::DetectVariants2::Classify::Loh;

use strict;
use warnings;

use Genome;
use Genome::Info::IUB;

class Genome::Model::Tools::DetectVariants2::Classify::Loh {
    is => 'Genome::Model::Tools::DetectVariants2::Result::Classify',
    has_input =>[
        control_result_id => {
            is => 'Text',
            doc => 'ID of the snv results considered "germline"',
        },
    ],
    has_param => [
        variant_type => {
            is_constant => 1,
            value => 'snv',
        },
    ],
    has => [
        control_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            id_by => 'control_result_id',
        },
    ],
};

sub _validate_inputs {
    my $self = shift;

    unless($self->control_result) {
        $self->error_message('No Control SNV result found.');
        return;
    }

    unless(-e (join('/', $self->control_result->output_dir, 'snvs.hq.bed'))) {
        $self->error_message('Could not find snvs file for control result.');
        return;
    }

    return $self->SUPER::_validate_inputs;
}

sub _classify_variants {
    my $self = shift;

    my $version = 2;
    my $control_variant_file = join('/', $self->control_result->output_dir, 'snvs.hq.bed');
    my $detected_snvs = join('/', $self->prior_result->output_dir, 'snvs.hq.bed');

    my $output_dir = $self->temp_staging_directory;
    my $somatic_output = $output_dir."/snvs.somatic.v".$version.".bed";
    my $loh_output = $output_dir."/snvs.loh.v".$version.".bed";

    $self->control_result->add_user(label => 'uses', user => $self);
    $self->prior_result->add_user(label => 'uses', user => $self);

    return $self->run_loh($control_variant_file, $detected_snvs, $somatic_output, $loh_output);
}

sub run_loh {
    my $self = shift;
    my ($control_variant_file,$detected_snvs,$somatic_output,$loh_output) = @_;

    my $somatic_fh = Genome::Sys->open_file_for_writing($somatic_output);
    my $loh_fh = Genome::Sys->open_file_for_writing($loh_output);

    my $normal_snp_fh = Genome::Sys->open_file_for_reading($control_variant_file);
    my $input_fh = Genome::Sys->open_file_for_reading($detected_snvs);

    #MAKE A HASH OF NORMAL SNPS!!!!!!!!!!!!!
    #Assuming that we will generally be doing this on small enough files (I hope). I suck. -- preserved in time from dlarson
    my %normal_variants;
    while(my $line = $normal_snp_fh->getline) {
        chomp $line;
        my ($chr, $start, $pos2, $ref,$var) = split /\t/, $line;
        my $var_iub;
        #Detect if ref and var columns are combined
        if($ref =~ m/\//){
            ($ref,$var_iub) = split("/", $ref);
        }
        else {
            $var_iub = $var;
        }
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
                $loh_fh->print("$line\n");
            }
            else {
                $somatic_fh->print("$line\n");
            }
        }
        else {
            $somatic_fh->print("$line\n");
        }
    }
    $input_fh->close;
    return 1;
}

sub _needs_symlinks_followed_when_syncing { 0 };
sub _working_dir_prefix { 'dv2-loh-result' };
sub resolve_allocation_disk_group_name { $ENV{GENOME_DISK_GROUP_MODELS} };

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $staged_basename = File::Basename::basename($self->temp_staging_directory);
    return join('/', 'build_merged_alignments', $self->id, 'dv2-classify-loh-' . $staged_basename);
};

sub path {
    my $self = shift;
    my ($str) = @_;

    if($str eq 'snvs.hq.bed') {
        return $self->SUPER::path('snvs.somatic.v2.bed');
    } else {
        return $self->SUPER::path(@_);
    }
}

1;
