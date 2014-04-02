package Genome::InstrumentData::Microarray::Result::Vcf;

use strict;
use warnings;
use Genome;
use Sys::Hostname;

class Genome::InstrumentData::Microarray::Result::Vcf {
    is => 'Genome::SoftwareResult::Stageable',
    has_input => [
        sample => {
            is => 'Genome::Sample',
        },
        known_sites_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
        },
        filters => {
            is => 'Text',
            is_many => 1,
        },
    ],
};

sub _error {
    my ($self, $msg) = @_;
    $self->error_message($msg);
    $self->delete;
    die $self->error_message;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_error("Failed to prepare staging directory") unless $self->_prepare_staging_directory;
    my $vcf = File::Spec->join($self->temp_staging_directory, $self->vcf_file_name);
    my $params = {
        sample => $self->sample,
        variation_list_build => $self->known_sites_build,
        output => $vcf,
    };
    if ($self->filters) {
        $params->{filters} = [$self->filters];
    }
    my $rv = Genome::Model::GenotypeMicroarray::Command::ExtractToVcf->execute(%{$params});
    unless ($rv and -s $vcf) {
        $self->_error("Could not get vcf file for ".Data::Dumper::Dumper($params));
    }

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;
    return $self;
}

sub vcf_file_name {
    my $self = shift;
    return "genotype.vcf";
}

sub vcf_path {
    my $self = shift;
    return File::Spec->join($self->output_dir, $self->vcf_file_name);
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;

    my $user = $ENV{'USER'};
    my $base_dir = sprintf("genotypevcf-%s-%s-%s-%s",           $hostname, $user, $$, $self->id);
    my $directory = join('/', 'build_merged_alignments',$self->id,$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}
1;

