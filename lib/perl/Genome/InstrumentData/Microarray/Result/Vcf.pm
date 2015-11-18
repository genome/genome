package Genome::InstrumentData::Microarray::Result::Vcf;

use strict;
use warnings;
use Genome;
use Sys::Hostname;

class Genome::InstrumentData::Microarray::Result::Vcf {
    is => 'Genome::SoftwareResult::StageableSimple::SingleFile',
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

sub _run {
    my $self = shift;

    my $vcf = $self->_temp_staging_file_path;
    my $params = {
        sample => $self->sample,
        variation_list_build => $self->known_sites_build,
        output => $vcf,
    };
    if ($self->filters) {
        $params->{filters} = [$self->filters];
    }
    my $rv = Genome::Model::GenotypeMicroarray::Command::ExtractToVcf->execute(%{$params});
    unless ($rv and $rv->result and -s $vcf) {
        $self->fatal_message("Could not get vcf file for ".Data::Dumper::Dumper($params));
    }
}

sub _file_name {
    my $self = shift;
    return "genotype.vcf";
}

sub vcf_path {
    my $self = shift;
    return $self->file_path;
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
    Genome::Config::get('disk_group_models');
}
1;

