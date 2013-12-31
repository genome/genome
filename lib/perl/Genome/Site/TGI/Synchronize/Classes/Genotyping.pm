package Genome::Site::TGI::Synchronize::Classes::Genotyping; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::Genotyping { # EXTERNAL 2877138689 ILLUMINA (INTERNAL) 2883765759
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    table_name => <<SQL
    (
        select g.seq_id id, g.status status,
         g.organism_sample_id sample_id, s.full_name sample_name,
	     lower(p.name) sequencing_platform, p.chip_type chip_name, p.version version,
         'external' import_source_name
	    from external_genotyping g
        join genotyping_platform p on p.genotyping_platform_id = g.genotyping_platform_id
        join organism_sample s on s.organism_sample_id = g.organism_sample_id
        union all
        select g.seq_id id, g.status status,
         g.organism_sample_id sample_id, s.full_name sample_name,
         'infinium' sequencing_platform, 'HumanOmniExpress' chip_name, '12v1_A' version,
         'wugc' import_source_name
        from illumina_genotyping g
        join organism_sample s on s.organism_sample_id = g.organism_sample_id
        where g.status = 'pass'
    ) genotyping
SQL
    ,
    id_by => [
        id => { is => 'Text', },
    ],
    has => [
        status => { is => 'Text', },
        sample_id => { is => 'Text', },
        sample_name => { is => 'Text', },
        sequencing_platform => { is => 'Text', },
        chip_name => { is => 'Text', },
        version => { is => 'Text', },
        import_source_name => { is => 'Text', },
    ],
    has_optional_transient => [
        genotype_file => { is => 'Text', },
        library => { is => 'Genome::Library', },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub entity_name { return 'instrument data microarray'; }

sub genome_class_for_create { return 'Genome::InstrumentData::Imported'; }

sub properties_to_copy {
    return ( 'id', properties_to_keep_updated() );
}

sub properties_to_keep_updated {
    return (qw/ 
        chip_name
        version
        import_source_name
        sequencing_platform
        /);
}

sub lims_property_name_to_genome_property_name {
    my ($class, $name) = @_;
    my %lims_to_genome = (
        platform_name => 'sequencing_platform',
    );
    return $lims_to_genome{$name} if exists $lims_to_genome{$name};
    return $name;
}

sub create_in_genome {
    my $self = shift;

    my $genotype_file = $self->genotype_file;
    if ( not $genotype_file or not -s $genotype_file ) {
        Carp::confess('No genotype file set!');
    }

    my $library_name = $self->sample_name.'-microarraylib';
    my ($library) = Genome::Library->get(name => $library_name, sample_id => $self->sample_id);
    if ( not $library ) {
        $library = Genome::Library->create(name => $library_name, sample_id => $self->sample_id);
        if ( not $library ) {
            Carp::confess('Failed to create genotype microarray library for sample: '.$self->sample_id);
        }
    }
    $self->library($library);

    my %params = $self->params_for_create_in_genome;
    return if not %params;

    my $genome_class = $self->genome_class_for_create;
    my $genome_object = $genome_class->create(%params); 
    return if not $genome_object;

    my $new_genotype_file = eval{ Genome::InstrumentData::Microarray->update_genotype_file($genome_object, $genotype_file); };
    if ( not $new_genotype_file ) {
        Carp::confess "$@\nFailed to update genotype_file: $genotype_file on instrument data: ".$genome_object->id;
    }

    return $genome_object;
}

sub params_for_create_in_genome {
    my $self = shift;

    my %params = $self->SUPER::params_for_create_in_genome;
    return if not %params;

    $params{import_format} = 'genotype file';
    $params{library} = $self->library;

    return %params;
}

1;

