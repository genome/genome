package Genome::InstrumentData::Command::Microarray::Extract;

use strict;
use warnings;

use Genome;

use Data::Dumper;

class Genome::InstrumentData::Command::Microarray::Extract {
    is => 'Command::V2',
    has_optional => [
        output => {
            is => 'Text',
            default_value => '-',
            doc => 'The output. Defaults to STDOUT.',
        },
        fields => {
            is => 'Text',
            is_many => 1,
            default_value => [qw/ chromosome position alleles /],
            valid_values => [qw/ chromosome position alleles id sample_id log_r_ratio gc_score cnv_value cnv_confidence allele1 allele2 /],
            doc => 'The fields to output in the genotype file.',
        },
        headers => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'Add the fields to the top of the genotype file.',
        },
        separator => {
            is => 'Text',
            default_value => 'tab',
            doc => 'Field separator of the output. Use "tab" for tab delineated.',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'The genotype instrument data to work with.',
        },
        sample => {
            is => 'Genome::Sample',
            doc => 'The sample instrument data to work with.',
        },
        use_default => {
            is => 'Boolean',
            default_value => 0,
            doc => 'If getting by sample, get the default genotype data, if available.',
        },
        use_external => {
            is => 'Boolean',
            default_value => 0,
            doc => 'If getting by sample, get the external genotype data.',
        },
        filters => {
            is => 'Text',
            is_many => 1,
            doc => "Filter genotypes. Give name and parameters, if required. Filters:\n gc_scrore => filter by min gc score (Ex: gc_score:min=0.7)\n invalid_iscan_ids => list of invalid iscan snvs compiled by Nate",
        },
        _filters => { is_transient => 1 },
        alleles => { is => 'Hash', is_transient => 1 },
        genotypes_loaded => { is => 'Number', is_transient => 1, default_value => 0, },
        genotypes_filtered => { is => 'Number', is_transient => 1, default_value => 0, },
        genotypes_kept => { is => 'Number', is_transient => 1, default_value => 0, },
        genotypes_output => { is => 'Number', is_transient => 1, default_value => 0, },

    ],
    has => [
        variation_list_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'Imported variation list build. Give id from command line. Commonly used:
                     ID          REFERENCE                   VERSION
                     106227442   dbSNP-NCBI-human-build36    130
                     106375969   dbSNP-g1k-human-build37     132',
        },
    ],
};

sub help_brief {
    return 'extract genotype data';
}

sub help_detail {
    return <<HELP;
HELP
}

sub execute {
    my $self = shift;

    my $variation_list_build = $self->variation_list_build;
    if ( not $variation_list_build ) {
        $self->error_message('No variation list build given!');
        return;
    }

    my $instrument_data = $self->_resolve_instrument_data;
    return if not $instrument_data;
    $self->status_message('Instrument data: '.$self->instrument_data->__display_name__);

    my $filters = $self->_create_filters;
    return if not $filters;

    my $output_ok = eval{ Genome::Sys->validate_file_for_writing_overwrite($self->output); };
    if ( not $output_ok ) {
        $self->error_message('Failed to validate output ('.$self->output.') for over write!');
        return;
    }

    my $genotypes = $self->_load_genotyopes;
    return if not $genotypes;

    my $annotated_genotypes = $self->_annotate_genotypes($genotypes);
    return if not $annotated_genotypes;

    my $ok = $self->_output_genotypes($annotated_genotypes);
    return if not $ok;

    if ( $self->genotypes_output == 1 ) {
    }

    return 1;
}

sub _resolve_instrument_data {
    my $self = shift;

    if ( $self->instrument_data ) {
        return 1;
    }

    my $sample = $self->sample;
    if ( not $sample ) {
        $self->error_message('No instrument data or sample given!');
        return;
    }

    if ( $self->use_default ) {
        my $default = $self->_resolve_instrument_data_from_sample_default_genotype_id($sample);
        return $default;
    }

    return $self->_resolve_instrument_data_from_library($sample);
}

sub _resolve_instrument_data_from_sample_default_genotype_id {
    my ($self, $sample) = @_;

    $self->status_message('Use default genotype data');

    if ( not $sample->default_genotype_data_id ) {
        $self->status_message('Sample ('.$sample->__display_name__.') does not have default genotype data!');
        return;
    }

    my $instrument_data = Genome::InstrumentData->get($sample->default_genotype_data_id);
    if ( not $instrument_data ) {
        $self->error_message('Failed to get default genotype data for id: '.$sample->default_genotype_data_id);
        return;
    }
    $self->instrument_data($instrument_data);

    return 1;
}

sub _resolve_instrument_data_from_library {
    my ($self, $sample) = @_;

    my $lib_name = $sample->name.'-microarraylib';
    my @libraries = Genome::Library->get(name => $lib_name);
    if (@library > 1) {
        $self->error_message("Mulitple libraries found with name $lib_name for sample ".$sample->__display_name__);
        return;
    } elsif ( ! @libraries ) {
        $self->error_message('Failed to get microarry library for sample: '.$sample->__display_name__);
        return;
    }

    my $library = $libraries[0];
    my %params = (
        library => $library,
        'import_source_name in' => ( $self->use_external )
                                    ? [qw/ BGI bgi Broad broad CSHL cshl external /]
                                    : [qw/ wugsc wugc wutgi tgi /],
    );
    my @instrument_data = Genome::InstrumentData::Imported->get(%params);
    if ( not @instrument_data ) {
        $self->error_message(
            'Failed to get instrument data for library ('.$library->__display_name__
            .') and import source ('.join(',', $params{'import_source_name in'}).')'
        );
        return;
    }
    elsif ( @instrument_data > 1 ) {
        $self->status_message(
            'Found multiple instrument data for library '.$library->__display_name__.
            ' and import source ('.join(',', @{$params{'import_source_name in'}}).
            '): '.join(',', map { $_->__display_name__ } @instrument_data).'. Using most recent.'
        );
    }
    $self->instrument_data($instrument_data[$#instrument_data]);

    return 1;
}

sub _create_filters {
    my $self = shift;

    return 1 if not $self->filters;

    $self->status_message('Filters...');

    my @filters;
    for my $filter_string ( $self->filters ) {
        $self->status_message('Filter: '.$filter_string);
        my ($name, $config) = split(':', $filter_string, 2);
        my %params;
        %params = map { split('=') } split(':', $config) if $config;
        my $filter_class = 'Genome::Model::GenotypeMicroarray::Filter::By'.Genome::Utility::Text::string_to_camel_case($name);
        my $filter = $filter_class->create(%params);
        if ( not $filter ) {
            $self->error_message("Failed to create filter for $filter_string");
            return;
        }
        push @filters, $filter;
    }

    $self->status_message('Filters...OK');

    $self->_filters(\@filters);
    return 1;
}

sub _open_output {
      my $self = shift;

    $self->status_message('Open output...');

    my $output = $self->output;
    unlink $output if -e $output;
    $self->status_message('Output file: '.$output);
    my $output_fh = eval{ Genome::Sys->open_file_for_writing($output); };
    if ( not $output_fh ) {
        $self->error_message("Failed to open output file ($output): $@");
        return;
    }

    $self->status_message('Open output...OK');

    return $output_fh;
}

sub _load_genotyopes {
    my $self = shift;

    $self->status_message('Open Genotype file...');

    my $instrument_data = $self->instrument_data;
    my $genotype_file_attr = $instrument_data->attributes(attribute_label => 'genotype_file');
    if ( not $genotype_file_attr ) {
        $self->error_message('No genotype file attribute for instrument data: '.$instrument_data->id);
        return;
    }

    my $genotype_file = $genotype_file_attr->attribute_value;
    if ( not -s $genotype_file ) {
        $self->error_message('Genotype file file does not exist! '.$genotype_file);
        return;
    }
    $self->status_message('Genotype file: '.$genotype_file);

    my $genotype_fh = eval{ Genome::Sys->open_file_for_reading($genotype_file); };
    if ( not $genotype_fh ) {
        $self->error_message("Failed to open reader for genotype file: $genotype_file): $@");
        return;
    }
    $self->status_message('Open genotype file...OK');

    my $snp_id_mapping = Genome::InstrumentData::Microarray->get_snpid_hash_for_variant_list(
        $instrument_data, $self->variation_list_build
    );
    $self->status_message(( $snp_id_mapping ? 'Got' : 'No' ).' snp id mapping...OK');

    $self->status_message('Load and filter genotypes...');
    my $header_line;
    do { $header_line = $genotype_fh->getline; } until not $header_line or $header_line =~ /,/;
    if ( not $header_line ) {
        $self->error_message('Failed to get header line for genotype file!');
        return;
    }
    chomp $header_line;
    my @headers = map { s/\s/_/g; s/_\-\_top$//i; lc } split(',', $header_line);
    $self->status_message('Genotype headers: '. join(', ', @headers));

    my (%genotypes, %alleles);
    my $total = 0;
    my $filters = $self->_filters;
    GENOTYPE: while ( my $line = $genotype_fh->getline ) {
        $total++;
        chomp $line;
        my %genotype;
        @genotype{@headers} = split(',', $line);
        # The id is from the snp mapping or the genotype's snp_name
        $genotype{id} = ( $snp_id_mapping and $snp_id_mapping->{ $genotype{snp_name} } )
        ? $snp_id_mapping->{ $genotype{snp_name} }
        : $genotype{snp_name};
        # Filter...
        for my $filter ( @$filters ) {
            next GENOTYPE if not $filter->filter(\%genotype);
            #print $filter->class." => $line\n" and next GENOTYPE if not $filter->filter(\%genotype);
        }
        if ( exists $genotypes{ $genotype{id} } ) {
            $self->error_message('Already have a genotype for snp id: '.Dumper(\%genotype, $genotypes{ $genotype{id} }));
            return;
        }
        $genotypes{ $genotype{id} } = \%genotype;
        $genotype{alleles} = $genotype{allele1}.$genotype{allele2};
        $alleles{ $genotype{allele1}.$genotype{allele2} }++;
    }
    $self->alleles(\%alleles);

    if ( not %genotypes ) {
        $self->error_message("None of the $total genotypes survived filtering!");
        return;
    }

    my $kept = keys(%genotypes);
    $self->genotypes_loaded($total);
    $self->genotypes_kept($kept);
    my $filtered = $total - $kept;
    $self->genotypes_filtered($filtered);
    $self->status_message("Kept $kept of $total genotypes. Filtered $filtered. Load...OK");

    return \%genotypes;
}

sub _annotate_genotypes {
    my ($self, $genotypes) = @_;

    Carp::confess('No genotypes!') if not $genotypes;

    my $variation_list_build = $self->variation_list_build;
    $self->status_message('Variant list name: '.$variation_list_build->model_name);
    $self->status_message('Variant list version: '.$variation_list_build->version);

    my $snvs_file = $variation_list_build->snvs_bed;
    if ( not $snvs_file ) {
        $self->error_message('No snvs file (snvs_bed) for build: '.$variation_list_build->__display_name__);
        return;
    }
    $self->status_message('snvs file: '.$snvs_file);

    my $dbsnp_fh = eval{ Genome::Sys->open_file_for_reading($snvs_file); };
    if ( not $dbsnp_fh ) {
        $self->error_message("Failed to open file: $snvs_file");
        return;
    }

    $self->status_message("Annotate genotypes...");
    my $variant_id_pos = ( $variation_list_build->version and $variation_list_build->version eq 130 ? 8 : 7 );
    my %annotated_genotypes;
    my $cnt = 0;
    while ( my $line = $dbsnp_fh->getline ) {
        #while ( my $line = $dbsnp_fh->getline and %$genotypes ) {
        chomp $line;
        my @tokens = split(/\s+/, $line);
        my $variant_id = $tokens[$variant_id_pos];
        if ( exists $annotated_genotypes{$variant_id} ) {
            if ( $annotated_genotypes{$variant_id}->{position} != $tokens[2] ) {
                $annotated_genotypes{$variant_id}->{ignore} = 1;
            }
            next;
        }
        my $genotype = delete $genotypes->{$variant_id};
        next if not $genotype;
        $genotype->{chromosome} = $tokens[0];
        $genotype->{position} = $tokens[2];
        $genotype->{order} = $cnt++;
        $annotated_genotypes{$variant_id} = $genotype;
    }

    if ( not %annotated_genotypes ) {
        $self->error_message("None of the ".keys(%$genotypes)." genotypes survived annotation!");
        return;
    }

    $self->status_message("Annotate ".keys(%annotated_genotypes)." genotypes...OK");
    return \%annotated_genotypes;
}

sub _output_genotypes {
    my ($self, $genotypes) = @_;
    Carp::confess('No genotypes!') if not $genotypes;
    $self->status_message('Output genotypes...');

    my $output_fh = $self->_open_output;
    return if not $output_fh;

    my $sep = ( $self->separator eq 'tab' ? "\t" : $self->separator );
    my @fields = $self->fields;
    $output_fh->print( join($sep, @fields)."\n" ) if $self->headers;

    my $cnt = 0;
    for my $genotype ( sort { $a->{order} <=> $b->{order} } values %$genotypes ) {
        next if $genotype->{ignore};
        $output_fh->print( join($sep, map { defined $genotype->{$_} ? $genotype->{$_} : 'NA' } @fields)."\n" );
        $cnt++;
    }
    $output_fh->flush;
    $output_fh->close;

    $self->genotypes_output($cnt);
    $self->status_message("Output $cnt of ".$self->genotypes_kept." genotypes. Ignored ".($self->genotypes_kept - $cnt).". Output...OK");
    if ( $cnt == 0 ) {
        $self->status_message("No genotypes output. Removing output file.");
        unlink $self->output if -e $self->output;
    }

    return 1;
}

1;

