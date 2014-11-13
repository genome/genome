package Genome::Model::PhenotypeCorrelation::Command::BurdenTestSummary;

# 3 annotation files:
#   things that passed the BT
#   things that got into but failed BT
#   everything else for that gene
# 1 big vcf

use Carp qw/confess/;
use Compress::Zlib;
use Data::Dumper;
use Genome;
use Genome::File::Vep::Reader;
use Genome::File::Vep::Writer;
use Genome::File::Vcf::Reader;

use Sort::Naturally qw/nsort/;
use Storable qw/dclone/;
use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::BurdenTestSummary {
    is => "Command::V2",
    has => [
        summary_file => {
            is => "Text",
            doc => "Path to burden analysis summary file",
        },
        burden_results_directory => {
            is => "Text",
            doc => "Path to the burden test results",
        },
        top_n => {
            is => "Text",
            doc => "Return the top n most significant genes for each given test",
        },
        test_names => {
            is => "Text",
            doc => "The names of the tests to produce data for (",
            is_many => 1,
        },
        trait => {
            is => "Text",
            doc => "The name of the trait to summarize",
        },
        clinical_data_file => {
            is => "Text",
            doc => "The path to the clinical data file used to partition samples into groups",
        },
        vcf_file => {
            is => "Text",
            doc => "The bgzipped, tabix indexed, vcf file to use",
        },
        per_site_report_file => {
            is => "Text",
            is_optional => 1,
            doc => "Path to bgzip-compressed per site callset metrics file (from joinx vcf-report)",
        },
        full_annotation_file => {
            is => "Text",
            doc => "Annotation file containing variant data and gene names",
        },
        burden_matrix_file => {
            is => "Text",
            doc => "The burden matrix",
        },
        output_directory => {
            is => "Text",
            doc => "The output path to write files to",
        },
        annotation_build_id => {
            is => "Text",
            doc => "The id of the ImportedAnnotation build to use",
            default_value => $ENV{GENOME_DB_ENSEMBL_DEFAULT_IMPORTED_ANNOTATION_BUILD},
        },
    ],
    has_transient_optional => [
        _trait_values => {
            is => "HASH",
        },
        _burden_summary => {
            is => "HASH",
        },
        _full_annotation => {
            is => "HASH",
            doc => "Internal: hashref of gene_name => [vep entries]",
        },
        _burden_input_variant_ids => {
            is => "HASH",
            doc => "Internal: Set of variant ids from the burden matrix file",
        },
    ]
};


# join path components by / and prepend the output directory
sub _output_path {
    my $self = shift;
    return join("/", $self->output_directory, @_);
}

sub execute {
    my $self = shift;

    # Load clinical data and get trait values for each sample
    my $cdata = Genome::Model::PhenotypeCorrelation::ClinicalData->from_file($self->clinical_data_file);
    my @sample_names = $cdata->sample_names;
    my %trait_values;
    @trait_values{@sample_names} = $cdata->attribute_values_for_samples($self->trait, @sample_names);
    $self->_trait_values(\%trait_values);

    # Load burden summary and find interesting genes
    $self->_burden_summary($self->_load_burden_summary($self->summary_file));
    my @genes = $self->_top_genes_for_tests();
    $self->status_message("Reporting on genes:\n\t" . join("\n\t", @genes) . "\n");

    # Load annotation and burden test data
    $self->_full_annotation($self->_read_annotation_for_genes($self->full_annotation_file, @genes));
    $self->_burden_input_variant_ids($self->_variant_ids_from_burden_matrix);

    # If we've gotten this far, we will need somewhere to write
    Genome::Sys->create_directory($self->output_directory)
        unless -d $self->output_directory;

    # Process each gene
    $self->_process_gene($_) for (@genes);

    return 1;
}

sub _process_gene {
    my ($self, $gene) = @_;
    $self->status_message("Processing gene $gene...");

    # Create subdirectory for this gene
    my $output_dir = $self->_output_path($gene);
    Genome::Sys->create_directory($output_dir);

    # All annotation for this gene
    my $full_annotation = $self->_full_annotation->{$gene};

    # Build a list of regions we can feed to tabix to create a subset of the big vcf
    my @all_locations = nsort keys %{{map {$_->{location} => undef} @$full_annotation}};
    my @all_regions = map { my @f = split(":"); "$f[0]:$f[1]-$f[1]"; } @all_locations;
    confess "Failed to find any variants to analyze for gene $gene" unless @all_regions;

    # Write output vcf, per-site report, and burden summary for this gene
    my $output_vcf = $self->_output_path($gene, "mutations.vcf");
    my $output_psr = $self->_output_path($gene, "per_site_report.txt");
    my $output_burden = $self->_output_path($gene, "burden.csv");
    $self->_write_vcf_subset($output_vcf, @all_regions);
    $self->_write_per_site_subset($output_psr, @all_regions) if $self->per_site_report_file;
    $self->_write_burden_summary_subset($gene, $self->trait, $output_burden);

    # Read in the vcf entries we just selected
    my $vcf = Genome::File::Vcf::Reader->new($output_vcf);
    my @vcf_entries;
    while (my $entry = $vcf->next) {
        next if $entry->is_filtered;
        push(@vcf_entries, $entry)
    }

    # Update full annotation information with information from the vcf
    # this adds COUNT and AF for case/control as well as vcf info fields
    $self->_update_annotation_with_vcf($full_annotation, \@vcf_entries, $vcf->header);

    # Partition annotation into 3 lists
    my %burden_output_variant_ids = $self->_variant_ids_from_burden_result($gene);
    my ($rare, $deleterious, $common) = $self->_partition_vep_entries(
        $full_annotation, # arrayref of vep entries
        \%burden_output_variant_ids, # hashref {id => undef} of things that came through the burden test
        $self->_burden_input_variant_ids # as above, but things that went into the burden test
        # since the set of output ids comes first, it takes priority.
        );

    unless (@$rare) {
        warn "Unable to find any rare variants for gene $gene, something must be wrong.";
    }

    # Write out the 3 lists of annotation
    my $rare_file = $self->_output_path($gene, "rare-deleterious.vep");
    my $deleterious_file = $self->_output_path($gene, "deleterious.vep"); # watch out for this file ;)
    my $common_file = $self->_output_path($gene, "common.vep");
    $self->_write_annotation($rare, $rare_file);
    $self->_write_annotation($deleterious, $deleterious_file);
    $self->_write_annotation($common, $common_file);

    # Create mutation-diagram plots
    # For now, we need to look at the annotation file produced when the burden matrix is created.
    $self->_make_mutation_diagrams($rare_file, $output_dir, 'case.', 'CASE_COUNT');
    $self->_make_mutation_diagrams($rare_file, $output_dir, 'control.', 'CONTROL_COUNT');
}

# _partition_vep_entries
# params:
#   $entries - an arrayref of Genome::File::Vep::Entry
#   @id_arrays - an array of sets of variant ids (hashrefs: id => undef)
# returns:
#   n+1 disjoint lists of entries where n = length of @id_arrays.
#   the last list will catch any entries whose ids are not in any of the supplied sets.
# notes:
#   sets that appear earlier in @id_arrays will take precedence over those appearing
#   later.
sub _partition_vep_entries {
    my ($self, $entries, @id_sets) = @_;

    # Create array of n+1 arrayrefs to hold results
    my @results = map {[]} 0..scalar(@id_sets);

    # Merge the sets into one hashref mapping ids to the index of the
    # first set they appear in. Then, the index of the output list for
    # a given entry is just $ids_to_indices{$id} (assuming $id exists).
    my %ids_to_indices;
    # roll through the list of sets in reverse order so things with lower
    # indices can clobber things that come later
    for my $idx (0..$#id_sets) {
        my $ridx = $#id_sets - $idx;
        my $id_set = $id_sets[$ridx];
        my @ids = keys %{$id_sets[$ridx]};
        @ids_to_indices{@ids} = ($ridx) x @ids;
    }

    for my $entry (@$entries) {
        my $idx = -1;
        if (exists $ids_to_indices{$entry->{uploaded_variation}}) {
            $idx = $ids_to_indices{$entry->{uploaded_variation}};
        }
        push(@{$results[$idx]}, $entry);
    }
    return @results;
}

# Add INFO,AF,COUNT for case/control to vep entries
sub _add_vcf_data_to_vep {
    my ($vep, $vcf, $cases, $controls) = @_;

    my $alt = $vep->{allele};
    my $alt_idx = $vcf->allele_index($alt);

    # Compute allele count and relative frequency for each of case, control
    my ($n_case_alleles, %case_allele_counts) = $vcf->allelic_distribution(@$cases);
    my ($n_control_alleles, %control_allele_counts) = $vcf->allelic_distribution(@$controls);

    my $case_count = $case_allele_counts{$alt_idx} || 0;
    my $control_count = $control_allele_counts{$alt_idx} || 0;

    my $case_af = $n_case_alleles ? $case_count / $n_case_alleles : 0;
    my $control_af = $n_control_alleles ? $control_count / $n_control_alleles : 0;

    $vep->set_extra_field("REF", $vcf->{reference_allele});
    $vep->set_extra_field("CASE_COUNT", $case_count);
    $vep->set_extra_field("CASE_AF", $case_af);
    $vep->set_extra_field("CONTROL_COUNT", $control_count);
    $vep->set_extra_field("CONTROL_AF", $control_af);

    # Add relevant info fields from vcf to vep entries
    my $info_fields = $vcf->info_for_allele($alt);
    for my $iname (keys %$info_fields) {
        $vep->set_extra_field($iname, $info_fields->{$iname});
    }
}

sub _update_annotation_with_vcf {
    my ($self, $annotation, $vcf_entries, $vcf_header) = @_;

    # Get sample names and column offsets
    my @samples = $vcf_header->sample_names;
    my %sample_indices;
    @sample_indices{@samples} = 0..$#samples;

    # We don't want any samples in the vcf that we don't have clinical data for
    my @samples_missing_cdata = grep {!exists $self->_trait_values->{$_}} @samples;
    if (@samples_missing_cdata) {
        my $t = $self->trait;
        my $f = $self->clinical_data_file;
        confess "No clinical data for trait $t found in clinical data file $f for samples:\n\t"
            . join("\n\t", @samples_missing_cdata);
    }

    # Dichotomize
    my @cases = @sample_indices{grep {$self->_trait_values->{$_} == 1} @samples};
    my @controls = @sample_indices{grep {$self->_trait_values->{$_} == 0} @samples};

    # build a hash of annotation entries keyed on loc,alt (e.g., 1:30,C <- snp at chr1:30, alt=C)
    # now we can look them up quickly given a vcf entry
    my %anno;
    for my $a (@{$annotation}) {
        my $key = join(",", $a->{location}, $a->{allele});
        $anno{$key} = [] if !exists $anno{$key};
        push(@{$anno{$key}}, $a);
    }

    # Roll through the vcf entries updating any matching annotation we can find
    for my $site (@$vcf_entries) {
        for my $alt (@{$site->{alternate_alleles}}) {
            my $key = sprintf("%s:%d,%s", $site->{chrom}, $site->{position}, $alt);
            next unless exists $anno{$key};
            _add_vcf_data_to_vep($_, $site, \@cases, \@controls) for (@{$anno{$key}});
        }
    }
}

sub _write_annotation {
    my ($self, $entries, $path) = @_;
    my $ofh = Genome::File::Vep::Writer->new($path);
    $ofh->write($_) for (@$entries);
}

sub _read_annotation_for_genes {
    my ($self, $path, @genes) = @_;
    my %genes = map {$_ => undef} @genes;

    $self->status_message("Loading annotation data from $path...");
    my $vep = Genome::File::Vep::Reader->new($path);
    my $result = {};
    while (my $entry = $vep->next) {
        # replace gene name by HGNC field if present
        if (exists $entry->{extra}{HGNC} && $entry->{gene}) {
            $entry->{gene} = $entry->{extra}{HGNC};
        }

        $entry->{uploaded_variation} = "$entry->{gene}_$entry->{uploaded_variation}";
        $entry->{uploaded_variation} =~ s/[^A-Za-z0-9]/_/g;

        my $gene = $entry->{gene};
        if (exists $genes{$gene}) {
            $result->{$gene} = [] unless exists $result->{$gene};
            push(@{$result->{$gene}}, $entry)
        }
    }

    return $result;
}

# Use tabix to grab interesting sites from the main vcf
sub _write_vcf_subset {
    my ($self, $output_file, @regions) = @_;

    my $input_file = $self->vcf_file;

    $self->status_message("Copying " . scalar(@regions) .
        " variants from $input_file to $output_file");
    my $tmp_unsorted = Genome::Sys->create_temp_file_path;
    my $cmd = Genome::Model::Tools::Tabix::Fetch->create(
        input_file => $input_file,
        regions => \@regions,
        output_file => $tmp_unsorted,
        print_header => 1,
        );
    if (!$cmd->execute) {
        confess "Failed to fetch vcf entries from file $input_file"
    }
    $cmd = Genome::Model::Tools::Joinx::Sort->create(
        input_files => [$tmp_unsorted],
        output_file => $output_file
    );
    if (!$cmd->execute) {
        confess "Failed to sort vcf file $tmp_unsorted => $output_file";
    }
}

# Use tabix to grab interesting sites from vcf-report (per-site)
sub _write_per_site_subset {
    my ($self, $output_file, @regions) = @_;
    my $input_file = $self->per_site_report_file;

    my $tmp_unsorted = Genome::Sys->create_temp_file_path;
    my $cmd = Genome::Model::Tools::Tabix::Fetch->create(
        input_file => $input_file,
        regions => \@regions,
        output_file => $tmp_unsorted,
        print_header => 1,
        );
    if (!$cmd->execute) {
        confess "Failed to fetch per site regions from file $input_file"
    }

    my $tmp_sorted = Genome::Sys->create_temp_file_path;
    $cmd = Genome::Model::Tools::Joinx::Sort->create(
        input_files => [$tmp_unsorted],
        output_file => $tmp_sorted
    );
    if (!$cmd->execute) {
        confess "Failed to sort per-site file $tmp_unsorted => $tmp_sorted";
    }

    # add the header back in
    my $ofh = Genome::Sys->open_file_for_overwriting($output_file);
    my $ifh = gzopen($input_file, "r");
    my $header;
    $ifh->gzreadline($header);
    $ofh->print($header);
    # append the tempfile
    $ifh = Genome::Sys->open_file_for_reading($tmp_sorted);
    my $rv = {};
    while (my $line = $ifh->getline) {
        $ofh->print($line);
        chomp $line;
        my @fields = split("\t", $line);
        my $chrom_pos_alt = join(":", @fields[0,1,3]);
        $rv->{$chrom_pos_alt} = \@fields;
    }
    return $rv;
}

# Make a mutation diagram for the specified Vep annotation file
sub _make_mutation_diagrams {
    my ($self, $annotation_file, $output_dir, $prefix, $frequency_field) = @_;
    Genome::Sys->create_directory($output_dir);
    my $b = Genome::Model::Build->get($self->annotation_build_id);
    my %params = (
        annotation => $annotation_file,
        annotation_format => 'vep',
        annotation_build_id => $self->annotation_build_id,
        output_directory => $output_dir,
        file_prefix => $prefix,
        vep_frequency_field => $frequency_field,
        );

    my $cmd = Genome::Model::Tools::Graph::MutationDiagram->create(%params);
    if (!$cmd->execute) {
        confess $self->error_message("Failed to create mutation diagram with params: " . Dumper(\%params));
    }
}

# Load the burden_summary.csv file from the burden test containing results for all traits/genes.
sub _load_burden_summary {
    my ($self, $file) = @_;

    my $fh = Genome::Sys->open_file_for_reading($file);
    my $header_line = $fh->getline;
    chomp $header_line;
    my @header_fields = split(",", $header_line);
    my $bs = {
        header => \@header_fields,
        column_indices => {},
        traits => {},
    };
    @{$bs->{column_indices}}{@header_fields} = 0..$#header_fields;

    while (my $line = $fh->getline) {
        chomp $line;
        my @fields = split(",", $line);
        next unless @fields;
        my ($trait, $gene) = @fields[0,1];
        next unless $trait eq $self->trait;
        $bs->{traits}{$trait} = {} unless exists $bs->{traits}{$trait};
        confess "Found multiple entries for trait $trait, gene $gene in $file!"
            if exists $bs->{traits}{$trait}{$gene};
        $bs->{traits}{$trait}{$gene} = \@fields;
    }

    return $bs;
}

sub _write_burden_summary_subset {
    my ($self, $gene, $trait, $output_file) = @_;

    my $bs = $self->_burden_summary;
    confess "Attempted to write burden summary info for trait $trait, gene $gene which does not exist"
        unless exists $bs->{traits}{$trait}{$gene};

    my $ofh = Genome::Sys->open_file_for_writing($output_file);
    $ofh->print(join(",", @{$bs->{header}}) . "\n");
    $ofh->print(join(",", @{$bs->{traits}{$trait}{$gene}}) . "\n");
}

# Find interesting genes (lowest n p-values for specified tests)
sub _top_genes_for_tests {
    my $self = shift;

    my @test_names = $self->test_names;
    my $top_n = $self->top_n;
    $self->status_message("Finding the $top_n most significant genes for tests:\n\t"
        . join("\n\t", @test_names) . "\n");

    my $bs = $self->_burden_summary;
    my @unknown_tests = grep {!exists $bs->{column_indices}->{$_}} @test_names;
    if (@unknown_tests) {
        confess "The following test names were not found in the header:\n\t"
            .join("\n\t", @unknown_tests) . "\n"
            ."The columns present are:\n\t"
            .join("\n\t", @{$bs->{header}});
    }

    my @test_indices = @{$bs->{column_indices}}{@test_names};
    my %results;
    for my $trait (keys %{$bs->{traits}}) {
        for my $gene (keys %{$bs->{traits}{$trait}}) {
            my $fields = $bs->{traits}{$trait}{$gene};
            unless (exists $results{$trait}) {
                $results{$trait} = {};
            }
            my @values = @{$fields}[@test_indices];
            for my $idx (0..$#test_names) {
                my $test_name = $test_names[$idx];
                next if $values[$idx] eq "NA";
                $results{$trait}{$test_name} = [] unless exists $results{$trait}{$test_name};
                _insert_result($results{$trait}{$test_name}, $gene, $values[$idx], $top_n);
            }
        }
    }

    my %genes = map {$_->[0] => 1} map {@$_} map {values %{$results{$_}} } keys %results;
    return sort keys %genes;
}

# Helper to build up list of interesting genes
sub _insert_result {
    my ($arr, $gene, $score, $max_size) = @_;
    my $idx = 0;
    while ($idx <= $#$arr && $score > $arr->[$idx]->[1]) {
        ++$idx;
    }
    if ($idx < $max_size) {
        splice(@$arr, $idx, 0, [$gene, $score]);
        pop(@$arr) if $#$arr >= $max_size;
    }
}

# Look up what variants the burden test used
sub _variant_ids_from_burden_result {
    my ($self, $gene) = @_;
    my $filename = $self->trait . "_$gene.single.csv";
    my $path = join("/", $self->burden_results_directory, $filename);

    my $fh = Genome::Sys->open_file_for_reading($path);
    $fh->getline;
    my %variants;
    while (my $line = $fh->getline) {
        chomp $line;
        my ($variant, @fields) = split(",", $line);
        # remove leading V that the burden test appends to variant names
        $variant =~ s/^V//g;
        $variants{$variant} = 1;
    }
    return %variants;
}

sub _variant_ids_from_burden_matrix {
    my $self = shift;
    my $fh = Genome::Sys->open_file_for_reading($self->burden_matrix_file);
    # discard header
    $fh->getline;
    my $result = {};
    while (my $line = $fh->getline) {
        chomp $line;
        next unless $line;
        my ($id, @rest) = split("\t", $line);
        $result->{$id} = 1;
    }
    return $result;
}

1;
