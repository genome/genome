package Genome::Model::PhenotypeCorrelation::Command::BurdenTestSummary;

use Carp qw/confess/;
use Compress::Zlib;
use Data::Dumper;
use Genome;
use Genome::File::Vep::Reader;
use Genome::File::Vep::Writer;

use Sort::Naturally qw/nsort/;
use Storable qw/dclone/;
use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::BurdenTestSummary {
    is => "Genome::Command::Base",
    has => [
        maximum_maf => {
            is => "Number",
            doc => "Maximum minor allele frequency cutoff to include",
            default_value => 0.01,
        },
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
            doc => "Return the top n most significant genes for each trait",
        },
        trait => {
            is => "Text",
            doc => "The name of the trait to summarize",
        },
        test_names => {
            is => "Text",
            doc => "The names of the tests to produce data for (",
            is_many => 1,
        },
        clinical_data_file => {
            is => "Text",
            doc => "The path to the clinical data file used to partition samples into groups (e.g., case/control)",
        },
        vcf_file => {
            is => "Text",
            doc => "The bgzipped vcf file to use",
        },
        per_site_report_file => {
            is => "Text",
            is_optional => 1,
            doc => "Path to bgzip-compressed per site callset metrics file",
        },
        annotation_file => {
            is => "Text",
            doc => "Annotation file containing variant data and gene names",
        },
        output_directory => {
            is => "Text",
            doc => "The output path to write files to",
        },
        annotation_build_id => {
            is => "Text",
            doc => "The id of the ImportedAnnotation build to use",
            default_value => $ENV{GENOME_DB_ENSEMBL_DEFAULT_IMPORTED_ANNOTATION_BUILD},
        }
    ],
    has_transient_optional => [
        _trait_values => {
            is => "HASH",
        },
        _burden_summary => {
            is => "HASH",
        }
    ]
};

sub execute {
    my $self = shift;

    # Load clinical data and get attribute values for each sample
    my $cdata = Genome::Model::PhenotypeCorrelation::ClinicalData->from_file($self->clinical_data_file);
    my @sample_names = $cdata->sample_names;
    my %trait_values;
    @trait_values{@sample_names} = $cdata->attribute_values_for_samples($self->trait, @sample_names);
    $self->_trait_values(\%trait_values);

    # Load burden summary and find interesting genes
    $self->_burden_summary($self->_load_burden_summary($self->summary_file));
    my @genes = $self->_top_genes_for_tests();
    $self->status_message("Reporting on genes:\n\t" . join("\n\t", @genes) . "\n");

    # If we've gotten this far we had better have somewhere to write
    Genome::Sys->create_directory($self->output_directory)
        unless -d $self->output_directory;

    # Process each gene
    $self->_process_gene($_) for (@genes);

    return 1;
}

# join path components by / and prepend the output directory
sub _output_path {
    my $self = shift;
    return join("/", $self->output_directory, @_);
}

sub _process_gene {
    my ($self, $gene) = @_;
    $self->status_message("Processing gene $gene...");

    my $output_dir = $self->_output_path($gene);
    Genome::Sys->create_directory($output_dir);

    # Load vep annotation
    #  returns { header => (header obj), data => { gene => [ (entry objs) ] } }
    my %snv_ids = $self->_snv_ids_for_gene($gene);
    my $gene_annotation = $self->_read_annotation_for_genes(genes => [$gene], include_ids => \%snv_ids);
    

    # Build a list of locations (chr:pos strings) we are interested in
    my %locations = map {$_->{location} => 1} @{$gene_annotation->{data}{$gene}};
    my @locations = nsort keys(%locations);
    confess "Unable to find any variants for gene $gene, something must be wrong." unless @locations;

    # Transform locations to ranges: 1:20 => 1:20-20
    my @regions = map {
        my @f = split(":");
        "$f[0]:$f[1]-$f[1]";
        } @locations;

    # Write out subsets of our input files
    my $output_vcf = $self->_output_path($gene, "mutations.vcf");
    my $output_psr = $self->_output_path($gene, "per_site_report.txt");
    my $output_burden = $self->_output_path($gene, "burden.csv");
    $self->_write_vcf_subset($output_vcf, @regions);
    $self->_write_per_site_subset($output_psr, @regions) if $self->per_site_report_file;
    $self->_write_burden_summary_subset($gene, $self->trait, $output_burden);

    # Read in the subset of the vcf we just created
    my $vcf = Genome::File::Vcf::Reader->create(name => $output_vcf);
    $vcf->open($output_vcf);
    my @entries;
    while (my $entry = $vcf->next) {
        push(@entries, $entry);
    }

    # Get sample names and column offsets
    my @samples = $vcf->header->sample_names;
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

    # Get allele frequencies and info fields for each of case, control.
    # TODO: refactor the next two blocks so they call one thing instead of repeating code
    my @case_data = map {
        my ($total, %counts) = $_->allelic_distribution(@cases);
        my $entry = $_;
        {
            total => $total,
            allele_counts => \%counts,
            chrompos => join(":", $_->chrom, $_->position),
            # TODO: make this next block less obnoxious later
            info_by_alt => {
                map {
                    $_ => ($entry->info_for_allele($_) || undef)
                } keys %counts
            },
        }
    } @entries;

    my @control_data = map {
        my ($total, %counts) = $_->allelic_distribution(@controls);
        my $entry = $_;
        {
            total => $total,
            allele_counts => \%counts,
            chrompos => join(":", $_->chrom, $_->position),
            # TODO: should probably make this next block less obnoxious later too
            info_by_alt => {
                map {
                    $_ => ($entry->info_for_allele($_) || undef)
                } keys %counts
            }
        }
    } @entries;

    # Write out new annotation files for each of case, control
    my $cases_annotation_file = $self->_output_path($gene, "cases.vep");
    my $controls_annotation_file = $self->_output_path($gene, "controls.vep");
    $self->_write_merged_annotation($gene_annotation, $gene, \@case_data, $cases_annotation_file);
    $self->_write_merged_annotation($gene_annotation, $gene, \@control_data, $controls_annotation_file);

    # Create mutation-diagram plots
    my $case_plots_dir = join("/", $output_dir, "case_plots");
    my $control_plots_dir = join("/", $output_dir, "control_plots");
    $self->_make_mutation_diagrams($cases_annotation_file, $output_dir, 'case.');
    $self->_make_mutation_diagrams($controls_annotation_file, $output_dir, 'control.');
}

sub _read_annotation_for_genes {
    my ($self, %params) = @_;
    # Params should be:
    # genes => [gene names],
    # include_ids => { valid_snv_id => 1, ... }
    # if include_ids is set, then annotation entries are limited to
    # those whose id (uploaded_variation) is found in the hashref
    
    my %genes = map {$_ => 1} @{$params{genes}};

    my $vep = Genome::File::Vep::Reader->new($self->annotation_file);
    
    my $rv = {
        header => $vep->{header},
        data => {}
    };

    # Read all the entries where entry->gene in %genes
    while (my $entry = $vep->next) {
        # if we are filtering ids, skip things that don't match our list
        if (exists $params{include_ids} && !exists $params{include_ids}->{$entry->{uploaded_variation}}) {
            next;
        }

        my $gene = $entry->{gene};
        if (exists $genes{$gene}) {
            $rv->{data}{$gene} = [] unless exists $rv->{data}{$gene};
            push(@{$rv->{data}{$gene}}, $entry);
        }
    }
    return $rv;
}

# Given some annotation data, and some extra data about the snps (vcf info fields, allele
# frequency, ...), extend the Vep annotation by appending this snp info to the EXTRA
# field and write the result to disk.
sub _write_merged_annotation {
    my ($self, $annotation, $gene, $snp_data, $output_file) = @_;
    my $vep = new Genome::File::Vep::Writer($output_file);

    # build a hash of annotation entries keyed on loc,alt (e.g., 1:30,C <- snp at chr1:30, alt=C)
    # now we can look them up quickly given a snp
    my %anno = map { join(",", $_->{location}, $_->{allele}) => $_ } @{$annotation->{data}{$gene}};

    # for each site
    for my $si (@{$snp_data}) {
        my @alleles = keys %{$si->{allele_counts}};
        # for each snp at the site
        for my $alt (@alleles) {
            # look up the annotation for this snp or skip it
            my $key = join(",", $si->{chrompos}, $alt);
            next unless exists $anno{$key};

            # we want to modify a copy of the entry
            my $entry = $anno{$key}->clone();

            my $count = $si->{allele_counts}{$alt};
            my $af = $count / $si->{total};

            # add our own extra annotation fields
            $entry->set_extra_field("AF", $af);
            $entry->set_extra_field("COUNT", $count);
            # copy info fields from vcf
            for my $info_name (keys %{$si->{info_by_alt}{$alt}}) {
                $entry->set_extra_field($info_name, $si->{info_by_alt}{$alt}{$info_name});
            }
            $vep->write($entry);
        }
    }
}

# Use tabix to grab interesting sites from the main vcf
sub _write_vcf_subset {
    my ($self, $output_file, @regions) = @_;

    my $input_file = $self->vcf_file;

    $self->status_message("Copying " . scalar(@regions) .
        " snvs from $input_file to $output_file");
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
    my ($self, $annotation_file, $output_dir, $prefix) = @_; 
    Genome::Sys->create_directory($output_dir);
    my $b = Genome::Model::Build->get($self->annotation_build_id);
    my %params = (
        annotation => $annotation_file,
        annotation_format => 'vep',
        annotation_build_id => $self->annotation_build_id,
        output_directory => $output_dir,
        file_prefix => $prefix
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
    print "Finding the $top_n most significant genes for tests:\n"
        . join("\n", @test_names) . "\n";

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

    my @genes = map {$_->[0]} map {@$_} map {values %{$results{$_}} } keys %results;
    return @genes;
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

# Look up what snvs the burden test used
sub _snv_ids_for_gene {
    my ($self, $gene) = @_;
    my $filename = $self->trait . "_$gene.single.csv";
    my $path = join("/", $self->burden_results_directory, $filename);

    my $fh = Genome::Sys->open_file_for_reading($path);
    $fh->getline;
    my %snvs;
    while (my $line = $fh->getline) {
        chomp $line;
        my ($snv, @fields) = split(",", $line);
        # remove leading V that the burden test appends to variant names
        $snv =~ s/^V//g;
        $snvs{$snv} = 1;
    }
    return %snvs;
}

1;
