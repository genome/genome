package Genome::Model::Tools::SmrtAnalysis::Consensus;

use strict;
use warnings;

use Genome;

use Workflow;
use Workflow::Simple;

use Statistics::Descriptive;

# thresholds which define different reference or genome scopes
# reference scopes are defined vs. total length of the reference in kbp
my %REFERENCE_SCOPES = (
    small => 10,
    large => 1_000,
    huge  => 10_000_000,
);

# thresholds which define different read scopes
# In the format scopeName:upperLimit, ...
# where upper limit is expressed in megabases of post-filtered
# sequence
# These scopes are used to classify the scope of the requested
# analysis---tested in order, first one wins

my %READ_SCOPES = (
    small => 3.6,
    large => 100,
    huge => 1_000_000,
);

my %SCOPE_TO_CHUNK = (
    'small' => -1,
    # 'large' => 10_000,
    'large' => 100_000,
    # 'huge' => 100_000,
    'huge' => 1_000_000,
);

class Genome::Model::Tools::SmrtAnalysis::Consensus {
    is =>  ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        job_directory => {
            is => 'Text',
            doc => 'The base job directory.',
        },
        cmp_hdf5_file => {
            is => 'Text',
            doc => 'An aligned reads cmp.h5 format file.',
        },
        alignment_summary_gff => {
            is => 'Text',
            doc => 'An alignment summary GFF3 format file.',
        },
    ],
    has_optional_input => [
        min_variant_quality => {
            is => 'Number',
            default_value => 0,
        },
        min_coverage => {
            is => 'Number',
            default_value => 4,
        },
        max_coverage => {
            is => 'Number',
            default_value => 500,
        },
        nproc => {
            is => 'Number',
            default_value => 1,
        },
        scope => {
            is => 'Text',
            valid_values => ['small','large','huge'],
        },
        base_map => {
            is => 'Text',
            doc => 'Comma-separted basemap string: dye numbers for (A,C,G,T)',
            default_value => '3,4,2,1',
        },
        decode_file => {
            is => 'Text',
            doc => 'File containing Decode matrix',
            is_optional => 1,
        },
        consensus_quality_path => {
            is => 'Text',
            is_optional => 1,
        },
    ],
    has_optional => {
        _contigs => {},
        _input_bp => {},
        _ref_len => {},
        _weighted_coverage => {},
    },
    has_optional_output => {
        variants_gff_file => { },
    },
    has_optional_param => [
        lsf_queue => {
            default_value => 'workflow',
        },
        lsf_resource => {
            default_value => '',
        },
    ],
};

sub help_brief {
    ''
}


sub help_detail {
    return <<EOS 

EOS
}

sub execute {
    my $self = shift;

    my $job_directory = $self->job_directory;
    my $output_directory = $job_directory .'/data';
    unless (-d $output_directory) {
        Genome::Sys->create_directory($output_directory);
    }

    unless ($self->decode_file) {
        my $decode_file = $self->seymour_home .'/analysis/etc/defaultDecode.params';
        $self->status_message('Using default decode file: '. $decode_file);
        $self->decode_file($decode_file);
    }
    unless ($self->consensus_quality_path) {
        my $consensus_quality_path = $self->seymour_home .'/analysis/etc/consensusQual.tsv';
        $self->status_message('Using default consensus quality file: '. $consensus_quality_path);
        $self->consensus_quality_path($consensus_quality_path);
    }

    my $contigs = $self->contigs;
    $self->resolve_scope;
    my @evi_cons_params;
    # TODO: Determine minimum confidence from min_variant_quality
    for my $contig (@{$contigs}) {
        my $evi_cons_params = $self->resolve_evi_cons_params_for_contig($contig);
        unless ($evi_cons_params) {
            $self->status_message('There is no need to run consensus on '. $contig .' with '. $self->contig_mean_coverage($contig) .' mean coverage! Skipping consensus!');
            next;
        }
        my $chunk_size = $self->resolve_chunk_size_for_contig($contig);
        my $contig_length = $self->contig_length($contig);
      CHUNK: for (my $i = 0; $i < $contig_length; $i += $chunk_size) {
            my $start = $i;
            my $end = ($i + $chunk_size) - 1;
            # Prevent the chunk from extending beyond the reference contig length(0-based)
            if ($end >= $contig_length) {
                $end = ($contig_length - 1);
            }
            my $chunk_params = $evi_cons_params .' --hdf5Reference='. $contig .' --refStart='. $start .' --refEnd='. $end;
            # This will be the only chunk representing the full reference contig
            if ($start == 0 and $end == ($contig_length - 1)) {
                push @evi_cons_params, $chunk_params;
                last CHUNK;
            }
            $chunk_params .= ' --subAlignment';
            push @evi_cons_params, $chunk_params;
        }
    }
    my $variants_gff_file = $output_directory .'/variants.gff';
    my %params = (
        alignment_summary_gff_file => $self->alignment_summary_gff,
        data_directory => $output_directory,
        consensus_params => \@evi_cons_params,
        cmp_hdf5_file => $self->cmp_hdf5_file,
        variants_gff_file => $variants_gff_file,
        gff_to_bed_purpose => 'variants',
        gff_to_bed_name => 'variants',
        gff_to_bed_description => 'PacBio: snps, insertions, and deletions derived from consensus calls against reference',
    );
    my $module_path = $self->get_class_object->module_path;
    my $xml_path = $module_path;
    $xml_path =~ s/\.pm/\.xml/;
    my $workflow = Workflow::Operation->create_from_xml($xml_path);
    my @errors = $workflow->validate;
    unless ($workflow->is_valid) {
        die('Errors encountered while validating workflow '. $xml_path ."\n". join("\n", @errors));
    }
    my $output = Workflow::Simple::run_workflow_lsf($xml_path,%params);
    unless (defined $output) {
        my @errors = @Workflow::Simple::ERROR;
        for (@errors) {
            print STDERR $_->error ."\n";
        }
        return;
    }
    $self->variants_gff_file($output->{variants_gff_file});
    my $directories = $output->{temp_directories};
    for my $directory (@{$directories}) {
        unless (File::Path::rmtree($directory)) {
            die('Failed to remove intermediate result_directories!  Specifically this one: '. $directory);
        }
    }
    return 1;
}

sub resolve_scope {
    my $self = shift;

    if (defined($self->scope)) { return 1; }

    my $input_mbp = ($self->input_bp / 1_000_000);
    my $ref_len_kbp = ($self->ref_len / 1_000);
    for my $scope ( sort {$REFERENCE_SCOPES{$a} <=> $REFERENCE_SCOPES{$b}} keys %REFERENCE_SCOPES) {
        $self->scope($scope);
        if ($ref_len_kbp < $REFERENCE_SCOPES{$scope}) {
            last;
        }
    }
    return 1;
}

sub contigs {
    my $self = shift;
    unless ($self->_contigs) {
        $self->_load_gff;
    }
    return $self->_contigs;
}

sub input_bp {
    my $self = shift;
    unless ($self->_input_bp) {
        $self->_load_gff;
    }
    return $self->_input_bp;
}

sub ref_len {
    my $self = shift;
    unless ($self->_ref_len) {
        $self->_load_gff;
    }
    return $self->_ref_len;
}

sub weighted_coverage {
    my $self = shift;
    unless ($self->_weighted_coverage) {
        $self->_load_gff;
    }
    return $self->_weighted_coverage;
}

sub _load_gff {
    my $self = shift;
    my $reader = Genome::Utility::IO::GffReader->create(
        input => $self->alignment_summary_gff,
    );
    my %weighted_coverage;
    my $input_bp;
    my $ref_length;
    while (my $data = $reader->next) {
        my $chr = $data->{chr};
        unless ( $weighted_coverage{$chr} ) {
            $weighted_coverage{$chr} = Statistics::Descriptive::Sparse->new();
        }
        my $attributes = $data->{attributes};
        my @attributes = split(';',$attributes);
        my %attributes;
        for my $attribute (@attributes) {
            unless ($attribute =~ /\s*(\S+)=(\S+)/) {
                die(Data::Dumper::Dumper($data));
            }
            $attributes{$1} = $2;
        }
        my $cov_2_attribute = $attributes{'cov2'};
        # TODO: verify that the second value is the stdev or variance...
        my ($mean, $stdev) = split(',',$cov_2_attribute);
        for ($data->{start} .. $data->{end}) {
            $weighted_coverage{$chr}->add_data($mean);
        }
    }
    $self->_weighted_coverage(\%weighted_coverage);
    my @contigs = keys %weighted_coverage;
    $self->_contigs(\@contigs);
    for my $contig (@contigs) {
        my $stats = $weighted_coverage{$contig};
        my $ref_len = $stats->count;
        $self->_ref_len($ref_len);
        my $input_bp = int($stats->mean * $ref_len);
        $self->_input_bp($input_bp);
    }
    return 1;
}

sub contig_length {
    my $self = shift;
    my $contig = shift;
    my $weighted_coverage = $self->weighted_coverage_for_contig($contig);
    return $weighted_coverage->count;
}

sub contig_mean_coverage {
    my $self = shift;
    my $contig = shift;
    my $weighted_coverage = $self->weighted_coverage_for_contig($contig);
    return $weighted_coverage->mean;
}

sub weighted_coverage_for_contig {
    my $self = shift;
    my $contig = shift;
    my $weighted_coverage = $self->weighted_coverage;
    if (defined($weighted_coverage->{$contig})) {
        return $weighted_coverage->{$contig};
    }
    return;
}

sub resolve_evi_cons_params_for_contig {
    my $self = shift;
    my $contig = shift;

    my $evi_cons_params;
    my $mean_coverage = $self->contig_mean_coverage($contig);

    my $plurality_opts = '--baseMap='. $self->base_map .' --threshold=0.0';
    my $evicons_opts = '--fastMode --baseMap='. $self->base_map .' --runDecode --decodeFile='. $self->decode_file .' --nproc='. $self->nproc;
    if ($mean_coverage == 0) {
        return;
    } elsif ($mean_coverage > $self->max_coverage) {
        $evi_cons_params = $plurality_opts;
    } elsif ($mean_coverage < $self->min_coverage) {
        $evi_cons_params = $plurality_opts;
    } elsif ($mean_coverage >= $self->min_coverage) {
        $evi_cons_params = $evicons_opts;
    }
    return $evi_cons_params;
}

sub resolve_chunk_size_for_contig {
    my $self = shift;
    my $contig = shift;
    my $chunk_size = $SCOPE_TO_CHUNK{$self->scope};
    if ($chunk_size == -1) {
        return $self->contig_length($contig);
    }
    return $chunk_size;
}

1;
