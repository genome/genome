package Genome::InstrumentData::VerifyBamIdResult;

use strict;
use warnings;
use Genome;
use Sys::Hostname;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Writer;
use Genome::File::Vcf::DbsnpAFParser;
use Memoize qw();

use constant AF_HEADER => '<ID=AF,Number=A,Type=Float,Description="Allele frequence for non-reference alleles">';

class Genome::InstrumentData::VerifyBamIdResult {
    is => ['Genome::SoftwareResult::Stageable', 'Genome::SoftwareResult::WithNestedResults'],
    has_input => [
        aligned_bam_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
        on_target_list => {
            is => "Genome::FeatureList",
            is_optional => 1,
        },
        sample => {
            is => "Genome::Sample",
        },
        known_sites_build => {
            is => "Genome::Model::Build::ImportedVariationList",
        },
    ],
    has_param => [
        genotype_filters => {
            is => 'Text',
            is_many => 1,
        },
        max_depth => {
            is => "Integer",
        },
        precise => {
            is => 'Boolean',
        },
        version => {
            is => "Text",
        },
        result_version => {
            is => "Integer",
        },
    ],
    has_metric => [
        freemix => {
            is => "UR::Value::Number",
        },
        chipmix => {
            is => "UR::Value::Number",
        },
        af_count => {
            is => "UR::Value::Number",
        },
    ],
};

sub _error {
    my ($self, $msg) = @_;
    $self->error_message($msg);
    die $self->error_message;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return if not $self;
    $self->_error("Failed to prepare staging directory") unless $self->_prepare_staging_directory;

    $self->_error("Failed to run verifyBamID") unless $self->_run_verify_bam_id;


    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    $self->_error("Failed to add metrics") unless $self->_add_metrics;

    return $self;
}

sub _add_metrics {
    my $self = shift;
    my $self_sm = File::Spec->join($self->output_dir, "output.selfSM");
    my $in = Genome::Utility::IO::SeparatedValueReader->create(
        separator => "\t",
        input => $self_sm,
    );
    my $metrics = $in->next;
    $self->freemix($metrics->{FREEMIX});
    $self->chipmix($metrics->{CHIPMIX});
    return 1;
}

sub _run_verify_bam_id {
    my $self = shift;

    my $bam_file = $self->_resolve_bam_file;
    my $vcf_file = $self->_resolve_vcf_file;
    my $out_prefix = File::Spec->join($self->temp_staging_directory, "output");
    my $cmd = Genome::Model::Tools::VerifyBamId->create(
        vcf => $vcf_file,
        bam => $bam_file,
        out_prefix => $out_prefix,
        max_depth => $self->max_depth,
        precise => $self->precise,
        version => $self->version
    );
    return $cmd->execute if $cmd;
    return;
}

sub _resolve_bam_file {
    my $self = shift;
    my $bam_result = $self->aligned_bam_result;

    my $path = $bam_result->bam_path;
    unless (-s $path) {
        $self->_error("Could not get bam file for ".$bam_result->id);
    }
    return $path;
}

sub _resolve_vcf_file {
    my $self = shift;
    my $genotype_vcf_result = $self->_resolve_genotype_vcf_result;
    my $vcf = $genotype_vcf_result->vcf_path;
    unless (-s $vcf) {
        $self->_error("Could not get vcf file for genotype vcf".$genotype_vcf_result);
    }
    return $self->_clean_vcf($vcf);
}

sub _resolve_genotype_vcf_result {
    my $self = shift;

    my $users = $self->_user_data_for_nested_results;
    $users->{uses} = $self;

    my %params = (
        sample => $self->sample,
        known_sites_build => $self->known_sites_build,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME},
        users => $users,
    );
    if ($self->genotype_filters) {
        $params{filters} = [$self->genotype_filters];
    }
    my $result = Genome::InstrumentData::Microarray::Result::Vcf->get_or_create(%params);
    $self->_error("Could not get or create genotype vcf result") unless $result;
    return $result;
}

sub _clean_vcf {
    my $self = shift;
    my $vcf_path = shift;

    if ($self->on_target_list) {
        $self->debug_message("Using on_target_list");
        my $on_target_path = Genome::Sys->create_temp_file_path;
        my $on_target_bed = $self->on_target_list->processed_bed_file;
        my $rv = Genome::Model::Tools::BedTools::Intersect->execute(
            input_file_a => $vcf_path,
            input_file_b => $on_target_bed,
            input_file_a_format => "bed",
            intersection_type => "unique",
            output_file => $on_target_path,
            header => 1,
        );
        $self->_error("Could not intersect with on target bed") unless $rv && $rv->result;
        $vcf_path = $on_target_path;
    }

    my $fixed_frequency_path = $self->_fix_allele_frequencies($vcf_path);

    $self->debug_message("Using cleaned vcf at $fixed_frequency_path");
    return $fixed_frequency_path;
}

sub _fix_allele_frequencies {
    my $self = shift;
    my $vcf_path = shift;

    my $af_count = 0;
    my $new_vcf_path = Genome::Sys->create_temp_file_path;
    my $reader = Genome::File::Vcf::Reader->new($vcf_path);
    my $header = $reader->header;
    $header->add_info_str(AF_HEADER);
    my $writer = Genome::File::Vcf::Writer->new($new_vcf_path, $header);

    while (my $entry = $reader->next) {
        if (defined $entry->info->{AF} or (defined $entry->info->{AC} and defined $entry->info->{AN})) {
            $af_count++;
        }
        elsif (defined $entry->info->{CAF}) {
            my $af = _convert_caf_to_af($entry);
            if ($af) {
                $entry->info->{AF} = $af;
                $af_count++;
            }
        }
        $writer->write($entry);
    }
    $self->af_count($af_count);
    $writer->close;
    return $new_vcf_path;
}

sub _caf_parser {
    my $header = shift;
    return Genome::File::Vcf::DbsnpAFParser->new($header);
}

Memoize::memoize('_caf_parser');

sub _convert_caf_to_af {
    my $entry = shift;
    my @fields;
    my $parser = _caf_parser($entry->{header});
    my $caf = eval {$parser->process_entry($entry);};
    my $error = $@;
    #Allow this error until we fix some other problems with CAF
    if ($error and !($error =~ /Frequency list and allele list differ in length/)) {
        die $error;
    }
    unless ($caf) {
        return undef;
    }
    for my $alt (@{$entry->{alternate_alleles}}) {
        push @fields, $caf->{$alt};
    }
    return join(",", @fields);
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;

    my $user = $ENV{'USER'};
    my $base_dir = sprintf("verifybamidresult-%s-%s-%s-%s",           $hostname, $user, $$, $self->id);
    my $directory = join('/', 'build_merged_alignments',$self->id,$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}
1;

