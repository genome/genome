package Genome::Model::Tools::BamQc::SummarizeInstrumentData;

use strict;
use warnings;

use File::Slurp qw(read_file);
use Carp qw(confess);
use Genome;

my $R_SCRIPT = sprintf("%s.R", __FILE__);

my %SAMPLE_PROPERTIES = (name => "sample_name", common_name => "common_name");
my %LIBRARY_PROPERTIES = (name => "library_name");
my %INSTRUMENT_DATA_PROPERTIES = (
      id => "id"
    , flow_cell_id => "flow_cell_id"
    , lane => "lane"
    , read_count => "read_count"
    );

class Genome::Model::Tools::BamQc::SummarizeInstrumentData {
    is => "Command::V2",
    has_input => [
        instrument_data => {
            is => "Genome::InstrumentData",
            is_many => 1,
            shell_args_position => 1,
            doc => "Instrument data to summarize"
        },

        output_tsv => {
            is => "String",
            doc => "Summary output file"
        },

        output_pdf => {
            is => "String",
            doc => "Summary pdf file",
            is_optional => 1
        }
    ],
};

sub _make_bamqc {
    my $obj = shift;

    my $path = $obj->disk_allocations->absolute_path;

    my $x = {
        obj => $obj,
        path => $path,
        get_path => sub {
            my $expr = shift;
            my $test = File::Spec->catfile($path, $expr);
            my @results = glob($test);
            return $results[0] if @results;
        },
    };

    return $x;
}

sub _get_duplication_rate {
    my $bam_qc = shift;
    my $path = $bam_qc->{get_path}->("*_fastqc/fastqc_data.txt");
    confess sprintf("Failed to find fastqc report for bamqc %s", $bam_qc->{obj}->id) unless $path;
    my @lines = grep {/^#Total Duplicate Percentage/} read_file($path);
    chomp @lines;
    confess "Failed to find duplicate percentage in fastqc file $path" unless @lines == 1;
    my ($x, $percent) = split("\t", $lines[0]);
    return sprintf("%.03f", $percent);
}

sub _find_bamqc {
    my $idata = shift;

    my @alignment_results = $idata->get_default_alignment_results;
    confess sprintf("No alignment results for instrument data %s", $idata->id) unless @alignment_results;
    my @ids = map {$_->id} @alignment_results;
    my @bam_qc = Genome::InstrumentData::AlignmentResult::Merged::BamQc->get(
        'alignment_result_id in' => \@ids);

    confess sprintf("No BamQc found for instrument data %s", $idata->id) unless @bam_qc;

    my @available = grep {!$_->disk_allocations->is_archived} @bam_qc;
    if (!@available) {
        confess sprintf("All BamQc results for instrument data %s are archived:\n\t%s",
            join("\n\t", map {
                    sprintf("BamQc id: %s, Allocation id: %s", $_->id, $_->disk_allocations->id)
                } @bam_qc)
            );
    }

    return _make_bamqc($bam_qc[0]);
}

sub _fields {
    return
          values %INSTRUMENT_DATA_PROPERTIES
        , values %LIBRARY_PROPERTIES
        , values %SAMPLE_PROPERTIES
        , "duplicate_percent"
        , "bamqc_path"
        ;
}

sub _generate_plots {
    my $self = shift;
    return 1 unless $self->output_pdf;
    my @cmd = ("Rscript", $R_SCRIPT, $self->output_tsv, $self->output_pdf);
    return Genome::Sys->shellcmd(
          cmd => join(" ", @cmd)
        , input_files => [$self->output_tsv]
        , output_files => [$self->output_pdf]
        );
}

sub execute {
    my $self = shift;

    my @idata = $self->instrument_data;
    my @bam_qcs = map {_find_bamqc($_)} @idata;

    $self->status_message(sprintf("Found %d bam qc results", scalar(@bam_qcs)));

    my @fields = _fields();

    my $tsv_fh = Genome::Sys->open_file_for_writing($self->output_tsv);
    $tsv_fh->print(join("\t", @fields), "\n");
    for my $i (0..$#idata) {
        my $id = $idata[$i];
        my $bqc = $bam_qcs[$i];
        my $sample = $id->sample;
        my $library = $id->library;

        my $dup_rate = _get_duplication_rate($bqc);
        my %properties = (
              (map {$LIBRARY_PROPERTIES{$_} => ($library->$_ // "NA")} keys %LIBRARY_PROPERTIES)
            , (map {$SAMPLE_PROPERTIES{$_} => ($sample->$_ // "NA")} keys %SAMPLE_PROPERTIES)
            , (map {$INSTRUMENT_DATA_PROPERTIES{$_} => ($id->$_ // "NA")} keys %INSTRUMENT_DATA_PROPERTIES)
            , duplicate_percent => $dup_rate
            , bamqc_path => $bqc->{path}
            );

        $tsv_fh->print(join("\t", @properties{@fields}), "\n");
    }
    $tsv_fh->close();

    return $self->_generate_plots();
}

1;
