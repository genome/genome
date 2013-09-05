package Genome::Model::Tools::Convergence::GatherRefalignSummaryData;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::Convergence::GatherRefalignSummaryData{
    is => 'Command',
    has => [
        build_id => {
            is => 'Text',
            is_input => 1,
            doc => 'the build for which to gather the summary report data',
        },
        _build => {
            is => 'Genome::Model::Build',
            id_by => 'build_id',
        }
    ],
    has_optional => [
        data => {
            is_output => 1,
            doc => 'after this tool has been execute()d, this holds the extracted data (any input is ignored)',
        }
    ],
};

sub help_brief {
    "extract data from a reference alignment",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt convergence gather-refalign-summary-data --build-id 101091022
EOS
}

sub help_detail {
    return <<EOS 
Reference Alignments have produced data that is only present in the reports.  This tool's
main purpose is to get it back into a convenient machine-usable form.  If this tool is
called from within another script, the data parameter holds a hash of the data extracted.

This tool attempts to be somewhat lenient about the contents of the summary report.txt,
but assumes that the data it seeks are in a certain order and have an exactly matching
label. (e.g. 'unfiltered SNP calls:' for the unfiltered SNP calls.)
EOS
}

sub execute {
    my $self = shift;

    my $build = $self->_build;
    unless($build) {
        $self->error_message('Could not get build for ' . $self->build_id);
        return;   
    }

    my $reports_directory = $build->resolve_reports_directory;
    unless($reports_directory) {
        $self->error_message('Could not get report directory');
        return;
    }

    #XML doesn't contain data due to deprecated report format--try text version
    my $summary_report = $reports_directory . "Summary/report.txt";
    $self->status_message('Using ' . $summary_report);

    unless(-e $summary_report and -f $summary_report) {
        $self->error_message('Summary report file ' . $summary_report . ' does not exist or is not a file.');
        return;   
    }
    
    unless($self->parse_file($summary_report) and $self->data) {
        $self->error_message('No data extracted.');
        return;
    }

    $self->output_values;
    
    return 1;
}

sub output_values {
    my $self = shift;

    print Data::Dumper::Dumper $self->data;
    
    return 1;
}

sub parse_file {
    my $self = shift;
    my $file = shift;

    my $fh = Genome::Sys->open_file_for_reading($file);
    
    unless($fh) {
        $self->error_message('Could not open file: ' . Genome::Sys->error_message);
        return;
    }

    my $phrases = $self->phrases;

    my %values;
    
    my ($key, $regex) = $self->next_phrase($phrases);

    while(<$fh>) {
        my $line = $self->trim($_);
        
        if(my ($value) = $line =~ m/^$regex$/) {
            $values{$key} = $value;

            last unless $self->has_next_phrase($phrases); #Found them all!

            ($key, $regex) = $self->next_phrase($phrases);
        } else {
            next; #just to be explicit
        }
    }

    $fh->close;

    if($self->has_next_phrase($phrases)) {
        $self->error_message('Failed to find all phrases in file (missing: ' . join(", ", @$phrases) . '.');
        return;
    }

    $self->data(\%values);
    return \%values;
}

sub trim {
    my $self = shift;
    my $line = shift;

    chomp($line);
    $line =~ s/^\s+//;

    return $line;
}


sub has_next_phrase {
    my $self = shift;
    my $phrases = shift;

    return scalar @$phrases;
}

sub next_phrase {
    my $self = shift;
    my $phrases = shift;

    $self->has_next_phrase($phrases)
        or die('Ran out of phrases');

    return shift(@$phrases), shift(@$phrases);
}

sub phrases {
    my $self = shift;
    return [ #An array to preserve order
        'build_id'                      => 'build id: (\d+)',
        'lane_count'                    => 'lane count: (\d+)',
        'haploid_coverage'              => 'haploid coverage: ((?:\\d+\\.\\d{1,3})|(?:Not Available))',
        'unfiltered_snp_calls'          => 'unfiltered SNP calls: ([\d,]+)',
        'filtered_snp_calls'            => 'filtered SNP calls: ([\d,]+)',
        'unfiltered_diploid_heterozygous_percentage' => 'unfiltered diploid heterozygous %: ((?:\d+\.\d{1,3})|Not Available)',
        'filtered_diploid_heterozygous_percentage' => 'filtered diploid heterozygous %: ((?:\d+\.\d{1,3})|Not Available)',
        'unfiltered_dbsnp_concordance'  => 'unfiltered dbsnp concordance: ((?:\d{1,2}\.\d{1,4}%?)|Not Available)',
        'filtered_dbsnp_concordance'    => 'filtered dbsnp concordance: ((?:\d{1,2}\.\d{1,4}%?)|Not Available)',
    ];
}
