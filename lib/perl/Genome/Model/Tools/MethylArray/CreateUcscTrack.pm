package Genome::Model::Tools::MethylArray::CreateUcscTrack;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::MethylArray::CreateUcscTrack {
    is => ['Genome::Model::Tools::MethylArray::Base'],
    has => {
        array_meta_csv_file => {
            is => 'Text',
            doc => 'The Illumina Methyl Array CSV file containing meta-data for each probe',
            default_value =>'/gscmnt/sata132/techd/solexa/jwalker/Methyl-Seq/Infinium_HumanMethylation27_BeadChip/HumanMethylation450_15017482_v.1.1_removeMeta.csv',
            is_optional => 1,
        },
        ucsc_track_file => {
            is => 'Text',
            doc => 'The output BED15 format microarray track for import to UCSC as a custom track.',
        },
        methylation_tsv_file => {
            is => 'Text',
            doc => 'A tab separated value file containing probes and associated methylation values.',
        },
        array_description => {
            is => 'Text',
            doc => 'The description of the array to display in the UCSC browser',
            default_value => 'Microarray custom track',
            is_optional => 1,
        },
        array_name => {
            is => 'Text',
            doc => 'The name of the array to display in the UCSC browser',
            default_value => 'Microarray',
            is_optional => 1,
        },
        experiment_step => {
            is => 'Text',
            doc => 'The experiment step value used for color gradients',
            default_value => '0.5',
            is_optional => 1,
        },
        experiment_scale => {
            is => 'Text',
            doc => 'The experiment scale value used for color gradients',
            default_value => '6',
            is_optional => 1,
        },
        experiment_color => {
            is => 'Text',
            doc => 'The experiment color scheme',
            default_value => 'redBlueOnWhite',
            valid_values => ['redGreen', 'redBlue', 'yellowBlue', 'redBlueOnWhite', 'redBlueOnYellow'],
            is_optional => 1,
        },
        track_visibility => {
            is => 'Integer',
            doc => 'Defines the initial display mode of the annotation track. Values for display_mode include: 0 - hide, 1 - dense, 2 - full, 3 - pack, and 4 - squish.',
            default_value => 2,
            is_optional => 1,
        },
    },
};

sub execute {
    my $self = shift;
    
    my $meta_csv_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->array_meta_csv_file,
    );
    my %probes;
    while (my $data = $meta_csv_reader->next){
        my $probe_name = $data->{IlmnID};
        my $hg19_coordinate = $data->{'MAPINFO'};
        if (!defined($hg19_coordinate) || $hg19_coordinate eq '') {
            next;
            #die(Data::Dumper::Dumper($data));
        }
        $probes{$probe_name}{start} = ($hg19_coordinate-1);
        $probes{$probe_name}{end} = $hg19_coordinate;
        $probes{$probe_name}{chr} = $data->{CHR};
    }
    my $methyl_reader = Genome::Utility::IO::SeparatedValueReader->create(
        is_regex => 1,
        input => $self->methylation_tsv_file,
        separator => '\s+',
    );

    my @headers = @{$methyl_reader->headers};
    my @sample_names;
    for my $header (@headers) {
        if ($header eq '') { next; }
        push @sample_names, $header;
    }
    my $exp_names = join(',',@sample_names) .',';

    my $bed15_fh = Genome::Sys->open_file_for_writing($self->ucsc_track_file);
    unless ($bed15_fh) {
        die('Failed to open output file:'. $self->ucsc_track_file);
    }

    print $bed15_fh 'track type="array" visibility='. $self->track_visibility .' expColor="'. $self->experiment_color .'" expScale='. $self->experiment_scale .' expStep='. $self->experiment_step .' expNames="'. $exp_names .'" name="'. $self->array_name .'" description="'. $self->array_description .'"' . "\n";

    my $samples = scalar(@sample_names);
    # The id string
    my $ids_string = '';
    for (0 .. ($samples - 1)) {
        $ids_string .= $_ .',';
    }
    while ( my $report_data = $methyl_reader->next) {
        my $probe_name = $report_data->{''};
        my $data = $probes{$probe_name};
        unless ($data) {
            next;
            #die('Missing meta-data for probe: '. $probe_name);
        }
        my $start = $probes{$probe_name}{start};
        my $end = $probes{$probe_name}{end};
        my $chr = $probes{$probe_name}{chr};

        my @scores;
        for my $sample (@sample_names) {
            my $short_num = sprintf("%.04f",$report_data->{$sample});
            push @scores, $short_num;
        }
        # The score string
        my $score_string = join(',',@scores) .',';

        my $string = 'chr'. $chr ."\t". #chrom
            $start ."\t". #chromStart
                $end ."\t". #chromEnd
                    $probe_name ."\t". #name
                        "500\t". #score
                            ".\t". #strand
                                $start ."\t". #thickStart
                                    $end ."\t". #thickEnd
                                        "0\t". #reserved
                                            "1\t". #blockCount
                                                ($end - $start) .",\t". #blockSizes
                                                    "0,\t". #chromStarts
                                                        $samples ."\t". #expCount
                                                            $ids_string ."\t". #expIds
                                                                $score_string ."\n"; #expScores
        print $bed15_fh $string;
    }
    $bed15_fh->close;
    return 1;
}

