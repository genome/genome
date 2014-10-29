package Genome::Model::Tools::Bed::MergeBed12Junctions;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::Bed::MergeBed12Junctions {
    is => ['Command::V2'],
    has_input => [
        input_files => {
            doc => 'The input BED12 format files of junctions.',
        },
        output_file => {
            is => 'Text',
            doc => 'The output junctions file in BED12 format.',
        },
        bedtools_version => {
            is => 'Text',
            doc => 'The version of BEDTools to use for sort function.',
            is_optional => 1,
            default_value => Genome::Model::Tools::BedTools->default_bedtools_version,
        },
    ],
};

sub help_brief {
    "Tool to merge BED12 files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed merge-bed12-junctions...
EOS
}

sub help_detail {
    return <<EOS
    This script will merge tophat's junction.bed files into a single file with cumulative scores and blockSizes
EOS
}

sub execute {
    my $self = shift;

    my $input_files = $self->input_files;
    unless  (ref($input_files) eq 'ARRAY') {
        my @input_files = split(',',$input_files);
        $input_files = \@input_files;
    }
    
    my $bed12_writer = Genome::Utility::IO::BedWriter->create(
        output => $self->output_file,
        headers => Genome::Utility::IO::BedReader->bed12_headers,
    );
    unless ($bed12_writer) {
        $self->error_message('Failed to open BED12 writer for merged output: '. $self->output_file);
        die($self->error_message);
    }
    my %junctions;
    for my $input_file (@$input_files) {
        my $bed12_reader = Genome::Utility::IO::BedReader->create(
            input => $input_file,
            headers => Genome::Utility::IO::BedReader->bed12_headers,
        );
        while (my $bed12_data = $bed12_reader->next) {
            my $blocks = $bed12_data->{'blockCount'};
            unless ($blocks == 2) {
                die('This merge relies on the BED12 format files being output from TopHat and the intervals represented being exon-exon splice junctions!');
            }
            my $start = $bed12_data->{'chromStart'};
            my @b_sizes = split(/,/, $bed12_data->{'blockSizes'});
            my @b_offsets = split(/,/, $bed12_data->{'blockStarts'});
            for (my $b = 1; $b < $blocks; $b++) {
                my $chr = $bed12_data->{'chrom'};
                my $left = $start + $b_offsets[$b-1] + $b_sizes[$b-1];
                my $right = $start + $b_offsets[$b] + 1;
                my $strand = $bed12_data->{'strand'};
                my $key = $chr .':'. $left .'-'. $right .'('. $strand .')';
                push @{$junctions{$key}}, $bed12_data;
            }
        }
    }
    my $tmp_unsorted_bed12_file = Genome::Sys->create_temp_file_path();
    my $tmp_unsorted_bed12_writer = Genome::Utility::IO::BedWriter->create(
        output => $tmp_unsorted_bed12_file,
        headers => Genome::Utility::IO::BedReader->bed12_headers,
    );
    for my $junction_id (keys %junctions) {
        my $read_count;
        my @matches = @{$junctions{$junction_id}};
        my %bed12_data;
        if (@matches > 1) {
            my ($largest_left_overlap) = sort _by_left_block_size (@matches);
            my ($left_size,undef) = split(',',$largest_left_overlap->{blockSizes});
            my ($largest_right_overlap) = sort _by_right_block_size (@matches);
            my (undef,$right_size) = split(',',$largest_right_overlap->{blockSizes});
            my $read_count;
            for my $match (@matches) {
                $read_count += $match->{score};
            }
            %bed12_data = (
                chrom => $largest_left_overlap->{chrom},
                chromStart => $largest_left_overlap->{chromStart},
                chromEnd => $largest_right_overlap->{chromEnd},
                #TODO: Add names like: JUNC00186640
                name => 'J',
                score => $read_count,
                strand  => $largest_left_overlap->{strand},
                thickStart => $largest_left_overlap->{thickStart},
                thickEnd  => $largest_right_overlap->{thickEnd},
                itemRgb  => $largest_left_overlap->{itemRgb},
                blockCount  => 2,
                blockSizes  => $left_size .','. $right_size,
                blockStarts => $largest_left_overlap->{blockStarts},
            );
        } else {
            %bed12_data = %{$matches[0]};
        }
        $tmp_unsorted_bed12_writer->write_one(\%bed12_data);
    }
    $tmp_unsorted_bed12_writer->output->close;

    #Sort the tmp BED12
    my $tmp_sorted_bed12_file = Genome::Sys->create_temp_file_path();
    my $sort_cmd = Genome::Model::Tools::BedTools::Sort->create(
        input_file => $tmp_unsorted_bed12_file,
        output_file => $tmp_sorted_bed12_file,
        use_version => $self->bedtools_version,
    );
    unless ($sort_cmd->execute) {
        $self->error_message('Failed to sort BED12 file!');
        die($self->error_message);
    }

    # Rename the junctions 
    my $tmp_sorted_bed12_reader = Genome::Utility::IO::BedReader->create(
        input => $tmp_sorted_bed12_file,
        headers => Genome::Utility::IO::BedReader->bed12_headers,
    );
    unless ($tmp_sorted_bed12_reader) {
        $self->error_message('Failed to open tmp sorted BED12 file: '.$tmp_sorted_bed12_file);
        die($self->error_message);
    }
    my $count = 1;
    while (my $data = $tmp_sorted_bed12_reader->next) {
        $data->{name} ='JUNC'. sprintf("%08d",$count++);
        $bed12_writer->write_one($data);
    }
    $bed12_writer->output->close;
    return 1;
}

sub _by_left_block_size {
    my ($a_left_block_size,undef) = split(',',$a->{blockSizes});
    my ($b_left_block_size,undef) = split(',',$b->{blockSizes});
    return $b_left_block_size <=> $a_left_block_size;
}

sub _by_right_block_size {
    my (undef,$a_right_block_size) = split(',',$a->{blockSizes});
    my (undef,$b_right_block_size) = split(',',$b->{blockSizes});
    return $b_right_block_size <=> $a_right_block_size;
}


1;

