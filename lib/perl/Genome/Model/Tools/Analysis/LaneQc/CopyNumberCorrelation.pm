package Genome::Model::Tools::Analysis::LaneQc::CopyNumberCorrelation;

use warnings;
use strict;
use Genome;
use IO::File;
use List::Util qw(sum);
use Statistics::Descriptive;

class Genome::Model::Tools::Analysis::LaneQc::CopyNumberCorrelation {
    is => 'Genome::Command::Base',
    has_optional => [
        output_file => {
            type => 'FilePath',
            doc => 'output filename',
        },
        copy_number_laneqc_file_glob => {
            type => 'FilePath',
            doc => 'glob string for grabbing copy-number laneqc files to compare',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'instrument data to correlate lane-qc for',
        },
        lane_qc_models => {
            is => 'Genome::Model::ReferenceAlignment',
            is_many => 1,
            doc => 'lane qc models to specifically use',
        },
    ],
};

sub help_brief {
    "Script to create a correlation matrix from copy-number lane-qc data."
}
sub help_detail {
    "Script to create a correlation matrix from copy-number lane-qc data."
}

sub resolve_copy_number_laneqc_files {
    my $self = shift;
    if    ( $self->copy_number_laneqc_file_glob && !$self->instrument_data && !$self->lane_qc_models) {
        my @files = sort glob($self->copy_number_laneqc_file_glob);
        return @files;
    }
    elsif (!$self->copy_number_laneqc_file_glob &&  $self->instrument_data && !$self->lane_qc_models) {
        my @files;
        for my $instrument_data ($self->instrument_data) {
            my $dir = $instrument_data->lane_qc_dir;
            print "Instrument data " . $instrument_data->__display_name__ . " (" . $instrument_data->id . ") " . ($dir ? 'has' : 'is missing') . " lane QC.\n";
            next unless ($dir);
            push @files, sort glob("$dir/*.cnqc");
        }
        print "\n";
        return @files;
    }
    elsif (!$self->copy_number_laneqc_file_glob && !$self->instrument_data &&  $self->lane_qc_models) {
        my @files;
        for my $model ($self->lane_qc_models) {
            unless ($model->is_lane_qc) {
                print $model->__display_name__ . ' is not a lane QC model.' . "\n";
                next;
            }

            my $build = $model->last_succeeded_build;
            unless ($build) {
                print $model->__display_name__ . ' does not have a succeeded build.' . "\n";
                next;
            }

            my $dir = $build->qc_directory;
            unless ($dir && -d $dir) {
                print $model->__display_name__ . ' does not have a qc direcotry.' . "\n";
                next;
            }

            push @files, sort glob("$dir/*.cnqc");
        }
        print "\n";
        return @files;
    }
    else {
        print STDERR "ERROR: Cannot provide both instrument data and files, please only use one option.\n";
    }
}

sub execute {
    my $self = shift;

    #parse inputs
    my @cnfiles = $self->resolve_copy_number_laneqc_files;
    my $num_files = $#cnfiles;

    my $outfh;
    if ($self->output_file) {
        my $outfile = $self->output_file;
        $outfh = new IO::File $outfile, "w";
        unless ($outfh) {
            die $self->error_message("Failed to open $outfile for writing.");
        }
    }
    else {
        $outfh = *STDOUT;
    }

    #print outfile headers
    print $outfh "File1\tFile2\tCommon_Probes\tCorrelation_coefficient(max=1)\n";

    #Check that files are reasonably similar
    my $standard_wc;
    for my $file (@cnfiles) {
        my $wc_call = `wc -l $file`;
        my ($wc) = $wc_call =~ m/^(\d+)\s+\w+$/;
        unless ($standard_wc) { 
            $standard_wc = $wc; 
            next;
        }
        my $wc_diff = $standard_wc - $wc;
        $wc_diff = abs($wc_diff);
        if ($wc_diff > 100) {
            $self->status_message("Files have largely varied wordcounts (diff>100) - just letting you know in case this is of concern.");
        }
    }

    #Load a hash with the values from all of the files
    my %data;
    for my $file (@cnfiles) {
        unless (-s $file) {
            $self->error_message("File has no size! ($file)");
            return;
        }
        my $fh = new IO::File $file,"r";
        while (my $line = $fh->getline) {
            next if $line =~ m/(^#|CHR)/;
            chomp $line;
            my ($chr,$pos,$rc,$cn) = split /\t/,$line;
            $data{$file}{$chr}{$pos} = $cn;
        }
    }

    #Loop through copy-number to write correlation output file
    my %corr_matrix;
    for my $i1 (0..$num_files) {
        my $f1 = $cnfiles[$i1];
        my $loop2index = $i1 + 1;
        for my $i2 ($loop2index..$num_files) {
            my $f2 = $cnfiles[$i2];
            next if $f1 eq $f2;

            #find common probes
            my (@f1_common,@f2_common);
            my $f1_common = \@f1_common;
            my $f2_common = \@f2_common;
            ($f1_common,$f2_common) = $self->find_common_probes(\%data,$f1,$f2,$f1_common,$f2_common);
            unless (@f1_common) {
                $self->error_message('Empty @f1_common.');
                return;
            }
            unless (@f2_common) {
                $self->error_message('Empty @f2_common.');
                return;
            }
            if ($#f1_common ne $#f2_common) {
                $self->error_message("Common probe numbers don't match for $f1 and $f2.");
                return;
            }

            #find means and standard deviations of common probes
            my $stats1 = Statistics::Descriptive::Full->new();
            my $stats2 = Statistics::Descriptive::Full->new();
            $stats1->add_data(@f1_common);
            $stats2->add_data(@f2_common);
            my $mean1 = $stats1->mean();
            my $mean2 = $stats2->mean();
            my $std1 = $stats1->standard_deviation();
            my $std2 = $stats2->standard_deviation();
            #correlation denominator = $std1*$std2
            my $corr_denominator = $std1 * $std2;

            #divide data from common probes by the means of the arrays
            #and multiply them together to start numerator calculation
            my @numerator_array;
            for (my $i=0; $i<@f1_common; $i++) {
                $f1_common[$i] -= $mean1;
                $f2_common[$i] -= $mean2;
                $numerator_array[$i] = $f1_common[$i] * $f2_common[$i];
                #This works since the arrays were required to be equal lengths above.
            }
            #finish numerator:
            my $corr_numerator = sum(@numerator_array);

            #print output:
            my $corr = $corr_numerator / $corr_denominator;
            my $num_common_probes = scalar @f1_common;
            $corr /= ($num_common_probes-1);
            my $outline = join("\t",$f1,$f2,$num_common_probes,"$corr\n");
            print $outfh $outline;

        }#end, f2 loop
    }#end, f1 loop

    return 1;
}

sub find_common_probes {
    my ($self,$data,$f1,$f2,$f1comref,$f2comref) = @_;
    #recall that $data{$file}{$chr}{$pos} = $cn;
    for my $chr (keys %{$data->{$f1}}) {
        for my $pos (keys %{$data->{$f1}{$chr}}) {
            if (exists $data->{$f2}{$chr}{$pos}) {
                push @$f1comref,$data->{$f1}{$chr}{$pos};
                push @$f2comref,$data->{$f2}{$chr}{$pos};
            }
        }
    }
    return ($f1comref,$f2comref);
}

1;
