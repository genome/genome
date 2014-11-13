package Genome::Model::Tools::Bmr::CombineClassSummaryFiles;

use warnings;
use strict;

use IO::File;
use Genome;

class Genome::Model::Tools::Bmr::CombineClassSummaryFiles {
    is => 'Command::V2',
    has_input => [
    class_summary_output_dir => {
        is => 'String',
        is_optional => 0,
        doc => 'Directory containing .class_summary files from batch-class-summary',
    },
    output_file => {
        is => 'String',
        is_optional => 0,
        doc => 'The final class summary file for the dataset',
    },
    ]
};

sub help_brief {
    "Combine results from batch-class-summary jobs."
}

sub help_detail {
    "Combine results from batch-class-summary jobs into a single file."
}

sub execute {
    my $self = shift;
    my $summary_dir = $self->class_summary_output_dir;
    $summary_dir = (( $summary_dir =~ m/\/$/ ) ? $summary_dir : "$summary_dir/" );
    my $outfile = $self-> output_file;

    #read dir
    opendir(SUM,$summary_dir);
    my @files = readdir(SUM);
    closedir(SUM);
    @files = grep { /\.class_summary$/ } @files;
    @files = map { $summary_dir . $_ } @files;

    #Merge together all the BMR data from the files
    my %DATA = ();
    for my $file (@files) {
        my $fh = new IO::File $file, "r";
        $fh->getline; #discard the header
        while (my $line = $fh->getline) {
            chomp $line;
            my ($class,$bmr,$cov,$muts) = split /\t/,$line;
            $DATA{$class}{'coverage'} += $cov;
            unless ($DATA{$class}{'mutations'}) {
                $DATA{$class}{'mutations'} = $muts;
            }
        }
        $fh->close;
    }

    #Store merged BMR data to file
    my $outfh = new IO::File $outfile, "w";
    print $outfh "Class\tBMR\tCoverage(Bases)\tNon_Syn_Mutations\n";
    for my $class (sort keys %DATA) {
        my $rate = $DATA{$class}{'mutations'} / $DATA{$class}{'coverage'};
        print $outfh "$class\t$rate\t$DATA{$class}{'coverage'}\t$DATA{$class}{'mutations'}\n";
    }
    $outfh->close;
    return 1;
}

1;
