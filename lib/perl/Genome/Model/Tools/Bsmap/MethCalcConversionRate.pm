# % Last Change: Thu Jan 15 01:21:18 PM 2015 CST
package Genome::Model::Tools::Bsmap::MethCalcConversionRate;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Bsmap::MethCalcConversionRate {
    is => 'Command',
    has => [
        snvs_file => {
            is => 'String',
            is_optional => 1,
            doc => 'Use snvs.hq file to calculate methylation conversion',
        },
        model_id => {
            is => 'String',
            is_optional => 1,
            doc => 'Use genome model ID to calculate methylation conversion',
        },
        output_file => {
            is => 'String',
            is_optional => 1,
            doc => 'Output methylation conversion',
        },
    ],
};

sub help_synopsis {
  return <<EOS
    gmt bsmap meth-calc-conversion-rate --model-id=394d6228a7b5487a9cb0ad0c448b5a44

    gmt bsmap meth-calc-conversion-rate --snvs-file=/gscmnt/gc9016/info/model_data/394d6228a7b5487a9cb0ad0c448b5a44/buildf92b6057072948c2a6056a3ee412d596/variants/snv/meth-ratio-2.74-d41d8cd98f00b204e9800998ecf8427e/MT/snvs.hq

EOS
}

sub help_brief {
    "calculate methylation conversion using Mitochondrial DNA or spiked lambda"
}

sub help_detail {
    "calculate methylation conversion using Mitochondrial DNA or spiked lambda"
}

sub bs_rate {
    my $count = 0;
    my $totalreads = 0;
    my $methreads = 0;

    open (my $fh, '<', $_[0]) or die("Could not open snvs file.");
    while (my $line = <$fh>) {
        if ( $count == 0 ) {
            $count = 1;
            next;
        }
        chomp($line);
        my @F = split("\t",$line);
        if(($F[2] eq "-" && $F[3] =~ /.CG../ ) || ($F[2] eq "+" && $F[3] =~ /..CG./)){
            $totalreads = $totalreads + $F[5];
            $methreads = $methreads + $F[6];
        }
    }
    close($_[0]);

    my $cfile;

    if (defined($_[2])){
        open($cfile, ">>$_[2]") || die "Can't open output file for writing.\n";
        if($_[1] eq "MT"){
            print $cfile "\nMethylation conversion based on mtDNA:\n";
        }
        if($_[1] eq "lambda"){
            print $cfile "\nMethylation conversion based on lambda:\n";
        }
        print $cfile "Meth reads\t=\t", $methreads, "\n";
        print $cfile "Total reads\t=\t", $totalreads, "\n";
        if ($totalreads != 0) {
            print $cfile "Bisulfite conversion (%)\t=\t", 100-($methreads/$totalreads*100), "\n\n";
        }
        close $cfile;
    }
    else {
        print "Meth reads\t=\t", $methreads, "\n";
        print "Total reads\t=\t", $totalreads, "\n";
        if ($totalreads != 0) {
            print "Bisulfite conversion (%)\t=\t", 100-($methreads/$totalreads*100), "\n\n";
        }
    }
}

sub execute {
    my $self = shift;
    my $snvs_file = $self->snvs_file;
    my $model_id = $self->model_id;
    my $output_file = $self->output_file;

    if (!defined($snvs_file) && !defined($model_id) ){
        die("Must provide snvs file OR model ID.\n");
    }

    if (defined($snvs_file) && defined($model_id) ){
        die("Must provide snvs file OR model ID.\n");
    }

    my $cfile;

    if (defined($snvs_file)){
        # snvs
        if (-s "$snvs_file"){
            if (defined($output_file)){
                open($cfile, ">$output_file") || die "Can't open output file for writing.\n";
                close($cfile);
                bs_rate($snvs_file, "MT", $output_file);
            }
            else {
                bs_rate($snvs_file, "MT");
            }
        }
        else {
            $self->error_message("can't find the snvs file");
        }
    }

    if (defined($model_id)){
        # get model IDs
        my $model = Genome::Model->get($model_id);
        # get the bam paths of the last succeeded build
        my $dir = $model->last_succeeded_build->data_directory;
        # flagstat 
        my @flagstat = glob("$dir/alignments/*flagstat");
        my (@field, $total, $duplicates, $mapped, $properly);
        for my $flagstat (@flagstat) {
            if (-s "$flagstat"){
                if (defined($output_file)){
                    open($cfile, ">$output_file") || die "Can't open output file for writing.\n";
                    print $cfile "\nMethylation alignment status:\n";
                    print $cfile $flagstat, "\n";

					my $flagstat_data = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat);
					print $cfile "Total read\t=\t", $flagstat_data->{total_reads}, "\n";

					print $cfile "Duplicates\t=\t", $flagstat_data->{reads_marked_duplicates}, "\n";
					if ($flagstat_data->{total_reads} != 0) {
						my $dupe_rate = $flagstat_data->{reads_marked_duplicates} / $flagstat_data->{total_reads} * 100;
						print $cfile "Duplicates rate (%)\t=\t", $dupe_rate, "\n";
					}

					print $cfile "Mapped read\t=\t", $flagstat_data->{reads_mapped}, "\n";
					if ($flagstat_data->{total_reads} != 0) {
						print $cfile "Mapped rate (%)\t=\t", $flagstat_data->{reads_mapped_percentage}, "\n";
					}

					print $cfile "Properly paired\t=\t", $flagstat_data->{reads_mapped_in_proper_pairs}, "\n";
					if ($flagstat_data->{total_reads} != 0) {
						print $cfile "Properly paired rate (%)\t=\t", $flagstat_data->{reads_mapped_in_proper_pairs_percentage}, "\n";
					}

                    close $cfile;
                }
                else {
                    print "\nMethylation alignment status:\n";

					my $flagstat_data = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat);
					print "Total read\t=\t", $flagstat_data->{total_reads}, "\n";

					print "Duplicates\t=\t", $flagstat_data->{reads_marked_duplicates}, "\n";
					if ($flagstat_data->{total_reads} != 0) {
						my $dupe_rate = $flagstat_data->{reads_marked_duplicates} / $flagstat_data->{total_reads} * 100;
						print "Duplicates rate (%)\t=\t", $dupe_rate, "\n";
					}

					print "Mapped read\t=\t", $flagstat_data->{reads_mapped}, "\n";
					if ($flagstat_data->{total_reads} != 0) {
						print "Mapped rate (%)\t=\t", $flagstat_data->{reads_mapped_percentage}, "\n";
					}

					print "Properly paired\t=\t", $flagstat_data->{reads_mapped_in_proper_pairs}, "\n";
					if ($flagstat_data->{total_reads} != 0) {
						print "Properly paired rate (%)\t=\t", $flagstat_data->{reads_mapped_in_proper_pairs_percentage}, "\n";
					}
                }
            }
            else {
                $self->error_message("can't find flagstat file");
            }
        }

		my %cases = (
			MT => { glob => "$dir/variants/snv/meth-ratio-*/MT/snvs.hq", name => "mtDNA" },
			lambda => { glob => "$dir/variants/snv/meth-ratio-*/gi_9626243_ref_NC_001416.1_/snvs.hq", name => "lambda" },
		);
		for my $chrom (keys %cases) {
			my ($glob, $name) = ($cases{$chrom}{glob}, $cases{$chrom}{name});

			my @file = glob($glob);
			for my $file (@file) {
				if (-s "$file"){
					if (defined($output_file)){
						bs_rate($file, $chrom, $output_file);
					}
					else {
						print "\nMethylation conversion based on $name:\n";
						bs_rate($file, $chrom);
					}
				}
				else {
					$self->error_message("can't find $name snvs file");
				}
			}
		}
    }
	return 1;

}

1;
