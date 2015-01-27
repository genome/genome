package Genome::Model::Tools::Bsmap::MethCalcConversionRate;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Bsmap::MethCalcConversionRate {
    is => 'Command',
    has => [
        snvs_file => {
            is => 'String',
            doc => 'Use snvs.hq file to calculate methylation conversion',
            shell_args_position => 1,
        },
        output_file => {
            is => 'String',
            is_optional => 1,
            doc => 'Output methylation conversion',
            shell_args_position => 2,
        },
    ],
};

sub help_synopsis {
  return <<EOS
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
    my $filename = shift;
    my $chrom = shift;
    my $output_file = shift;

    my $count = 0;
    my $totalreads = 0;
    my $methreads = 0;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        headers => [qw(chr pos strand context ratio eff_CT_count C_count CT_count rev_G_count rev_GA_count CI_lower CI_upper)],
        separator => "\t",
        input => $filename,
    );
    while (my $data = $reader->next()) {
        if(($data->{strand} eq "-" && $data->{context} =~ /.CG../ ) || ($data->{strand} eq "+" && $data->{context} =~ /..CG./)){
            $totalreads = $totalreads + $data->{eff_CT_count};
            $methreads = $methreads + $data->{C_count};
        }
    }
    $reader->input->close();

    my $cfile = $output_file;
    if($chrom eq "MT"){
        print $cfile "\nMethylation conversion based on mtDNA:\n";
    }
    if($chrom eq "lambda"){
        print $cfile "\nMethylation conversion based on lambda:\n";
    }
    print $cfile "Meth reads\t=\t", $methreads, "\n";
    print $cfile "Total reads\t=\t", $totalreads, "\n";
    if ($totalreads != 0) {
        print $cfile "Bisulfite conversion (%)\t=\t", 100-($methreads/$totalreads*100), "\n\n";
    }
}

sub execute {
    my $self = shift;
    my $snvs_file = $self->snvs_file;
    my $output_file = $self->output_file || '-';

    my $cfile = Genome::Sys->open_file_for_writing($output_file);

    # snvs
    if (-s "$snvs_file"){
        bs_rate($snvs_file, "MT", $cfile);
    }
    else {
        die $self->error_message("can't find the snvs file");
    }

    return 1;
}

1;
