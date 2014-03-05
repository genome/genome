package Genome::Db::Tgi::Command::ImportCancerGeneCensus;

use strict;
use warnings;

use Genome;
use Genome::File::Fetch qw(fetch);

class Genome::Db::Tgi::Command::ImportCancerGeneCensus {
    is => 'Command',
    has => [
        output_file => {
            is => "Text",
        },
        data_url => {
            is => 'Text',
            default => "cancer.sanger.ac.uk/cancergenome/assets/cancer_gene_census.tsv",
        },
    ],
};

sub execute {
    my $self = shift;

    my $out = Genome::Sys->open_file_for_writing($self->output_file);

    my $infile = fetch($self->data_url);

    my $transformed_infile = "$infile.transformed";
    my $cmd = "cat $infile | tr  '\n' > $transformed_infile";
    Genome::Sys->shellcmd(cmd => $cmd);

    my $in = Genome::Utility::IO::SeparatedValueReader->create(
        input => $transformed_infile,
        separator => "\t",
    );

    my %genes;

    while (my $line = $in->next()) {
        $genes{$line->{Symbol}}++;
    }

    for my $gene (sort keys %genes) {
        $out->print("$gene\n");
    }
    $out->close;

    return 1;
}

1;

