package Genome::Db::Tgi::Command::ImportMitelman;

use strict;
use warnings;
use Genome;

class Genome::Db::Tgi::Command::ImportMitelman {
    is => 'Command',
    has => [
        data_url => {
            is => 'Text',
            default => "ftp://ftp1.nci.nih.gov/pub/CGAP/mitelman.tar.gz",
        },
        output_dir => {
            is => 'Text',
        },
    ],
};

sub execute {
    my $self = shift;
    my $temp_dir = Genome::Sys->create_temp_directory;
    my $tarball = "$temp_dir/mitel.tar.gz";
    my $cmd = sprintf(q(curl --silent --insecure --output %s %s), $tarball, $self->data_url);
    Genome::Sys->shellcmd(cmd => $cmd);

    $cmd = "tar -xzf $tarball --directory $temp_dir";
    Genome::Sys->shellcmd(cmd => $cmd);
    my $fusion_gene_list = "$temp_dir/molclingene.dat";
    unless (-e $fusion_gene_list) {
        $self->error_message("No fusion gene list found at $fusion_gene_list");
        return;
    }

    my $in = Genome::Utility::IO::SeparatedValueReader->create(
        input => $fusion_gene_list,
        separator => "\t",
    );

    my %fp;
    my %tp;
    while (my $line = $in->next()) {
        if ($line->{GENE} =~ /\//) {
            my ($five_prime, $three_prime) = split(/\//, $line->{GENE});
            $fp{$five_prime}++;
            $tp{$three_prime}++;
        }
    }

    my $five_prime_file = $self->output_dir."/MitelmanFivePrime.txt";
    my $three_prime_file = $self->output_dir."/MitelmanThreePrime.txt";
    my $fp_out = Genome::Sys->open_file_for_writing($five_prime_file);
    my $tp_out = Genome::Sys->open_file_for_writing($three_prime_file);
    for my $gene (sort keys %fp) {
        $fp_out->print("$gene\n");
    }
    for my $gene (sort keys %tp) {
        $tp_out->print("$gene\n");
    }

    $fp_out->close;
    $tp_out->close;
}

1;

