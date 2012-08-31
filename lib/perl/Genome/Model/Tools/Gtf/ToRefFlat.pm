package Genome::Model::Tools::Gtf::ToRefFlat;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gtf::ToRefFlat {
    is => ['Genome::Model::Tools::Gtf::Base'],
    has => [
        output_file => {
            is => 'Text',
            doc => 'The output refFlat format file.',
        },
    ],
};

sub execute {
    my $self = shift;

    my $tmp_gene_pred = Genome::Sys->create_temp_file_path();
    unless (Genome::Model::Tools::Gtf::ToGenePred->execute(
        input_gtf_file => $self->input_gtf_file,
        extended => 1,
        output_file => $tmp_gene_pred,
    )) {
        $self->error_message('Failed to convert GTF to genePredExt: '. $self->input_gtf_file);
        return;
    }
    my $gene_pred_fh = Genome::Sys->open_file_for_reading($tmp_gene_pred);
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    while (my $line = $gene_pred_fh->getline) {
        chomp($line);
        my @entry = split("\t",$line);
        my @new_entry;
        push @new_entry, $entry[11];
        for (0 .. 9) {
            push @new_entry, $entry[$_];
        }
        my $string = join("\t",@new_entry);
        print $output_fh $string ."\n";
    }
    $output_fh->close;
    return 1;
};


1;
