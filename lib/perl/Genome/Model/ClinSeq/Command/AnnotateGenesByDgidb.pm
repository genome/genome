package Genome::Model::ClinSeq::Command::AnnotateGenesByDgidb;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::AnnotateGenesByDgidb {
    is => 'Command::V2',
    has_input => [
        input_file => {
            is  => 'FilesystemPath',
            doc => 'tsv formatted file to added DGIDB annotations to',
        },
        gene_name_column => {
            is  => 'Text',
            doc => 'name of column containing gene names/symbols'
        }
    ],
    has_output => [
        output_file => {
            is  => 'FilesystemPath',            
            doc => 'result file of DGIDB output',
        }
    ],
    doc => 'take a tsv file as input, extract gene list and annotate genes against DGIDB',
};

sub help_synopsis {
    return <<EOS
genome model clin-seq annotate-genes-by-dgidb --input-file=/gscmnt/ams1108/info/model_data/2888708572/build134369422/AML103/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv --gene-name-column='mapped_gene_name'
EOS
}


sub help_detail {
    return <<EOS
Generate a DGIDB gene annotation output
EOS
}


sub execute {
    my $self   = shift;
    my $infile = $self->input_file;
    my $gene_name_column = $self->gene_name_column;

    $self->status_message("Adding DGIDB gene annotation to $infile using genes in this gene name column: $gene_name_column");

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input     => $infile,
        separator => "\t",
    );

    my $gene_list = __PACKAGE__->convert($reader, $gene_name_column);
    die $self->error_message("gene list is empty. check $infile") unless $gene_list;

    my $cmd =Genome::Model::Tools::Dgidb::QueryGene->create(
        output_file => $self->output_file,
        genes       => $gene_list
    );
    
    return 1;
}


sub convert {
    my ($class, $reader, $column_name) = @_;
    my %list;

    while (my $data = $reader->next) {
        my $item = $data->{$column_name};
        die "$column_name not found in file" unless $item;
        $list{$item} = 1;
    }

    return join ',', sort keys %list;
}


1;

