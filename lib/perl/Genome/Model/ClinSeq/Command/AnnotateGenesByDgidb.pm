package Genome::Model::ClinSeq::Command::AnnotateGenesByDgidb;

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::ClinSeq::Command::AnnotateGenesByDgidb {
    is => 'Command::V2',
    has_input => [
        input_file => {
            is  => 'FilesystemPath',
            doc => 'tsv formatted file to added DGIDB annotations to',
        },
        gene_name_regex => {
            is  => 'Text',
            doc => 'regular expression of column names containing gene names/symbols'
        },
    ],
    has_output => [
        output_file => {
            is  => 'FilesystemPath',            
            doc => 'result file of DGIDB output',
        },
    ],
    doc => 'take a tsv file as input, extract gene list and annotate genes against DGIDB',
};

sub help_synopsis {
    return <<EOS
genome model clin-seq annotate-genes-by-dgidb --input-file=/gscmnt/ams1108/info/model_data/2888708572/build134369422/AML103/snv/wgs_exome/snvs.hq.tier1.v1.annotated.compact.readcounts.tsv --gene-name-column=mapped_gene_name --outdir=/myhome/dir
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
    my $gene_name_regex = $self->gene_name_regex;

    $self->status_message("Adding DGIDB gene annotation to $infile using genes in this gene name column: $gene_name_regex");

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input     => $infile,
        separator => "\t",
    );

    my $gene_list = __PACKAGE__->convert($reader, $gene_name_regex);
    die $self->error_message("gene list is empty. check $infile") unless $gene_list;

    my ($outfile_name, $dir) = fileparse($infile);
    $outfile_name .= '.dgidb';
    my $output_file  = $dir . $outfile_name;

    my $cmd =Genome::Model::Tools::Dgidb::QueryGene->create(
        output_file => $output_file,
        genes       => $gene_list
    );
    
    unless ($cmd->execute) {
        die $self->error_message("Failed to run gmt dgidb query-gene on gene list: $gene_list");
    }

    $self->output_file($output_file);
    return 1;
}


sub convert {
    my ($class, $reader, $column_regex) = @_;
    my %list;

    while (my $data = $reader->next) {
        my $flag = 0;
        for my $column_name (keys %$data) {
            if ($column_name =~ /$column_regex/) {
                $list{$data->{$column_name}} = 1;
                $flag++;
            }
        }
        die "$column_regex not found in file" unless $flag;
    }

    return join ',', sort keys %list;
}


1;

