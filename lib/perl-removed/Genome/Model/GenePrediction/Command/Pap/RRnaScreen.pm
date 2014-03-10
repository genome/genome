#$Id$

package Genome::Model::GenePrediction::Command::Pap::RRnaScreen;

use strict;
use warnings;

use Workflow;

use Bio::Annotation::DBLink;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;

use Compress::Bzip2;
use English;
use File::Basename;
use File::Spec;
use File::Temp;
use IPC::Run;


class Genome::Model::GenePrediction::Command::Pap::RRnaScreen {
    is  => ['Command::V1'],
    has => [
        fasta_file      => { 
                            is  => 'SCALAR', 
                            doc => 'fasta file name',
                           },
        blast_db => {
            is => 'SCALAR',
            doc => 'blast database for rRNA',
            is_optional => 1,
            example_values => ["/gscmnt/278/analysis/HGMI/rRNA_testing/16s_23srnadb"],
        },
        blast_report => {
                         is          => 'SCALAR',
                         is_optional => 1,
                         doc         => 'instance of File::Temp pointing to raw blast output'
                        },
        bio_seq_feature => { 
                            is          => 'ARRAY',  
                            is_optional => 1,
                            doc         => 'array of Bio::Seq::Feature' 
                           },
        report_save_dir => {
                            is          => 'SCALAR',
                            is_optional => 1,
                            doc         => 'directory to save a copy of the blast report to',
                           },
        dead_genes => {
                        is          => 'ARRAY',
                        is_optional => 1,
                        doc         => 'array of sequence (query) names seen in the results, to mark as \'dead\'',
                       },
    ],
};

1;
