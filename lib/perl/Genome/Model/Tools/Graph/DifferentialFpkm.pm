package Genome::Model::Tools::Graph::DifferentialFpkm;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Graph::DifferentialFpkm {
    is => 'Command',
    has => [
    fpkm_matrix => {
        type => 'String',
        doc => "ClinSeq FPKM Matrix output",
    },
    list_of_genes => {
        type => 'String',
        doc => 'one gene per line file of genes of interest',
        default => 'NCBI-human.combined-annotation/54_36p_v2',
    },
    output_prefix   => { 
        type => 'String',  
        doc => "output prefix for output graphs. 6 genes per graph", 
        is_optional => 1
    },
    png_as_well => {
        type=>"Boolean",
        doc=>"--png as well if you want pdfs and pngs, not just pdfs",
        default=>0,
    },

    ],
    has_optional => [
    ],
};

sub help_brief {
    "Graph Tumor/Normal Differential FPKM",
}

sub help_synopsis {
    return <<"EOS"
gmt graph differential-fpkm --fpkm-matrix=CuffLinks_GeneLevel.tsv --list_of_genes=gene_list --output-prefix=\$PWD/graph
EOS
}

sub help_detail {
    return <<"EOS"
EOS
}

sub execute {
    $DB::single = 1;
    my $self = shift;
    my $gene_fh = Genome::Sys->open_file_for_reading($self->list_of_genes);
    my ($temp_fh, $temp_filename);
    my $count = 0;
    my @genes;


   my $r_script = $INC[0] . "/Genome/Model/Tools/Graph/differential_fpkm.R";
    while(my $line = $gene_fh->getline) {
        if($temp_fh && $temp_fh->opened) {
            $temp_fh->print($line);
        }
        else {
            ($temp_fh, $temp_filename)=Genome::Sys->create_temp_file();
            $temp_fh->print($line);
        }
        chomp($line);
        push @genes, $line;

        $count++;
        if($count == 6) {
            $temp_fh->close;
            $count=0;
            my $out_filename = $self->output_prefix . "_" . join("_", @genes) . ".pdf";
            @genes = ();
            my $fpkm_matrix = $self->fpkm_matrix;
            my $r_cmd = "R --slave --args $temp_filename $fpkm_matrix $out_filename < $r_script";
            $self->status_message("Running $r_cmd");
            Genome::Sys->shellcmd(cmd=>$r_cmd);
            if($self->png_as_well) {
                my $pdf_name = $out_filename;
                $out_filename =~s/pdf/png/;
                my $convert_cmd="convert -density 150 $pdf_name $out_filename";
                `$convert_cmd`;
            }
        }
    }
    if(@genes) {   #flush the last 1-5.  don't judge me, i'll fix this later
        $temp_fh->close;
        my $out_filename = $self->output_prefix . "_" . join("_", @genes) . ".pdf";
        my $fpkm_matrix = $self->fpkm_matrix;
        my $r_cmd = "R --slave --args $temp_filename $fpkm_matrix $out_filename < $r_script";
        $self->status_message("Running $r_cmd");
        Genome::Sys->shellcmd(cmd=>$r_cmd);
        if($self->png_as_well) {
            my $pdf_name = $out_filename;
            $out_filename =~s/pdf/png/;
            my $convert_cmd="convert -density 150 $pdf_name $out_filename";
            `$convert_cmd`;
        }
    }
    return 1; 
}


        1;


