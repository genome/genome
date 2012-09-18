package Genome::Model::Tools::Cufflinks::Spia;

use strict;
use warnings;

use Genome;

use Genome::Model::ClinSeq::Util;

class Genome::Model::Tools::Cufflinks::Spia {
    is => 'Genome::Model::Tools::Cufflinks',
    has_input => [
        differential_expression_path => {
            doc => 'A gene_exp.diff file output from cuffdiff.',
        },
        output_tsv_file => {
            doc => 'The output from SPIA as a tsv file.',
        },
        species => {
            default_value => 'human',
            valid_values => ['human','mouse'],
        },
    ],
};

sub execute {
    my $self = shift;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->differential_expression_path,
        separator => "\t",
    );
    unless ($reader) {
        die('Failed to open input file: '. $self->differential_expression_path);
    }
    my @diff_headers = @{$reader->headers};
    my @output_headers = ('entrez_id',@diff_headers);
    my $tmp_cuffdiff_file = Genome::Sys->create_temp_file_path();
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $tmp_cuffdiff_file,
        separator => "\t",
        headers => \@output_headers,
    );
    unless ($writer) {
        die('Failed to open output file: '. $self->output_path);
    }
    my $edata = Genome::Model::ClinSeq::Util::loadEntrezEnsemblData(-species => $self->species);

    my $symbols_map = $edata->{'symbols'};
    my $ensembl_map = $edata->{'ensembl_ids'};

    my %entrez_data;
    my %not_found;
    while (my $data = $reader->next) {
        my @genes = split(',',$data->{gene});
        for my $gene (@genes) {
            my $gene_name = Genome::Model::ClinSeq::Util::fixGeneName(
                -gene => $gene,
                -verbose => 0,
                -entrez_ensembl_data => $edata
            );
            my @entrez_ids;
            if ($symbols_map->{$gene_name}) {
                my $entrez_ids = $symbols_map->{$gene_name}->{entrez_ids};
                @entrez_ids = keys %{$entrez_ids};
            } elsif ($ensembl_map->{$gene_name}) {
                my $entrez_ids = $ensembl_map->{$gene_name}->{entrez_ids};
                @entrez_ids = keys %{$entrez_ids};
            } else {
                push @{$not_found{$gene_name}}, $data;
            }
            for my $entrez_id (@entrez_ids) {
                my %new_data = %{$data};
                $new_data{gene} = $gene_name;
                $new_data{entrez_id} = $entrez_id;
                push @{$entrez_data{$entrez_id}}, \%new_data;
            }
        }
    }
    
    for my $entrez_id (sort {$a <=> $b} keys %entrez_data) {
        my @lines = @{$entrez_data{$entrez_id}};
        # TODO: Create a better filter for the "most significant" comparison per gene
        # TODO: Alternatively, how do we deal with more than a 1 x 1 comparison
        my $largest_fold_change;
        for my $line (@lines) {
            # TODO: Figure out workaround for infinity log2(fold_change)
            if ($line->{'log2(fold_change)'} =~ /1\.79769e\+308/)  { next; }
            if (!defined($largest_fold_change)) {
                $largest_fold_change = $line;
            } elsif ( abs($largest_fold_change->{'log2(fold_change)'}) < abs($line->{'log2(fold_change)'}) ) {
                $largest_fold_change = $line;
            }
        }
        if ($largest_fold_change) {
            $writer->write_one($largest_fold_change);
        }
    }
    $writer->output->close;
    
    $self->status_message('Not Found : '. scalar(keys %not_found));
    $self->status_message('Entrez Ids : '. scalar(keys %entrez_data));

    my $kegg_species;
    if ($self->species eq 'human') {
        $kegg_species = 'hsa';
    } elsif ($self->species eq 'mouse') {
        $kegg_species = 'mmu';
    }
    my ($r_fh,$r_file) = Genome::Sys->create_temp_file();
    print $r_fh 'library(SPIA)' ."\n";
    print $r_fh 'cuff <- read.delim("'. $tmp_cuffdiff_file .'")'."\n";
    print $r_fh 'sig <- cuff[cuff$significant=="yes",]'."\n";
    print $r_fh 'DE <- sig$log2.fold_change.'."\n";
    print $r_fh 'names(DE) <- sig$entrez_id'."\n";
    print $r_fh 'ALL <- cuff$entrez_id'."\n";
    print $r_fh 'res = spia(de = DE, all = ALL, organism = "'. $kegg_species .'",nB = 2000, plots = FALSE, beta = NULL)'."\n";
    print $r_fh 'write.table(res,file="'. $self->output_path .'",sep="\t",row.names=FALSE,quote=FALSE)' ."\n";
    $r_fh->close;
    my $r_cmd = 'Rscript '. $r_file;
    Genome::Sys->shellcmd(cmd => $r_cmd);
    return 1;
}


1;
