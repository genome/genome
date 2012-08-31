package Genome::Model::MutationalSignificance::Command::CompileValidationList;

use strict;
use warnings;
use Genome;

class Genome::Model::MutationalSignificance::Command::CompileValidationList {
    is => ['Command::V2'],
    has_input => [
        validation_list_output_format => {
            is => 'Text',
            default_value => 'nimblegen',
            doc => 'Desired format for the validation list',
            valid_values => ['nimblegen'],
        },
        somatic_variation_builds => {
            is => 'Genome::Model::Build::SomaticVariation',
            is_many => 1,
        },
        tiers_to_use => {
            is => 'Number',
            is_many => 1,
            doc => 'A list of tiers to be included in the variant list',
            default_value => [1],
        },
        #genes_to_include => {
        #    type => 'Path',
        #    doc => 'A list of genes to include in the validation list',
        #},
        exon_bed => {
            type => 'Path',
            doc => 'A bed file containing the coordinates of genes that may be of interest',
        },
        snv_indel_span => {
            type => 'Integer',
            default => 0,
            doc => "The number of bases to span the region upstream and downstream of a SNV or Indel locus",
        },
        sv_span => {
            type => 'Integer',
            default => 200,
            doc => "The number of bases to span the region upstream and downstream of an SV locus",
        },
        include_mitochondrial_sites => {
            type => 'Boolean',
            default => 0,
            doc => "Whether or not to remove sites on the mitochondria or non-chromosomal contigs",
        },
        include_unplaced_contig_sites => {
            type => 'Boolean',
            default => 1,
            doc => "Whether or not to remove sites on the unplaced contigs of the chromosome",
        },
        include_y_chrom_sites => {
            type => 'Boolean',
            default => 1,
            doc => "Whether or not to include sites on the Y chromosome in the output (if cases are all female)",
        },
        reference_sequence_build => {
            type => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => "Reference sequence in use, to check chromosomal bounds (E.g. GRCh37-lite-build37)",
        },
    ],
    has_optional_input => [
        significantly_mutated_gene_list => {
            is => 'Path',
            doc => 'Significantly mutated genes to include in validation',
        },
        fdr_cutoff => {
            is => 'Number',
            doc => 'Only include genes with fdr less than this value',
            default => 1,
        },
    ],
    has_input_output => [
        significant_variant_list => {
            is => 'Path',
            doc => 'File containing list of compiled variants',
        },
    ],
};

sub execute {
    my $self = shift;
    my $genes_to_include_bed_file = Genome::Sys->create_temp_file_path;
    if ($self->significantly_mutated_gene_list) {
        my $fdr_cutoff = $self->fdr_cutoff;
        my $fh = Genome::Sys->open_file_for_reading($self->significantly_mutated_gene_list);
        my %gene_hash;
        my %headers;
        while(my $gene = <$fh>) {
            chomp $gene;
            my @fields = split /\t/, $gene;
            if ($gene =~ /^#/) {
                my $count = 0;
                foreach my $header (@fields) {
                    $headers{$header} = $count;
                    $count++;
                }
                next;
            }
            if ($fields[$headers{"FDR LRT"}] <= $fdr_cutoff) {
                $gene_hash{$fields[0]} = 1;
            }
        }
        $fh->close;
        my $exon_bed_fh = Genome::Sys->open_file_for_reading($self->exon_bed);
        my ($gene_file, $gene_file_path) = Genome::Sys->create_temp_file;
        while (my $bed_line = <$exon_bed_fh>) {
            chomp $bed_line;
            my @fields = split /\t/, $bed_line;
            if ($gene_hash{$fields[3]}) {
                my $start = $fields[1]-1;
                $gene_file->print($fields[0]."\t".$start."\t".$fields[2]."\n");
            }
        }
        $gene_file->close;

        my $sorted_gene_bed = Genome::Sys->create_temp_file_path;
        my $sort_rv = Genome::Model::Tools::Joinx::Sort->execute(
            input_files => [$gene_file_path],
            unique => 1,
            output_file => $sorted_gene_bed,
        );
        unless ($sort_rv) {
            $self->error_message("Joinx sort failed");
            return;
        }

        my $sorted_merged_gene_bed = Genome::Sys->create_temp_file_path;
        my $merge_rv = Genome::Model::Tools::BedTools::Merge->execute(
            input_file => $sorted_gene_bed,
            output_file => $genes_to_include_bed_file,
        );
        unless ($merge_rv) {
            $self->error_message("MergeBed failed");
            return;
        }
    }
    if ($self->validation_list_output_format eq 'nimblegen') {
        my @snv_files;
        my @indel_files;
        my @sv_files;
        my $input_file = Genome::Sys->create_temp_file_path;

        foreach my $build ($self->somatic_variation_builds) {
            foreach my $tier ($self->tiers_to_use){
                push @snv_files, $build->data_set_path("effects/snvs.hq.tier$tier",$tier,"annotated.top");
                push @indel_files, $build->data_set_path("effects/indels.hq.tier$tier",$tier,"annotated.top");
                push @sv_files, $build->data_set_path("effects/svs.hq.tier$tier",$tier,"annotated.top");
            }
        }

        if (-s $genes_to_include_bed_file) {
            push @indel_files, $genes_to_include_bed_file;
        }

        my $fh = Genome::Sys->open_file_for_writing($input_file);

        $fh->print("snvs\n");
        foreach my $file (@snv_files) {
            if (-s $file) {
                $fh->print("$file\n");
            }
        }

        $fh->print("indels\n");
        foreach my $file (@indel_files) {
            if (-s $file) {
                $fh->print("$file\n");
            }
        }

        $fh->print("svs\n");
        foreach my $file (@sv_files) {
            if (-s $file) {
                $fh->print("$file\n");
            }
        }

        $fh->close;

        my $nimblegen_design = Genome::Model::Tools::Nimblegen::DesignFromFiles->execute(
            input_file => $input_file,
            output_file => $self->significant_variant_list,
            snv_indel_span => $self->snv_indel_span,
            sv_span => $self->sv_span,
            include_mitochondrial_sites => $self->include_mitochondrial_sites,
            include_unplaced_contig_sites => $self->include_unplaced_contig_sites,
            include_y_chrom_sites => $self->include_y_chrom_sites,
            reference => $self->reference_sequence_build->name,
        );
        unless($nimblegen_design) {
            $self->error_message("Failed to generate nimblegen design");
            return;
        };
    }
    return 1;
}

1;

