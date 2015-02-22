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
        use_tier_1 => {
            is => 'Boolean',
            doc => 'Include tier 1 variants in the variant list',
            default_value => 1,
        },
        use_tier_2 => {
            is => 'Boolean',
            doc => 'Include tier 2 variants in the variant list',
            default_value => 0,
        },
        use_tier_3 => {
            is => 'Boolean',
            doc => 'Include tier 3 variants in the variant list',
            default_value => 0,
        },
        use_tier_4 => {
            is => 'Boolean',
            doc => 'Include tier 4 variants in the variant list',
            default_value => 0,
        },
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
        validate_svs => {
            type => 'Boolean',
            default => 1,
            doc => "Whether or not to include svs in validation list",
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
        regions_of_interest => {
            is => 'Genome::FeatureList',
            doc => 'Lists of regions to include in the validation list',
            is_many => 1,
        },
        gene_black_lists => {
            is => 'File',
            doc => 'Lists of genes to exclude from the validation list.  Gene symbols must match symbols in annotation file',
            is_many => 1,
        },
        variant_black_lists => {
            is => 'Genome::FeatureList',
            is_many => 1,
            doc => 'Lists of variants in the bed files to exclude from the validation list.',
        },
        additional_snv_lists => {
            is => 'Genome::FeatureList',
            doc => 'Lists of additional snv variants to include in the validation list',
            is_many => 1,
        },
        additional_indel_lists => {
            is => 'Genome::FeatureList',
            doc => 'Lists of additional indel variants to include in the validation list',
            is_many => 1,
        },
        additional_sv_lists => {
            is => 'Genome::FeatureList',
            doc => 'Lists of additional structural variants to include in the validation list',
            is_many => 1,
        },
    ],
    has_input_output => [
        significant_variant_list => {
            is => 'Path',
            doc => 'File containing list of compiled variants',
        },
    ],
    has_param => [
        lsf_resource => {
            default => "-M 14000000 -R 'select[mem>14000] rusage[mem=14000]'",
        },
    ],
};

sub execute {
    my $self = shift;
    my $genes_to_include_bed_file = Genome::Sys->create_temp_file_path;
    my %black_list;
    my %variant_black_list;
    if ($self->gene_black_lists) {
        foreach my $list ($self->gene_black_lists) {
            my $in = Genome::Sys->open_file_for_reading($list);
            while (my $line = <$in>) {
                chomp $line;
                my @fields = split /\t/, $line;
                $black_list{$fields[0]} = 1;
            }
        }
    }
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
            if ($fields[$headers{"FDR LRT"}] <= $fdr_cutoff and !$black_list{$fields[$headers{"#Gene"}]}) {
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
        my $sort_cmd = Genome::Model::Tools::Joinx::Sort->execute(
            input_files => [$gene_file_path],
            unique => 1,
            output_file => $sorted_gene_bed,
        );
        unless ($sort_cmd and $sort_cmd->result) {
            $self->error_message("Joinx sort failed");
            return;
        }

        my $merge_cmd = Genome::Model::Tools::BedTools::Merge->execute(
            input_file => $sorted_gene_bed,
            output_file => $genes_to_include_bed_file,
        );
        unless ($merge_cmd and $merge_cmd->result) {
            $self->error_message("MergeBed failed");
            return;
        }
    }
    if ($self->validation_list_output_format eq 'nimblegen') {
        my @snv_files;
        my @indel_files;
        my @sv_files;
        my $input_file = Genome::Sys->create_temp_file_path;
        my $bed_version = 2;
        my $anno_version = 1;

        my @tiers_to_use;
        if ($self->use_tier_1) {
            push @tiers_to_use, 1;
        }
        if ($self->use_tier_2) {
            push @tiers_to_use, 2;
        }
        if ($self->use_tier_3) {
            push @tiers_to_use, 3;
        }
        if ($self->use_tier_4) {
            push @tiers_to_use, 4;
        }

        foreach my $build ($self->somatic_variation_builds) {
            foreach my $tier (@tiers_to_use){
                my $anno = $build->data_set_path("effects/snvs.hq.tier$tier",$anno_version,"annotated.top");
                my $bed = $build->data_set_path("effects/snvs.hq.novel.tier$tier",$bed_version,"bed");
                if (-e $anno) {
                    my $filtered_file = $self->_filtered_file_based_on_black_list(\%black_list, $anno);
                    push @snv_files, $filtered_file;
                }
                elsif (-e $bed) {#TODO: need to shift to 1-based start?
                    my $filtered_file = $self->_filtered_file_based_on_variant_black_list($bed);
                    push @snv_files, $filtered_file;
                }
                $anno = $build->data_set_path("effects/indels.hq.tier$tier",$anno_version,"annotated.top");
                $bed = $build->data_set_path("effects/indels.hq.novel.tier$tier",$bed_version,"bed");
                if (-e $anno) {
                    my $filtered_file = $self->_filtered_file_based_on_black_list(\%black_list, $anno);
                    push @indel_files, $filtered_file;
                }
                elsif (-e $bed) {
                    my $filtered_file = $self->_filtered_file_based_on_variant_black_list($bed);
                    push @indel_files, $filtered_file;
                }
            }
            if ($self->validate_svs) {
                my $bed = $build->data_set_path("variants/svs.hq",$bed_version,'bed');
                if (-e $bed) {
                    push @sv_files, $bed;
                }
            }
        }



        if ($self->additional_snv_lists) {
            foreach my $additional_snvs ($self->additional_snv_lists) {
                push @snv_files, $additional_snvs;
            }
        }
        if ($self->additional_indel_lists) {
            foreach my $additional_indels ($self->additional_indel_lists) {
                push @indel_files, $additional_indels;
            }
        }
        if ($self->additional_sv_lists) {
            foreach my $additional_svs ($self->additional_sv_lists) {
                push @sv_files, $additional_svs;
            }
        }

        if (-s $genes_to_include_bed_file) {
            push @indel_files, $genes_to_include_bed_file;
        }
        if ($self->regions_of_interest) {
            foreach my $feature_list ($self->regions_of_interest) {
                push @indel_files, $feature_list->file_path;
            }
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
        unless($nimblegen_design and $nimblegen_design->result) {
            $self->error_message("Failed to generate nimblegen design");
            return;
        };
    }
    return 1;
}

sub _filtered_file_based_on_variant_black_list {
    my ($self, $file_name) = @_;
    my $filtered_file;
    if ($self->variant_black_lists) {
        foreach my $variant_black_list ($self->variant_black_lists) {
            $filtered_file = Genome::Sys->create_temp_file_path;
            my $intersect_cmd = Genome::Model::Tools::Joinx::Intersect->execute(
                input_file_a => $file_name,
                input_file_b => $variant_black_list->file_path,
                miss_a_file => $filtered_file,
            );
            unless ($intersect_cmd and $intersect_cmd->result) {
                $self->error_message("Failed to filter black listed variants");
                return;
            }
            $file_name = $filtered_file;
        }
    }
    else {
        $filtered_file = $file_name;
    }
    return $filtered_file;
}

sub _filtered_file_based_on_black_list {
    my ($self, $hash, $file_name) = @_;
    my ($filtered_fh, $filtered_file) = Genome::Sys->create_temp_file;
    my $anno_in = Genome::Sys->open_file_for_reading($file_name);
    while(my $line = <$anno_in>) {
        chomp $line;
        if ($line =~ /^#/) {
            next;
        }
        my @fields = split /\t/, $line;
        unless ($hash->{$fields[6]}) {
            $filtered_fh->print($line."\n");
        }
    }
    $filtered_fh->close;
    return $filtered_file;
}
1;

