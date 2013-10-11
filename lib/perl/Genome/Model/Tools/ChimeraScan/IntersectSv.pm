package Genome::Model::Tools::ChimeraScan::IntersectSv;

use strict;
use warnings;

use Sort::Naturally;

use Genome;
use Genome::File::BedPe::Entry;
use Genome::File::BedPe::Reader;

class Genome::Model::Tools::ChimeraScan::IntersectSv {
    is  => 'Command',
    has => [
        filtered_bedpe_file => {
            is => 'Text',
        },
        sv_output_file => {
            is => 'Text',
        },
        output_file => {
            is => 'Text',
        },
        annotation_build_id => {
            is => 'Text',
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            id_by => 'annotation_build_id',
        },
    ],
};


sub execute {
    my $self = shift;
    my $filtered_bedpe_file  = $self->filtered_bedpe_file;
    my $sv_output_file       = $self->sv_output_file;
    my $annot_build          = $self->annotation_build;
    
    #cache transcript
    my %params = (
        reference_build_id => $annot_build->reference_sequence_id,
        data_directory     => $annot_build->_annotation_data_directory,
    );
    my @cached_transcripts = Genome::Transcript->get(%params);

    #gather chimera filter info, chrom1 -> chrom5p, chrom2 -> chrom3p
    my @fusion_headers        = qw(name chrom1 start1 end1 chrom2 start2 end2);
    my @fusion_custom_headers = qw(fusion total_frag spanning_frag);
    my @fusion_all_headers    = (@fusion_headers, @fusion_custom_headers);

    my (%fusion_gene_5p, %fusion_gene_3p);
    my $fusion_reader = Genome::File::BedPe::Reader->new($filtered_bedpe_file);

    while (my $entry = $fusion_reader->next) {
        my ($id_5p, $id_3p) = map{$entry->{custom}->[$_]}qw(3 4);

        my %gene_5p = $self->_get_genes($id_5p, \%params);
        my %gene_3p = $self->_get_genes($id_3p, \%params);
                
        my (%hash, %custom_hash);
        @hash{@fusion_headers} = map{$entry->{$_}}@fusion_headers;
        @custom_hash{@fusion_custom_headers} = map{$entry->{custom}->[$_]}qw(0 5 6);
        my %fusion_info = (%hash, %custom_hash);

        for my $gene_5p_name (keys %gene_5p) {
            my $gene_5p = $gene_5p{$gene_5p_name};
            $gene_5p->{info} = \%fusion_info;
            push @{$fusion_gene_5p{$entry->{chrom1}}}, $gene_5p;
        }

        for my $gene_3p_name (keys %gene_3p) {
            my $gene_3p = $gene_3p{$gene_3p_name};
            $gene_3p->{info} = \%fusion_info;
            push @{$fusion_gene_3p{$entry->{chrom2}}}, $gene_3p;
        }
    }

    my $sv_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input     => $sv_output_file,
        separator => "\t",
    );

    my @sv_headers  = qw(chrA bpA chrB bpB event geneA ensemblIdA geneB ensemblIdB);
    my @all_headers = (@sv_headers, @fusion_all_headers, 'both_partners');

    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        headers   => \@all_headers,
        separator => "\t",
        output    => $self->output_file,
    );

    while (my $line = $sv_reader->next) {
        my %sv_info;
        @sv_info{@sv_headers} = map{$line->{$_}}@sv_headers;
        my ($chrA, $bpA, $chrB, $bpB) = map{$line->{$_}}qw(chrA bpA chrB bpB);

        my %genes = (
            A_5p => $fusion_gene_5p{$chrA},
            A_3p => $fusion_gene_3p{$chrA},
            B_5p => $fusion_gene_5p{$chrB},
            B_3p => $fusion_gene_3p{$chrB},
        );

        my (%pick, %type);

        for my $type (keys %genes) {
            my $bp = $type =~ /^A/ ? $bpA : $bpB;
            for my $gene (@{$genes{$type}}) {
                if ($bp >= $gene->{start} and $bp <= $gene->{stop}) {
                    #name is unique fusion chimerascan cluster ID
                    my $id = $gene->{info}->{name};
                    $type{$id}->{$type}++;
                    $pick{$id} = $gene->{info} unless $pick{$id};
                }
            }
        }

        if (%pick) {
            for my $id (keys %pick) {
                my ($fusion_info, $type) = ($pick{$id}, $type{$id});
                my $both_partners = 'N';
                if (($type->{A_5p} and $type->{B_3p}) or ($type->{A_3p} and $type->{B_5p})) {
                    $both_partners = 'Y';
                }
                my %content = (
                    %sv_info, 
                    %$fusion_info, 
                    both_partners => $both_partners,
                );
                $writer->write_one(\%content);
            }
        }
    }
    return 1;
}


sub _get_genes {
    my ($self, $trans_ids, $params) = @_;
    my @id_info = split /,/, $trans_ids;
    my %gene_info;
   
    for my $id_info (@id_info) {
        my ($trans_id) = $id_info =~ /^(E\S+)\:/;
        if ($trans_id) {
            my $trans = Genome::Transcript->get(transcript_name => $trans_id, %$params);
            if ($trans) {
                my $gene = $trans->gene;
                if ($gene) {
                    my $gene_id = $trans->gene_id;
                    unless ($gene_info{$gene_id}) {
                        $gene_info{$gene_id} = {
                            name  => $gene_id,
                            start => $gene->gene_start,
                            stop  => $gene->gene_stop,
                        };
                    }
                }
                else {
                    $self->warning_message("Failed to get gene from trans id: $trans_id");
                }
            }
            else {
                $self->warning_message("Failed to get transcript from id: $trans_id");
            }
        }
        else {
            $self->warning_message("Failed to get transcript id from $id_info");
        }
    }
    return %gene_info;
}

1;

