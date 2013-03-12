package Genome::Model::MutationalSignificance::Command::MergeAnnotations;

use strict;
use warnings;
use Genome;

class Genome::Model::MutationalSignificance::Command::MergeAnnotations {
    is => ['Command::V2'],
    has => [
        tgi_anno_file => {
            is => 'String',
        },
        regulome_db_file => {
            is => 'String',
        },
        annovar_file => {
            is => 'String',
        },
        annovar_columns_to_check => {
            is => 'String',
            is_many => 1,
        },
        output_file => {
            is => 'String',
        },
        annotation_build_id => {
            is => 'Number',
            default => 124434505,
        },
        annotation_build => {
            id_by => 'annotation_build_id',
            is => "Genome::Model::Build::ImportedAnnotation",
        },
    ],
};

sub execute {
    my $self = shift;

    my $tgi = Genome::Sys->open_file_for_reading($self->tgi_anno_file);
    my $reg_db = Genome::Sys->open_file_for_reading($self->regulome_db_file);
    my $annovar = Genome::Sys->open_file_for_reading($self->annovar_file);

    my $out = Genome::Sys->open_file_for_writing($self->output_file);

    my $err_out = Genome::Sys->open_file_for_writing($self->output_file.".errs");
    $err_out->print("Errors\n");
    my @egis = Genome::ExternalGeneId->get(data_directory => $self->annotation_build->data_directory."/annotation_data", reference_build_id => $self->annotation_build->reference_sequence_id, id_type => ["ensembl_default_external_name", "ensembl", "EntrezGene"]);

    my %variants;
    #read in the tgi file first
    while (my $line = <$tgi>) {
        chomp $line;
        my @fields = split /\t/, $line;
        my $chr = $fields[0];
        my $start = $fields[1];
        if (defined $variants{$chr}{$start}->{tgi}) {
            $self->warning_message("Duplicate variant: $chr $start");
        }
        $variants{$chr}{$start}->{tgi} = $line;
    }
    $tgi->close;

    #read in the reg_db file and fill in the info for the entries in the tgi hash
    while (my $line = <$reg_db>) {
        chomp $line;
        next if $line =~ /^#/;
        my @fields = split /\t/, $line;
        my $chr = $fields[0];
        $chr =~ s/chr//;
        my $start = $fields[1]+1;
        my $score = $fields[4];
       
        if (defined $variants{$chr}{$start}) {
            $variants{$chr}{$start}->{regdb_score} = $score;
            if ($score =~ /1/) {
            $err_out->print("regdb_score is 1\n");
                my @hits = split /,\s/, $fields[3];
                my @eQTLs;
                foreach my $hit (@hits) {
                $err_out->print("$hit\n");
                    if ($hit =~ /eQTL/) {
                        $err_out->print("Hit something with an eQTL\n");
                        my @parts = split /\|/, $hit;
                        my $gene_name = $parts[1];
                        my $entrez_id_obj = Genome::ExternalGeneId->get(data_directory => $self->annotation_build->data_directory."/annotation_data", reference_build_id => $self->annotation_build->reference_sequence_id, id_type => "EntrezGene", id_value => $gene_name);
                        if ($entrez_id_obj) {
                            my $default_name_obj = Genome::ExternalGeneId->get(data_directory => $self->annotation_build->data_directory."/annotation_data", reference_build_id => $self->annotation_build->reference_sequence_id, id_type => "ensembl_default_external_name", gene_id => $entrez_id_obj->gene_id);
                            if ($default_name_obj) {
                                $err_out->print("Translating $gene_name to ");
                                $gene_name = $default_name_obj->id_value;
                                $err_out->print("$gene_name\n");
                            }
                            else {
                            $err_out->print("Could not translate $gene_name\n");
                                $self->warning_message("Could not translate $gene_name");
                            }
                        }
                        else {
                        $err_out->print("Could not find gene with entrez name $gene_name\n");
                            $self->warning_message("Could not find gene with entrez name $gene_name for translation.");
                        }
                        $err_out->print("Adding $gene_name to eQTLS\n");
                        push @eQTLs, $gene_name;
                    }
                }
                $variants{$chr}{$start}->{regdb_eQTL} = \@eQTLs;
            }
        }
    }
    $reg_db->close;



    #read in the annovar and fill in further info in the tgi hash
    my $header = <$annovar>;
    chomp $header;
    my @header_fields = split /\t/, $header;
    my %header_hash;
    my $count = 0;
    foreach my $header_field (@header_fields) {
        $header_hash{$header_field} = $count;
        $count++;
    }
    while (my $line = <$annovar>) {
        chomp $line;
        my @fields = split /\t/, $line;
        foreach my $column_name ($self->annovar_columns_to_check) {
            $variants{$fields[0]}{$fields[1]}->{"annovar_$column_name"} = $fields[$header_hash{$column_name}];
        }
    }
    $annovar->close;

    
    #did we get a regulomedb score and annovar for every variant?
    foreach my $chr (keys %variants) {
        foreach my $start (keys %{$variants{$chr}}) {
            my %var = %{$variants{$chr}{$start}};
            unless ($var{regdb_score}) {
                $self->warning_message("No regulome db score for variant $chr:$start");
            }

            my $line = $var{tgi};
            my @fields = split /\t/, $line;
            #build up a list of gene annotations
            my %annotations;
            unless ($fields[21] eq "-") {
                $annotations{$fields[21]}->{sources} = [$fields[22]];
                $annotations{$fields[21]}->{trv_type} = [$fields[13]];
                $annotations{$fields[21]}->{ensembl_id} = $fields[23];
            }
            if ($var{regdb_score} =~ /1/) {
                foreach my $eQTL (@{$var{eQTLs}}){
                    my $ensembl_id = "-";
                    my $other_id_obj = Genome::ExternalGeneId->get(data_directory => $self->annotation_build->data_directory."/annotation_data", reference_build_id => $self->annotation_build->reference_sequence_id, id_type => "ensembl_default_external_name", id_value => $eQTL);
                    if ($other_id_obj) {
                        my $ensembl_id_obj = Genome::ExternalGeneId->get(data_directory => $self->annotation_build->data_directory."/annotation_data", reference_build_id => $self->annotation_build->reference_sequence_id, id_type => "ensembl", gene_id => $other_id_obj->gene_id);
                        if ($ensembl_id_obj) {
                            $ensembl_id = $ensembl_id_obj->id_value;
                        }
                        else {
                            $err_out->print("Could not translate $eQTL into ensembl id\n");
                            $ensembl_id = "-";
                        }
                    }
                    else {
                        $err_out->print("Could not find gene with default name $eQTL\n");
                        $ensembl_id = "-";
                    }
                    push @{$annotations{$eQTL}->{sources}}, "regulomeDB_eQTL";
                    push @{$annotations{$eQTL}->{trv_type}}, "regulatory";
                    if ((not defined $annotations{$eQTL}->{ensembl_id}) or $annotations{$eQTL}->{ensembl_id} eq "-") {
                        $annotations{$eQTL}->{ensembl_id} = $ensembl_id;
                    }
                    else {
                        unless ($annotations{$eQTL}->{ensembl_id} eq $ensembl_id) {
                            $self->warning_message("Ensembl id does not match: ".$annotations{$eQTL}->{ensembl_id}." and $ensembl_id");
                        }
                    }
                }
            }

            if ($var{regdb_score} =~ /1/ || $var{regdb_score} =~ /2/) {
                foreach my $column_name ($self->annovar_columns_to_check) {
                    my $annot = $var{"annovar_$column_name"};
                    unless ($annot) {
                        $self->warning_message("No annovar_$column_name for variant $chr:$start");
                    }
                    next if ($annot =~ "-");
                    $annot =~ s/^Name=//;
                    my @gene_names = split /,/, $annot;
                    foreach my $gene_name (@gene_names) {
                        my $ensembl_id;
                        if ($gene_name =~ /^ENSG/) {
                            $err_out->print("Translating $gene_name\n");
                            $gene_name =~ s/\.[0-9]+$//;
                            $ensembl_id = $gene_name;
                            my $ensembl_id_obj = Genome::ExternalGeneId->get(data_directory => $self->annotation_build->data_directory."/annotation_data", reference_build_id => $self->annotation_build->reference_sequence_id, id_type => "ensembl", id_value => $gene_name);
                            if ($ensembl_id_obj) {
                                my $default_name_obj = Genome::ExternalGeneId->get(data_directory => $self->annotation_build->data_directory."/annotation_data", reference_build_id => $self->annotation_build->reference_sequence_id, id_type => "ensembl_default_external_name", gene_id => $ensembl_id_obj->gene_id);
                                if ($default_name_obj) {
                                    $err_out->print("Translating $gene_name ");
                                    $gene_name = $default_name_obj->id_value;
                                    $err_out->print(" to $gene_name\n");
                                }
                                else {
                                    $err_out->print("Could not translate $gene_name\n");
                                    $self->warning_message("Could not translate $gene_name");
                                }
                            }
                            else {
                                $err_out->print("Could not find gene with ensembl name $gene_name\n");
                                $self->warning_message("Could not find gene with ensembl name $gene_name for translation.");
                            }
                        }
                        else {
                            my @other_id_obj = Genome::ExternalGeneId->get(data_directory => $self->annotation_build->data_directory."/annotation_data", reference_build_id => $self->annotation_build->reference_sequence_id, id_type => "ensembl_default_external_name", id_value => $gene_name);
                            if (@other_id_obj) {
                                my $count = 0;
                                my $str;
                                foreach my $other_obj (@other_id_obj) { #TODO: one will be arbitrarily chosen
                                    my $ensembl_id_obj = Genome::ExternalGeneId->get(data_directory => $self->annotation_build->data_directory."/annotation_data", reference_build_id => $self->annotation_build->reference_sequence_id, id_type => "ensembl", gene_id => $other_obj->gene_id);
                                    if ($ensembl_id_obj) {
                                        $ensembl_id = $ensembl_id_obj->id_value;
                                        $count++;
                                        $str .= $ensembl_id;
                                    }
                                    else {
                                        $err_out->print("Could not translate $gene_name into ensembl id\n");
                                        $ensembl_id = "-";
                                    }
                                }
                                if ($count > 1) {
                                    $self->warning_message("More than one ensembl_id found for gene name $gene_name: $str.  Arbitrarily chose $ensembl_id");
                                }
                            }
                            else {
                                $err_out->print("Could not find gene with default name $gene_name\n");
                                $ensembl_id = "-";
                            }
                        }
                        push @{$annotations{$gene_name}->{sources}}, "annovar_$column_name";
                        push @{$annotations{$gene_name}->{trv_type}}, "regulatory";
                        if ((not defined $annotations{$gene_name}->{ensembl_id}) or $annotations{$gene_name}->{ensembl_id} eq "-") {
                            $annotations{$gene_name}->{ensembl_id} = $ensembl_id;
                            $err_out->print("Setting ensembl_id to $ensembl_id\n");
                        }
                        else {
                            unless ($annotations{$gene_name}->{ensembl_id} eq $ensembl_id) {
                                $self->warning_message("Ensembl id does not match: ".$annotations{$gene_name}->{ensembl_id}." and $ensembl_id");
                            }
                        }
                    }
                }
            }
            unless (%annotations) {
                $annotations{"-"}->{sources} = ["-"];
                $annotations{"-"}->{trv_type} = [$fields[13]];
                $annotations{"-"}->{ensembl_id} = "-";
            }
            foreach my $annot (keys %annotations) {
                $fields[6] = $annot;
                $fields[21] = $annot;
                $fields[22] = join(",", @{$annotations{$annot}->{sources}});
                $fields[13] = join(",", @{$annotations{$annot}->{trv_type}});
                $fields[23] = $annotations{$annot}->{ensembl_id};
                $out->print(join("\t", @fields)."\n");
            }
        }
    }
    $out->close;
    $err_out->print("closing err_out\n");
    $err_out->close;
    return 1;
}

1;

