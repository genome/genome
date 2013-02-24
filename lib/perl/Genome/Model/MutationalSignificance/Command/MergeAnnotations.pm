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
    ],
};

sub execute {
    my $self = shift;

    my $tgi = Genome::Sys->open_file_for_reading($self->tgi_anno_file);
    my $reg_db = Genome::Sys->open_file_for_reading($self->regulome_db_file);
    my $annovar = Genome::Sys->open_file_for_reading($self->annovar_file);

    my $out = Genome::Sys->open_file_for_writing($self->output_file);

    my %variants;
    #read in the tgi file first
    while (my $line = <$tgi>) {
        chomp $line;
        my @fields = split /\t/, $line;
        my $chr = $fields[0];
        my $start = $fields[1];
        if (defined $variants{$chr}{$start}->{tgi}) {
            $self->error_message("Duplicate variant: $chr $start");
            return;
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
                my @hits = split /,\s/, $fields[3];
                my @eQTLs;
                foreach my $hit (@hits) {
                    if ($hit =~ /eQTL/) {
                        my @parts = split /|/, $hit;
                        push @eQTLs, $parts[1];
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
                $self->error_message("No regulome db score for variant $chr:$start");
                return;
            }

            my $line = $var{tgi};
            my @fields = split /\t/, $line;
            #build up a list of gene annotations
            my %annotations;
            unless ($fields[21] eq "-") {
                $annotations{$fields[21]}->{sources} = [$fields[22]];
                $annotations{$fields[21]}->{trv_type} = [$fields[13]];
            }
            if ($var{regdb_score} =~ /1/) {
                foreach my $eQTL (@{$var{eQTLs}}){
                    push @{$annotations{$eQTL}->{sources}}, "regulomeDB_eQTL";
                    push @{$annotations{$eQTL}->{trv_type}}, "regulatory";
                }
            }

            if ($var{regdb_score} =~ /1/ || $var{regdb_score} =~ /2/) {
                foreach my $column_name ($self->annovar_columns_to_check) {
                    my $annot = $var{"annovar_$column_name"};
                    unless ($annot) {
                        $self->error_message("No annovar_$column_name for variant $chr:$start");
                        return;
                    }
                    next if ($annot =~ "-");
                    $annot =~ s/^Name=//;
                    my @gene_names = split /,/, $annot;
                    foreach my $gene_name (@gene_names) {
                        push @{$annotations{$gene_name}->{sources}}, "annovar_$column_name";
                        push @{$annotations{$gene_name}->{trv_type}}, "regulatory";
                    }
                }
            }
            unless (%annotations) {
                $annotations{"-"}->{sources} = ["-"];
                $annotations{"-"}->{trv_type} = [$fields[13]];
            }
            foreach my $annot (keys %annotations) {
                $fields[6] = $annot;
                $fields[21] = $annot;
                $fields[22] = join(",", @{$annotations{$annot}->{sources}});
                $fields[13] = join(",", @{$annotations{$annot}->{trv_type}});
                $out->print(join("\t", @fields)."\n");
            }
        }
    }
    $out->close;
    return 1;
}

1;

