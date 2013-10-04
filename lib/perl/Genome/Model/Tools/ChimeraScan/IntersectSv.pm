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
    
    my $fusion_reader = Genome::File::BedPe::Reader->new($filtered_bedpe_file);
    my %fusion_transcripts;

    #gather chimera filter info
    my @fusion_headers        = qw(name chrom1 start1 end1 chrom2 start2 end2);
    my @fusion_custom_headers = qw(fusion total_frag spanning_frag);
    my @fusion_all_headers    = (@fusion_headers, @fusion_custom_headers);

    while (my $entry = $fusion_reader->next) {
        my ($id_5p, $id_3p) = map{$entry->{custom}->[$_]}qw(3 4);
        
        for my $id ($id_5p, $id_3p) {
            my @info = split /,/, $id;
            for my $info (@info) {
                my ($transcript_id) = $info =~ /^(\S+)\:/;

                my (%hash, %custom_hash);
                @hash{@fusion_headers} = map{$entry->{$_}}@fusion_headers;
                @custom_hash{@fusion_custom_headers} = map{$entry->{custom}->[$_]}qw(0 5 6);

                $fusion_transcripts{$transcript_id} = {%hash, %custom_hash};
            }
        }
    }

    #gather sv out info
    my $sv_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $sv_output_file,
        separator => "\t",
    );

    my %breakpoint_list;
    my @sv_headers = qw(chrA bpA chrB bpB event geneA ensemblIdA geneB ensemblIdB);

    while (my $line = $sv_reader->next) {
        my ($chrA, $chrB) = ($line->{chrA}, $line->{chrB});
        my %hash;
        @hash{@sv_headers} = map{$line->{$_}}@sv_headers;
        push @{$breakpoint_list{$chrA}}, \%hash;
        push @{$breakpoint_list{$chrB}}, \%hash unless $chrA eq $chrB;
    }

    my @all_headers = (@sv_headers, @fusion_all_headers);

    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        headers   => \@all_headers,
        separator => "\t",
        output    => $self->output_file,
    );

    my $build = $self->annotation_build
        or die "failed to get annotation build ".$self->annotation_build_id;

    my %keys;

    for my $chr (nsort keys %breakpoint_list) {
        my $chr_breakpoint_list = $breakpoint_list{$chr};
        my $transcript_iterator = $build->transcript_iterator(chrom_name => $chr);
        die "transcript iterator not defined for chr $chr" unless $transcript_iterator;

        while (my $transcript = $transcript_iterator->next) {
            my $transcript_name = $transcript->transcript_name; 
            my $fusion_item     = $fusion_transcripts{$transcript_name};
            next unless $fusion_item;

            for my $item (@$chr_breakpoint_list) {
                my $key = join '--', map{$item->{$_}}qw(chrA bpA chrB bpB event);
                next if $keys{$key};

                my $flag = 0;
                if ($item->{chrA} eq $chr) {
                    $flag++ if $transcript->within_transcript($item->{bpA});
                }
                if ($item->{chrB} eq $chr) {
                    $flag++ if $transcript->within_transcript($item->{bpB});
                }

                if ($flag) {
                    $keys{$key} = 1;
                    my %sv_info;
                    @sv_info{@sv_headers} = map{$item->{$_}}@sv_headers;
                    my %fusion_info;
                    @fusion_info{@fusion_all_headers} = map{$fusion_item->{$_}}@fusion_all_headers;
                    my %content = (%sv_info, %fusion_info);
                    $writer->write_one(\%content);
                }
            }
        }
    }

    return 1;
}

1;

