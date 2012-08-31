package Genome::Model::Tools::Annotate::TranscriptVariantsStream;

use strict;
use warnings;

use Genome; 

use Command;
use Data::Dumper;
use IO::File;
use Genome::Info::IUB;

class Genome::Model::Tools::Annotate::TranscriptVariantsStream{
    is => 'Genome::Model::Tools::Annotate',
    has_optional => [
        # Transcript Params
        annotation_filter => {
            is => 'String',
            is_optional => 1,
            default => 'gene',
            doc => 'The type of filtering to use on the annotation results. There are currently 3 valid options:
                    "none" -- This returns all possible transcript annotations for a variant. All transcript status are allowed including "unknown" status.
                    "gene" -- This returns the top transcript annotation per gene. This is the default behavior.
                    "top" -- This returns the top priority annotation for all genes. One variant in, one annotation out.',
        },
        flank_range => {
            is => 'Integer', 
            is_optional => 1,
            default => 50000,
            doc => 'Range to look around for flanking regions of transcripts',
        },
        reference_transcripts => {
            is => 'String',
            is_optional => 1, 
            doc => 'provide name/version number of the reference transcripts set you would like to use ("NCBI-human.combined-annotation/54_36p_v2").  Leaving off the version number will grab the latest version for the transcript set, and leaving off this option and build_id will default to using the latest combined annotation transcript set. Use this or --build-id to specify a non-default annoatation db (not both)'
        },
        build_id =>{
            is => "Number",
            is_optional => 1,
            doc => 'build id for the imported annotation model to grab transcripts to annotate from.  Use this or --reference-transcripts to specify a non-default annotation db (not both)',
        },
        build => {
            is => "Genome::Model::Build",
            id_by => 'build_id',
            is_optional => 1, 
        },
        extra_details => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'enabling this flag produces an additional four columns: flank_annotation_distance_to_transcript, intron_annotation_substructure_ordinal, intron_annotation_substructure_size, and intron_annotation_substructure_position',
        },
    ], 
};

############################################################

sub annotate { 
    my $self = shift;
    my ($chromosome, $start, $stop, $reference, $variant_bases) = @_;
    
    my $variant;
    $variant->{chromosome_name} = $chromosome;
    $variant->{start} = $start;
    $variant->{stop} = $stop;
    $variant->{reference} = $reference;
    $variant->{variant} = $variant_bases;
    $variant->{type} = Genome::Model::Tools::Annotate->infer_variant_type($variant);

    if ($self->build and $self->reference_transcripts){
        $self->error_message("Please provide a build id OR a reference transcript set name, not both");
        return;
    }
    
    #if no build is provided, use the v0 of our generic NCBI-human-36 imported annotation model
    unless ($self->build){
        if ($self->reference_transcripts){
            my ($name, $version) = split(/\//, $self->reference_transcripts);
            my $model = Genome::Model->get(name => $name);
            unless ($model){
                $self->error_message("couldn't get reference transcripts set for $name");
                return;
            }
            if (defined($version)){
                my $build = $model->build_by_version($version);
                unless ($build){
                    $self->error_message("couldn't get version $version from reference transcripts set $name");
                    return;
                }
                $self->build($build);
            }else{ 
                my $build = $model->last_complete_build;  #TODO latest by version
                unless ($build){
                    $self->error_message("couldn't get last complete build from reference transcripts set $name");
                    return;
                }
                $self->build($build);
            }
        }else{
            my $model = Genome::Model->get(name => 'NCBI-human.combined-annotation');
            my $build = $model->build_by_version(0);

            unless ($build){
                $self->error_message("couldn't get build v0 from 'NCBI-human.combined-annotation'");
                return;
            }
            $self->build($build);
        }
    }
            
    my $transcript_iterator = $self->build->transcript_iterator(chrom_name => $variant->{chromosome_name});
    unless ($transcript_iterator){
        die "Transcript iterator could not be created for chromosome " . $variant->{chromosome_name};
    }

    my $transcript_window =  Genome::Utility::Window::Transcript->create (
        iterator => $transcript_iterator, 
        range => $self->flank_range
    );
    unless ($transcript_window){
        $self->error_message("Couldn't create a transcript window from iterator for chromosome " . $variant->{chromosome_name});
        die;
    }

    my $annotator = Genome::Transcript::VariantAnnotator->create(
        transcript_window => $transcript_window,
    );
    unless ($annotator){
        $self->error_message("Couldn't create annotator for chromosome " . $variant->{chromosome_name});
        die;
    }

    # If we have an IUB code, annotate once per base... doesnt apply to things that arent snps
    # TODO... unduplicate this code
    my $annotation_filter = lc $self->annotation_filter;
    if ($variant->{type} eq 'SNP') {
        my @variant_alleles = Genome::Info::IUB->variant_alleles_for_iub($variant->{reference}, $variant->{variant});
        for my $variant_allele (@variant_alleles) {
            # annotate variant with this allele
            $variant->{variant} = $variant_allele;

            # get the data and output it
            my $annotation_method;
            if ($annotation_filter eq "gene") {
                # Top annotation per gene
                $annotation_method = 'prioritized_transcripts';
            } elsif ($annotation_filter eq "top") {
                # Top annotation between all genes
                $annotation_method = 'prioritized_transcript';
            } elsif ($annotation_filter eq "none") {
                # All transcripts, no filter
                $annotation_method = 'transcripts';
            } else {
                $self->error_message("Unknown annotation_filter value: " . $annotation_filter);
                return;
            }

            my @transcripts = $annotator->$annotation_method(%$variant);
            return @transcripts;
        }
    } else {
        # get the data and output it
        my @transcripts;
        if ($annotation_filter eq "gene") {
            # Top annotation per gene
            @transcripts = $annotator->prioritized_transcripts(%$variant);
        } elsif ($annotation_filter eq "top") {
            # Top annotation between all genes
            @transcripts = $annotator->prioritized_transcript(%$variant);
        } elsif ($annotation_filter eq "none") {
            # All transcripts, no filter
            @transcripts = $annotator->transcripts(%$variant);
        } else {
            $self->error_message("Unknown annotation_filter value: " . $annotation_filter);
            return;
        }
        
        return @transcripts;
    }

    # Shouldnt get here
    return 0;
}

sub transcript_attributes{
    my $self = shift;
    my @attrs = $self->SUPER::transcript_attributes;
    if ($self->extra_details){
        push @attrs, (qw/ flank_annotation_distance_to_transcript intron_annotation_substructure_ordinal intron_annotation_substructure_size intron_annotation_substructure_position/);
    }
    return @attrs;
}

1;
