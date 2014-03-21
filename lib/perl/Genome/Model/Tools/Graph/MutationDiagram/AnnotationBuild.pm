package Genome::Model::Tools::Graph::MutationDiagram::AnnotationBuild;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Graph::MutationDiagram::AnnotationBuild {
    is => "Genome::Model::Tools::Graph::MutationDiagram::DomainProvider",
    has => [
        build => {
            is => "Genome::Model::Build::ImportedAnnotation",
        },
    ],
};

sub get_domains {
    my ($self, $transcript_name) = @_;
    my $build = $self->build;
    my @features;
    my $transcript;
    for my $data_directory ($build->determine_data_directory){
        my $t = Genome::Transcript->get(
            data_directory => $data_directory,
            transcript_name => $transcript_name,
            reference_build_id => $build->reference_sequence_id,
            );
        $transcript = $t;
        next unless $t;
        push(@features, Genome::InterproResult->get(
            data_directory => $data_directory,
            transcript_name => $transcript_name,
            chrom_name => $t->chrom_name
            ));
    }
    if (!defined $transcript) {
        warn "No transcript found for $transcript_name";
        return;
    }
$DB::single=1;
    my @domains = $self->SUPER::get_domains;
    for my $feature (@features) {
        my ($source, @domain_name_parts) = split(/_/, $feature->name);
        # Some domain names are underbar delimited, but sources aren't.
        # Reassemble the damn domain name if necessary
        my $domain_name;
        if (scalar (@domain_name_parts) > 1){
            $domain_name = join("_", @domain_name_parts);
        }else{
            $domain_name = pop @domain_name_parts;
        }
        push @domains, {
            name => $domain_name,
            source => $source,
            start => $feature->start,
            end => $feature->stop
        };
    }
    return @domains;
}

sub get_amino_acid_length {
   my $self = shift;
   my $transcript_name = shift;
   for my $data_directory ($self->build->determine_data_directory){
       my $t = Genome::Transcript->get(
           data_directory => $data_directory,
           transcript_name => $transcript_name,
           reference_build_id => $self->build->reference_sequence_id,
       );
       if ($t) {
         return $t->amino_acid_length;
       }
   }
   die $self->error_message("Could not find transcript $transcript_name in annotation build ".$self->build->id);
}

1;

