package Genome::Model::Tools::Annotate::DescribeAnnotationRegions;

use strict;
use warnings;

use Genome; 
use Genome::Info::UCSCConservation;
use List::MoreUtils qw/ uniq /;

class Genome::Model::Tools::Annotate::DescribeAnnotationRegions{
    is => 'Genome::Model::Tools::Annotate',
       has => [
           variant_file => {
               is => 'Text',   
               doc => "File of dummy variants. Tab separated columns: chromosome_name start stop reference variant",
           },
           output_file => {
               is => 'Text',
               doc => "Store annotation in the specified file. Defaults to STDOUT if no file is supplied.",
               default => "STDOUT",
           },
           reference_transcripts => {
               is => 'String',
               doc => 'provide name/version number of the reference transcripts set you would like to use ("NCBI-human.combined-annotation/0").  Leaving off the version number will grab the latest version for the transcript set, and leaving off this option and build_id will default to using the latest combined annotation transcript set. Use this or --build-id to specify a non-default annoatation db (not both). See full help output for a list of available reference transcripts. Defaults to NCBI-human.combined-annotation/54_36p_v2',
               default => 'NCBI-human.combined-annotation/54_36p_v2',
           },
       ],
       has_optional => [
           flank_range => {
               is => 'Integer', 
               is_optional => 1,
               default => 50000,
               doc => 'Range to look around for flanking regions of transcripts',
           },
           annotation_build_version => {
                is => "Text",
                is_optional => 1,
           }
       ]
};

#NOTE: This is made of bits of Genome/Transcript/VariantAnnotator.pm and
#Genome/Model/Tools/Annotate/TranscriptVariants.pm that have been stripped
#down and pasted together.  It basically treats the "variant" file as a 
#list of regions to run the annotator on, but with all variant specific (ie
#change specific) bits stripped off.  If either of these parent modules
#diverge in these areas of functionality, this is likely to break.  
#YOU HAVE BEEN WARNED!!!  

sub execute {
    my $self = shift;

    my $variant_file = $self->variant_file;

    # preserve additional columns from input if desired 
    my @columns = ($self->variant_output_attributes);
    my $variant_svr = Genome::Utility::IO::SeparatedValueReader->create(
        input => $variant_file,
        headers => \@columns,
        separator => "\t",
        is_regex => 1,
        ignore_extra_columns => 1,
    );
    unless ($variant_svr) {
        $self->error_message("error opening file $variant_file");
        return;
    }

    # establish the output handle for the transcript variants
    my $output_fh;
    my $output_file = $self->output_file;
    if ($self->output_file =~ /STDOUT/i) {
        $output_fh = 'STDOUT';
    }
    else {
        $output_fh = $self->_create_file($output_file);
    }
    $self->_transcript_report_fh($output_fh);

    my $ref = $self->reference_transcripts;
    $ref = "NCBI-human.combined-annotation/54_36p_v2" unless defined $ref;
    my ($name) = split(/\//, $ref); # For now, version is ignored since only v2 is usable
                                    # This will need to be changed when other versions are available

    my $model = Genome::Model->get(name => $name);
    unless ($model){
        $self->error_message("couldn't get reference transcripts set for $name");
        return;
    }

    my $version = "54_37g_v2" if $name =~ /mouse/i;
    $version = "54_36p_v2" if $name =~ /human/i;
    unless (defined $version) {
        $self->error_message("Couldn't determine latest version for model $name");
        return;
    }

    my $build = $model->build_by_version($version);
    unless ($build){
        $self->error_message("couldn't get build from reference transcripts set $name");
        return;
    }

    my $full_version = $build->version; 
    my ($version_number) = $full_version =~ /^\d+_(\d+)[a-z]/;
    my %ucsc_versions = Genome::Info::UCSCConservation->ucsc_conservation_directories;
    my $ucsc_dir = $ucsc_versions{$version_number};
    $self->annotation_build_version($full_version);

    $output_fh->print(join("\t", $self->transcript_report_headers), "\n" );

    while ( my $variant = $variant_svr->next) {
        my $transcript_iterator;
            $transcript_iterator = $build->transcript_iterator(chrom_name => $variant->{chromosome_name});
        unless ($transcript_iterator){
            $self->error_message("Couldn't get transcript_iterator for chromosome " . $variant->{chromosome_name});
            die;
        }

        my $transcript_window =  Genome::Utility::Window::Transcript->create (
                iterator => $transcript_iterator, 
                range => $self->flank_range
        );
        unless ($transcript_window){
            $self->error_message("Couldn't create a transcript window from iterator for chromosome " . $variant->{chromosome_name});
            die;
        }

        my @transcripts_to_annotate = $transcript_window->scroll($variant->{start});
        my @annotations;
        for my $transcript (@transcripts_to_annotate){
            my %annotation = $self->_transcript_annotation($transcript, $ucsc_dir, %$variant) or next;
            push @annotations, \%annotation;
        }
        
        $self->_print_annotation($variant, \@annotations);

    }

    $output_fh->close unless $output_fh eq 'STDOUT';
}

sub _transcript_report_fh {
    my ($self, $fh) = @_;
    $self->{_transcript_fh} = $fh if $fh;
    return $self->{_transcript_fh};
}

sub transcript_report_headers {
    my $self = shift;
    return ($self->variant_output_attributes, $self->transcript_attributes);
}

sub _transcript_annotation{
    my ($self, $transcript, $ucsc_dir, %variant)  = @_;
    my $conservation = $self->_ucsc_conservation_score(\%variant, $transcript);
     my $all_protein_domains = $self->_protein_domain($transcript, \%variant);

    return (
        transcript_error => $transcript->transcript_error,
        transcript_name => $transcript->transcript_name, 
        transcript_status => $transcript->transcript_status,
        transcript_source => $transcript->source,
        transcript_species=> $transcript->species,
        transcript_version => $transcript->version,
        strand => $transcript->strand,
        gene_name  => $transcript->gene->name,
        amino_acid_length => $transcript->amino_acid_length,
        ucsc_cons => $conservation,
        all_domains => $all_protein_domains,
    );
}

sub transcript_attributes{
    my $self = shift;
    my @attrs = qw/gene_name	transcript_name	transcript_species	transcript_source	transcript_version	strand	transcript_status	ucsc_cons	all_domains	 transcript_error/;
    return @attrs;
}

sub variant_output_attributes{
    return qw/chromosome_name start   stop/;
}

sub _print_annotation {
    my ($self, $snp, $transcripts) = @_;

    # Basic SNP Info 
    my $snp_info_string = join
    (
        "\t", 
        map { $snp->{$_} } ($self->variant_output_attributes),
    );

    # If we have no transcripts, print the original variant with dashes for annotation info
    unless( @$transcripts ) {
        $self->_transcript_report_fh->print
        (
            join
            (
                "\t",                   
                $snp_info_string,
                map({ '-' } $self->transcript_attributes),
            ), 
            "\n",
        );
        return 1;
    }

    # Otherwise, print an annotation line for each transcript we have
    for my $transcript ( @$transcripts )
    {
        $self->_transcript_report_fh->print
        (
            join
            (
                "\t",                   
                $snp_info_string,
                map({ $transcript->{$_} ? $transcript->{$_} : '-' } $self->transcript_attributes),
            ), 
            "\n",
        );
    }
    return 1;
}

sub _ucsc_conservation_score {
    my ($self, $variant, $transcript) = @_;
    return 'NULL' if $variant->{chromosome_name} =~ /^[MN]T/;

    my $range = [ $variant->{start}..$variant->{stop} ] ;
    my $conservation_score_lookup = Genome::Model::Tools::Annotate::LookupConservationScore->execute(chromosome => $variant->{chromosome_name}, coordinates => $range, species => $transcript->species, version => $self->annotation_build_version);
    my $ref = $conservation_score_lookup->conservation_scores_results;
    my @ret;
    foreach my $item (@$ref)
    {
        push(@ret,sprintf("%.3f",$item->[1]));
    }
    return join(":",@ret);

}

sub _protein_domain {
    my ($self, $transcript, $variant) = @_;
    return 'NULL' unless defined $transcript and defined $variant;

    my @all_domains = Genome::InterproResult->get(
            transcript_name => $transcript->transcript_name,
            data_directory => $transcript->data_directory,
            chrom_name => $variant->{chromosome_name},
            );
    return 'NULL' unless @all_domains;

    my @all_domain_names;
    for my $domain (@all_domains) {
        push @all_domain_names, $domain->{name};
    }

    return join(",", uniq @all_domain_names);
}

1;
