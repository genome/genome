package Genome::Model::Tools::Sv::ConvertTranscriptToAnnotationBam;

use strict;
use warnings;

use Genome; 

use Command;
use Data::Dumper;
use IO::File;

class Genome::Model::Tools::Sv::ConvertTranscriptToAnnotationBam {
    is => 'Genome::Model::Tools::Annotate',
    has => [ 
    ], 
    has_optional => [
        reference_transcripts => {
            is => 'String',
            is_optional => 0,
            doc => 'provide name/version number of the reference transcripts set you would like to use ("NCBI-human.combined-annotation/0").  Leaving off the version number will grab the latest version for the transcript set, and leaving off this option and build_id will default to using the latest combined annotation transcript set. Use this or --build-id to specify a non-default annoatation db (not both). See full help output for a list of available reference transcripts.'
        },
        build => {
            is => "Genome::Model::Build",
            id_by => 'build_id',
            is_optional => 1, 
        },
        file_name => {
            is => 'String',
            is_optional => 0,
            doc => "output file name. .bam will be added automatically",
        },
        samtools_version => {
            is => 'String',
            is_optional => 0,
            default => 'r982',
            doc => "Version of samtools to use",
        },
        reference_sequence_index => {
            is => "String",
            is_optional => 0,
            doc => 'samtools fai index of the reference sequence the annotation database maps to',
        },
    ], 
};


############################################################

sub help_synopsis { 
}

sub help_detail {
    my $self = shift;

    #Generate the currently available annotation models on the fly
    my @currently_available_models = Genome::Model->get(type_name => "imported annotation");
    my $currently_available_builds; 
    foreach my $model (@currently_available_models) {
        next unless $model;
        foreach my $build ($model->builds) {
            if($build and $build->status eq 'Succeeded' and $build->name =~ /combined-annotation/) {  #probably implicit in the loops, but in case we get undefs in our list
                 $currently_available_builds .= "\t" . $model->name . "/" . $build->version . "\n" if $build->version !~ /old/ and $model->name and $build->version;
            }
        }
    }

    return <<EOS 
This dumps a TGI annotator database into a bam file for use with the pairoscope program for visualizing SV.    

COMMONLY USED REFERENCE TRANSCRIPTS WITH VERSIONS
$currently_available_builds
EOS
}

############################################################

sub execute { 
    my $self = shift;

    # establish the output handle for the transcript variants
    #if no build is provided, use the v0 of our generic NCBI-human-36 imported annotation model
    my $ref = $self->reference_transcripts;
    my ($name, $version) = split(/\//, $ref); # For now, version is ignored since only v2 is usable
    # This will need to be changed when other versions are available

    my $model = Genome::Model->get(name => $name);
    unless ($model){
        $self->error_message("couldn't get reference transcripts set for $name");
        return;
    }


    my $build = $model->build_by_version($version);
    unless ($build){
        $self->error_message("couldn't get build from reference transcripts set $name");
        return;
    }
    $self->build($build);


    #grab a temp file to write the .sam to
    my ($tfh,$temp_file) = Genome::Sys->create_temp_file();


    # annotate all of the input variants
    # it would probably be better to derive these names from the reference fasta then to have it hard-coded
    for my $chrom_name (1..22,"X","Y","MT") {
    #for my $chrom_name (7) {
    #    $DB::single=1;
        my $transcript_iterator = $self->build->transcript_iterator(chrom_name => $chrom_name);
        unless ($transcript_iterator){
            $self->error_message("Couldn't get transcript_iterator for chromosome $chrom_name!");
            next;
        }
        while(my $transcript = $transcript_iterator->next) {
            #       $DB::single = 1 if $transcript->transcript_name eq 'NM_006854';
            #my $structure_ordinal = 0;
            my @ordered_substructures = grep {$_->structure_type eq 'cds_exon' || $_->structure_type eq 'utr_exon'} $transcript->ordered_sub_structures;
            if($transcript->strand eq '-1') {
                @ordered_substructures = reverse @ordered_substructures;
            }
            my $number_structures = scalar(@ordered_substructures);
            my ($max_ordinal) = sort { $b <=> $a} map { $_->ordinal } @ordered_substructures;
            my $status = $transcript->transcript_status;
            my $source = $transcript->source;
            my $gene = $transcript->gene->name;
            my $length = 0;
            foreach my $substructure (@ordered_substructures) {
                $length += $substructure->length;
            }
            my $i = 0;
            while($i < $number_structures) {
                my $substructure = $ordered_substructures[$i];
                my $structure_size = $substructure->length;
                my $next_subtructure;
                my $j = $i+1;
                while(my $temp_substructure = $ordered_substructures[($j) % $number_structures]) {
                    if($j < $number_structures && $temp_substructure->ordinal == $substructure->ordinal) {
                        #these are in the same exon and need to be merged
                        $structure_size += $temp_substructure->length;
                        $j++;
                    }
                    else {
                        $next_subtructure = $temp_substructure;
                        $i = $j;
                        last;
                    }
                }

                my $next_pos = $next_subtructure->structure_start;
                my $structure_ordinal = $substructure->ordinal - 1;
                printf $tfh "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$transcript->transcript_name,$transcript->strand eq '-1' ? 0x0010 : 0x0000, $chrom_name, $substructure->structure_start, 255,$structure_size . "M", '*',0,0,'=' x $structure_size,'*',"NH:i:$max_ordinal\tIH:i:$max_ordinal\tHI:i:$structure_ordinal\tCC:Z:$chrom_name\tCP:i:$next_pos\tYG:Z:$gene\tYS:Z:$source\tYT:Z:$status\tYL:i:$length\n"; #this is just required fields
                #$structure_ordinal++;
            }
        }
    }
    $tfh->close;

    #convert to bam with samtools
    my $samtools = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
    my $reference_index = $self->reference_sequence_index;
    my $file_name = $self->file_name;
    my $cmd = "$samtools view -b -t $reference_index $temp_file | $samtools sort -m 20000000 - $file_name"; 
    Genome::Sys->shellcmd(cmd => $cmd);
    my $index_cmd = "$samtools index $file_name.bam";
    Genome::Sys->shellcmd(cmd => $index_cmd);

    return 1;
}
1;
