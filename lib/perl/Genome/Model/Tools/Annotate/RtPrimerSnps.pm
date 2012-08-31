package Genome::Model::Tools::Annotate::RtPrimerSnps;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Bio::DB::Fasta;


class Genome::Model::Tools::Annotate::RtPrimerSnps {
    is => 'Command',
    has => [
    'snp_file' => {
        type => 'String',
        is_optional => 0,
        doc => 'maq cns2snp output or, minimally, tab separated chr and pos',
    },
    'transcript' => {
        type => 'String',
        is_optional => 0,
        doc => 'transcript name',
    },
    'primer_sequence' => {
        type => 'String',
        is_optional => 0,
        doc => 'sequence of primer to search for snps within',
    },
    reference_transcripts => {
        is => 'String',
        is_optional => 1, 
        doc => 'provide name/version number of the reference transcripts set you would like to use ("NCBI-human.combined-annotation/0").  Leaving off the version number will grab the latest version for the transcript set, and leaving off this option and build_id will default to using the latest combined annotation transcript set. Use this or --build-id to specify a non-default annoatation db (not both)'
    },
    build_id =>{
        is => "Number",
        is_optional => 1,
        doc => 'build id for the imported annotation model to grab transcripts to annotate from.  Use this or --reference-transcripts to specify a non-default annotation db (not both)',
    },
	organism => {
	    type  =>  'String',
	    doc   =>  "provide the organism either mouse or human; default is human",
	    is_optional  => 1,
	    default => 'human',
	},
	version => {
	    type  =>  'String',
	    doc   =>  "provide the imported annotation version; default for human is 54_36p_v2 and for mouse is 54_37g_v2",
	    is_optional  => 1,
	    default => '54_36p_v2',
	},

    ]
};

sub execute {
    my $self = shift;

    unless(-e $self->snp_file) {
        $self->error_message("Input file does not exist");
        return;
    }
    my $build;
    if ($self->reference_transcripts){
        my ($name, $version) = split(/\//, $self->reference_transcripts);
        my $model = Genome::Model->get(name => $name);
        unless ($model){
            $self->error_message("couldn't get reference transcripts set for $name");
            return;
        }
        if (defined($version)){
            $build = $model->build_by_version($version);
            unless ($build){
                $self->error_message("couldn't get version $version from reference transcripts set $name");
                return;
            }
        }else{ 
            $build = $model->last_complete_build;  #TODO latest by version
            unless ($build){
                $self->error_message("couldn't get last complete build from reference transcripts set $name");
                return;
            }
        }
    }else{
        my $model = Genome::Model->get(name => 'NCBI-human.combined-annotation');
        $build = $model->build_by_version(0);
	
        unless ($build){
            $self->error_message("couldn't get build v0 from 'NCBI-human.combined-annotation'");
            return;
        }
    }
    my $reference_build_id = $build->reference_sequence_id;
    ###   my $build_id =$build->build_id;  ###Rather than build id now to get the transcript we need the data directory
    
    ###my $t = Genome::Transcript->get( transcript_name => $self->transcript, build_id => $build_id );
    my $organism = $self->organism;
    my $version = $self->version;
    my $transcript = $self->transcript;
    
    if ($organism eq "mouse") { if ($version eq "54_36p_v2") { $version = "54_37g_v2";}}
    
    my ($ncbi_reference) = $version =~ /\_([\d]+)/;
    my $eianame = "NCBI-" . $organism . ".ensembl";
    my $gianame = "NCBI-" . $organism . ".genbank";
    my $build_source = "$organism build $ncbi_reference version $version";
    
    my $ensembl_build = Genome::Model::ImportedAnnotation->get(name => $eianame)->build_by_version($version);
    my ($ensembl_data_directory) = $ensembl_build->determine_data_directory;
    my $genbank_build = Genome::Model::ImportedAnnotation->get(name => $gianame)->build_by_version($version);
    my ($genbank_data_directory) = $genbank_build->determine_data_directory;
    
    my $t;
    if ($transcript =~/^ENS/){ #ENST for Human ENSMUST
	($t) = Genome::Transcript->get( transcript_name =>$transcript, reference_build_id => $reference_build_id, data_directory => $ensembl_data_directory);
    }else{
	($t) = Genome::Transcript->get( transcript_name =>$transcript, reference_build_id => $reference_build_id, data_directory => $genbank_data_directory)
    }



    if($t) {
        $self->status_message("Found transcript ". $t->transcript_name." on ".substr($t->strand,0,1)." strand for gene ".$t->gene_name."\n\n");
    }
    else {
        $self->error_message("Couldn't find specified transcript: ".$self->transcript);
        return;
    }
    my %transcript_coords;

    my @substructures = grep {$_->structure_type eq 'cds_exon' || $_->structure_type eq 'utr_exon'} $t->ordered_sub_structures;
    my $total_substructures = @substructures;
    my $current_transcript_position = 1;
    my $tseq = "";
    for my $structure (@substructures) {
        if($t->strand == -1) {
            #store starts as appropriate
            $transcript_coords{$current_transcript_position} = $structure->{structure_stop};
            $tseq .= $structure->nucleotide_seq;
        }
        else {
            #must be on + strand
            $transcript_coords{$current_transcript_position} = $structure->{structure_start};
            $tseq .= $structure->nucleotide_seq;
        }
        #reset the offset into the transcript for the next structure start
        #This should be correct for instance with a 1 bp first exon the next start should be 2
        $current_transcript_position += $structure->{structure_stop} - $structure->{structure_start} + 1;
    }

    #at this point the offsets should all be stored
    #now find the primer sequence within the transcript
#    my $tseq = $t->cds_full_nucleotide_sequence;
    my $primer_seq = $self->primer_sequence;
    unless($primer_seq =~ /[ACTG]/i) {
        $self->error_message("Primer Sequence can only contain ACTG. Passed $primer_seq");
        return;
    }
    my ($primer_start, $primer_stop) = (0,0); 
    my $found = 0;
    while($tseq =~ /$primer_seq/gi) {
        unless($found) {
            #looking for uncomplemented sequence in transcript
            $primer_stop = pos($tseq);
            $primer_start = $primer_stop - length($primer_seq) + 1;
            $self->status_message("Found primer from $primer_start to $primer_stop of transcript\n    Primer: $primer_seq\nTranscript: ".substr($tseq, $primer_start-1,$primer_stop - $primer_start + 1) . "\n\n" );
            $found = 1;
        }
        else {
            $self->error_message("Primer found in multiple places");
            return;
        }
    }
    $primer_seq =~ tr/ACTGactg/TGACtgac/;
    $primer_seq = reverse $primer_seq;
    while($tseq =~ /$primer_seq/gi) {
        #looking for uncomplemented sequence in transcript
        unless($found) {
            $primer_start = pos($tseq);
            $primer_stop = $primer_start - length($primer_seq) + 1;
            $self->status_message("Found primer from $primer_start to $primer_stop of transcript\n    Primer: $primer_seq\nTranscript: ".substr($tseq, $primer_stop-1,$primer_start - $primer_stop + 1) . "\n\n" );
            $found = 1;
            #convert to always be transcript orientation
            ($primer_start, $primer_stop) = ($primer_stop, $primer_start);
        }
        else {
            $self->error_message("Primer found in multiple places");
            return;
        }
    }
    unless($found) {
        $self->error_message("Couldn't find primer sequence in transcript");
        return;
    }

    #Do conversion
    #grab any substructures that the primer MIGHT overlap in transcript coordinates
    my @transcript_coords = grep { $_ <= $primer_stop } sort { $a <=> $b } keys %transcript_coords;
    my @genomic_coords = @transcript_coords{@transcript_coords};
    my @primer_genomic_alignments;
    my $current_tcoord = pop @transcript_coords;
    my $current_gcoord = pop @genomic_coords;
    my $current_primer_align_end = $primer_stop;
    while($current_tcoord <= $current_primer_align_end) {
        #convert end to genomic coord
        my $genomic_end;
        my $structure_offset = $current_primer_align_end - $current_tcoord;
        if($t->strand == -1) {
            #offset negative direction from start
            $genomic_end = $current_gcoord - $structure_offset;
        }
        else {
            $genomic_end = $current_gcoord + $structure_offset;
        }

        my $genomic_start = -1;
        if($current_tcoord <= $primer_start) {
            if($t->strand == -1) {
                $genomic_start = $current_gcoord - ($primer_start - $current_tcoord) ;
            }
            else {
                $genomic_start = $primer_start - $current_tcoord + $current_gcoord;
            }
            $current_primer_align_end = -1; #no need to go through the loop again
        }
        else {
            #primer falls across multiple exons
            $current_primer_align_end = $current_tcoord - 1;
            $genomic_start = $current_gcoord;
            $current_tcoord = pop @transcript_coords;
            $current_gcoord = pop @genomic_coords;
        }
        push @primer_genomic_alignments, [$genomic_start, $genomic_end];
    }
    #print converted genomic coordinates
    $self->status_message("Primer genomic coordinates:");
    my $genomic_coords = q{};
    for my $aligned_block (@primer_genomic_alignments) {
        my ($start, $end) = @$aligned_block;
        $genomic_coords .= "$start...$end\t";
    }
    $self->status_message($genomic_coords."\n\n");
    #print retrieved sequence

    my $RefDir;
    if ($organism eq "human"){
	$RefDir = "/gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/";
    } else {
	$RefDir = "/gscmnt/sata147/info/medseq/rmeyer/resources/MouseB37/";
    }
    my $refdb = Bio::DB::Fasta->new($RefDir);

    $self->status_message("Sequence for genomic coordinates:\n");
    my $coords = q{};

    for my $aligned_block (@primer_genomic_alignments) {
        my ($start, $end) = @$aligned_block;
        my $ref_seq =  $refdb->seq($t->chrom_name, $start => $end); 

        $coords .= "$ref_seq\t";
    }
    $self->status_message($coords."\n\n");

    #now search the snp file for SNPs
    $self->status_message("Searching for variants within primer location...");
    my $fh = IO::File->new($self->snp_file,"r");
    unless($fh) {
        $self->error_message("Unable to open snp file");
    }
    my $test = $self->is_within_primer($t->chrom_name, @primer_genomic_alignments);
    while(my $line = $fh->getline) {
        my ($chr, $pos) = split /\t/, $line;
        print $line if $test->($chr,$pos);
    }
    return 1;
}

sub is_within_primer {
    my ($self, $chr, @primer_genomic_alignments) = @_;
    my @ordered_genomic_alignments;
    for my $block (@primer_genomic_alignments) {
        my ($start, $stop) = @$block;
        if($start > $stop) {
            ($start, $stop) = ($stop, $start);
        }
        push @ordered_genomic_alignments, [$start, $stop];
    }
    return sub {
        my ($snp_chr, $snp_pos) = @_;
        if($snp_chr eq $chr) {
            for my $block (@ordered_genomic_alignments) {
                my ($start, $stop) = @$block;
                if($snp_pos >= $start && $snp_pos <= $stop) {
                    return 1;
                }
            }
            return;
        }
        else {
            return;
        }
    };
}

1;

sub help_brief {
    return "This module searches for a primer sequence within a transcript and then reports any snps in the snp file that fall within it";
}

sub help_detail {
    return <<DOC
This module searches for an exact match of the primer sequence within the transcript sequence and uses the resulting transcript to calculate the genomic position of the primer. This is useful for RT primers which may be short and may span exons. 

Once genomic coordinates are determined, the passed snp_file is scanned for SNPs contained within the primer sequence. The snp_file is minimally a tab-separated file of chromosome and position. Upon finding a SNP within the primer the line from the snp_file is printed to stdout. If no snps are found, output to stdout is produced.";
DOC
}
