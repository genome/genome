package Genome::Model::Tools::RefCov::BuildReferenceFromRoi;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RefCov::BuildReferenceFromRoi {
    is => ['Command'],
    has_input => [
        input_fasta_file => {
            is => 'Text',
            doc => 'The input reference FASTA file.',
        },
        bed_file => {
            is => 'Text',
            doc => 'The input BED file of regions-of-interest(ROI) to create new reference sequences for.',
            is_output => 1,
        },
        output_fasta_file => {
            is => 'Text',
            doc => 'The output reference FASTA file',
        },
    ],
    has_optional_input => [
        stitch_level => {
            default_value => 'gene',
            valid_values => ['gene','transcript','none'],
            doc => 'The gene level stitching requires a squashed version of the BED file reduced to a single longest representation.  The transcript level stitching requires the complete BED file without squashing.',
        },
        delimiter => {
            is => 'Text',
            doc => 'The character that delimits the ROI names like GENE, TRANSCRIPT, TYPE, ORDINAL and STRAND',
            default_value => ':',
            valid_values => [':','_','.'],
        },
        coordinates => {
            is => 'Boolean',
            doc => 'Report a comma delimitted list of chromosome and coordinates for each stiched block of sequence in the FASTA description',
            default_value => 0,
        },
        reverse_complement => {
            is => 'Boolean',
            doc => 'Reverse complement the negative strand ROI',
            default_value => 0,
        },
        junctions_bed_file => {
            is => 'Text',
            doc => 'An output BED format file with coordinates of donor-acceptor sites',
        },
        junction_padding => {
            is => 'Text',
            default_value => 2,
            doc => 'The number of base pair to pad on either side of splice junction.',
        },
    ],
};

sub help_brief {
    "Given regions of interest belonging to a larger sequence this command will build a smaller reference sequence limited to those regions (ex. Transcriptome).",
}

sub help_detail {
return <<EOS
A FASTA entry per ROI is generated based on the coordinates of the features
in the reference sequence. The BED file must contain gene and transcript
names as defined by the gtf-to-bed converter.  Those names allow this command
to build full length gene or transcript sequences including a BED file of
the splice sites.

This command is intented to provide a means for producing a transcriptome file.
Other tools or commands like BEDTools or the Cufflinks command 'gffread' also
provide this function.
EOS
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }

    require Bio::DB::Sam;
    return $self;
}

sub execute {
    my $self = shift;

    # Load FASTA index
    my $input_fasta_file = $self->input_fasta_file;
    unless (-f $input_fasta_file .'.fai') {
        die('Failed to find fai index for fasta file '. $input_fasta_file);
    }
    my $fai = Bio::DB::Sam::Fai->load($input_fasta_file);
    unless ($fai) {
        die('Failed to load fai index for fasta file '. $input_fasta_file);
    }
    
    # Load ROI Regions
    my $sorted_bed_file = Genome::Sys->create_temp_file_path;
    unless (Genome::Model::Tools::BedTools::Sort->execute(
        use_version => '2.14.3',
        input_file => $self->bed_file,
        output_file => $sorted_bed_file,
    )) {
        die('Failed to sort BED file'. $self->bed_file);
    }
    my $regions = Genome::Model::Tools::RefCov::ROI::Bed->create(
        file => $sorted_bed_file,
    );
    unless ($regions) {
        die('Failed to load BED regions-of-interest file: '. $self->bed_file);
    }

    # Open files for output
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_fasta_file);

    my $junctions_fh;
    if ($self->junctions_bed_file) {
        $junctions_fh = Genome::Sys->open_file_for_writing($self->junctions_bed_file);
    }

    my %sequences;
    # Iterate over all ROI
    while (my $region = $regions->next_region) {
        my $name = $region->{name};
        my $id = $region->{id};

        # Get DNA sequence string
        my $dna_string = $fai->fetch($id);
        if ( !defined($dna_string) || ( length($dna_string) != $region->{length} )  ) {
            die('Failed to fetch the proper length ('. $region->{length} .') dna.  Fetched '. length($dna_string) .' bases for region: '. $region->{id});
        }

        # Parse ROI names for annotation info
        my ($gene,$transcript,$type,$ordinal,$strand) = split($self->delimiter,$name);
        my $sequence_name;
        if ($self->stitch_level eq 'gene') {
            $sequence_name = $gene;
        } elsif ($self->stitch_level eq 'transcript') {
            $sequence_name = $transcript;
        } elsif ($self->stitch_level eq 'none') {
            $sequence_name = $name;
        }

        # Create unique key based on name chrom and strand
        my $key = $sequence_name .'|'. $region->{chrom}.'|'.$strand;

        if ($sequences{$key}) {
            # Append ROI info to existing data
            $sequences{$key}{dna_string} .= $dna_string;
            unless ($sequences{$key}{chrom} eq $region->{chrom}) {
                warn('Invalid chromosome for region!');
                warn(Data::Dumper::Dumper($region));
                warn(Data::Dumper::Dumper($sequences{$key}));
            }
            push @{$sequences{$key}{exons}}, $region->{length};
            push @{$sequences{$key}{coords}}, $region->{start} .'-'.$region->{end};
            unless ($sequences{$key}{strand} eq $strand) {
                warn('Invalid strand for region!');
                warn(Data::Dumper::Dumper($region));
                warn(Data::Dumper::Dumper($sequences{$key}));
            }
        } else {
            # Create new ROI data
            $sequences{$key}{name} = $sequence_name;
            $sequences{$key}{strand} = $strand;
            $sequences{$key}{dna_string} = $dna_string;
            $sequences{$key}{chrom} = $region->{chrom};
            push @{$sequences{$key}{coords}}, $region->{start} .'-'. $region->{end};
            push @{$sequences{$key}{exons}}, $region->{length};
        }
    }
    my $junction_padding = $self->junction_padding;
    for my $key (keys %sequences) {
        # Print the header portion of the FASTA entry
        print $output_fh '>'. $key;
        # Print the genome coordinates of exon positions
        if ($self->coordinates) {
            if ( $self->reverse_complement && (($sequences{$key}{strand} eq '-') || ($sequences{$key}{strand} eq 'rev')) ) {
                my $dna_string = $sequences{$key}{dna_string};
                $sequences{$key}{dna_string} = &reverse_complement_string($dna_string);
                $sequences{$key}{exons} = reverse(@{$sequences{$key}{exons}});
                $sequences{$key}{coords} = reverse(@{$sequences{$key}{coords}});
            }
            my $desc .= 'COORDS:'. join(',',@{$sequences{$key}{coords}});
            print $output_fh ' '. $desc;
        }
        print $output_fh "\n";

        # Print the junction coordinates to a separate BED file relative to new FASTA reference position
        if ($junctions_fh) {
            my $sum = 0;
            my $count = 1;
            for my $exon (@{$sequences{$key}{exons}}) {
                if ($sum) {
                    print $junctions_fh $key ."\t".
                        ($sum - $junction_padding) ."\t".
                            ($sum + $junction_padding) ."\t".
                                $key .'|'. $count++ ."\t.\t".
                                    $sequences{$key}{strand} ."\n";
                }
                $sum += $exon;
            }
        }

        # Print the sequence portion of the FASTA entry
        my $dna_length = length( $sequences{$key}{dna_string} );
        my $pos    = 0;
        my $columns = 60;
        my $lines  = int( $dna_length / $columns) + 1;
        for (my $i = 1; $i <= $lines; $i++) {
            my $line = substr( $sequences{$key}{dna_string}, $pos, $columns );
            if ($line) {
                print $output_fh $line ."\n";
            }
            $pos = $pos + $columns;
        }
    }
    # Clean up filehandles
    if ($junctions_fh) {
        $junctions_fh->close;
    }
    $output_fh->close;
    return 1;
}

sub reverse_complement_string {
    my $dna_string = shift;

    my $rev_dna_string = reverse($dna_string);
    $rev_dna_string =~ tr/ACGTacgt/TGCAtgca/;
    return $rev_dna_string;
}


1;
