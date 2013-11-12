package Genome::Model::Tools::FastTier::AddTiers;

use warnings;
use strict;
use Genome;
use IO::File;

class Genome::Model::Tools::FastTier::AddTiers {
    is => 'Command',
    has => [
        input_file => {
            is => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'Input file (First 3 columns are Chr, Start, Stop, where Start is 1-based)',
        },

        input_is_maf => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'Assumes the input file is a MAF, and uses columns 5,6,7 for Chr, Start, Stop',
        },

        output_file => {
            is => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'Output file is equivalent to the input file with tier info added as an additional column',
        },

        tier_file_location =>{
            is => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'annotation directory containing tiering files',
        },

        tiering_version =>{
            is => 'String',
            is_optional => 1,
            is_input => 1,
            doc => 'Use this to override the default (v3) tiering version',
        },
    ]
};

sub help_brief {
    "Add tiering information to an existing file. Use --input-is-maf if input is a MAF file."
}

sub help_detail {
    #get a list of currrently available builds and ouput in the help
    my $listing;
    my @currently_available_models = Genome::Model->get(subclass_name => 'Genome::Model::ImportedAnnotation');
    my @currently_available_builds;
    my @tiering_dirs;
    foreach my $model (@currently_available_models) {
        next unless $model;
        foreach my $build ($model->builds) {
            if($build and $build->status eq 'Succeeded' and $build->name =~ /ensembl/) {  #probably implicit in the loops, but in case we get undefs in our list
                if ($build->version !~ /old/ and $model->name and $build->version && ($build->data_directory !~ /gscarchive/)){
                    push(@currently_available_builds, $model->name . "/" . $build->version . "\n");
                    push(@tiering_dirs, $build->data_directory . "annotation_data/tiering_bed_files_v3/");
                }
            }
        }
    }

    my $i = 0;
    for($i=0;$i<@tiering_dirs;$i++){
        my $name = $currently_available_builds[$i];
        chomp($name);
        my $dir = $tiering_dirs[$i];
        chomp($dir);

        $listing .= $name . "\n    " . $dir . "\n";
    }

    "Add tiering information to an existing file. Use --input-is-maf if input is a MAF file.\n\nCurrently Available Builds and Dirs\n----------------------------------\n$listing\n\n"
}

sub execute {

    my $self = shift;
    my $input_file = $self->input_file;
    my $output_file = $self->output_file;
    my $input_is_maf = $self->input_is_maf;
    my $tier_file_location = $self->tier_file_location;
    my $tiering_version = $self->tiering_version;

    # Create a tmp dir where we'll do the tiering
    my $tempdir = Genome::Sys->create_temp_directory();
    ( -e $tempdir ) or die "Unable to create temporary directory $!";

    # Parse the input file, and create a temporary bed file to use with the fast-tiering tool
    my $bedFh = IO::File->new( "$tempdir/temp.bed", ">" ) or die "Can't open temp file for writing ($tempdir/temp.bed)\n";
    my $inFh = IO::File->new( $input_file ) or die "Can't open input file\n";
    while( my $line = $inFh->getline )
    {
        next if( $line =~ /^(#|Hugo_Symbol|Chr|chromosome)/ );
        chomp( $line );
        my @F = split( /\t/, $line );

        if( $input_is_maf ){
            $bedFh->print( join( "\t", ( $F[4], $F[5]-1, $F[6] )) . "\n" );
        }
        else {
            $bedFh->print( join( "\t", ( $F[0], $F[1]-1, $F[2] )) . "\n" );
        }
    }
    $inFh->close;
    $bedFh->close;

    # Sort the bed file and run the fast-tier tool on it
    `joinx sort -s -i $tempdir/temp.bed -o $tempdir/temp.bed.sorted`;

    # Tier the loci in the bed file
    my $cmd = "gmt fast-tier fast-tier --tier-file-location $tier_file_location --variant-bed-file $tempdir/temp.bed.sorted --skip-line-count";
    if( defined( $tiering_version )){
        $cmd = $cmd . " --tiering-version $tiering_version";
    }
    my $return = Genome::Sys->shellcmd( cmd => "$cmd" );
    unless($return) {
        $self->error_message("Failed to execute: Returned $return");
        die $self->error_message;
    }

    # Now read in the tier files, match them up with the original bed
    my %tierhash = ();
    # Read in the tier files
    my @tiers = ( "tier1", "tier2", "tier3", "tier4" );
    foreach my $tier ( @tiers ){
        my $tierFh = IO::File->new( "$tempdir/temp.bed.sorted.$tier" ) or die "can't open $tier file\n";
        while( my $line = $tierFh->getline ) {
            chomp( $line );
            my @F = split( /\t/, $line );
            my $key = $F[0] . ":" . ( $F[1] + 1 ) . ":" . $F[2];
            $tierhash{$key} = $tier;
        }
        $tierFh->close;
    }

    # Match up the tiers with the original file, and write an output file annotated with tier info
    my $outFh = IO::File->new( $output_file, ">" ) or die "Can't open file for writing ($output_file)\n";
    $inFh = IO::File->new( $input_file ) or die "Can't open input file\n";
    while( my $line = $inFh->getline ) {

        chomp( $line );        
        # Print out headers with tier appended
        if( $line =~ /^(#|Hugo_Symbol|Chr|chromosome)/ ){
            if($input_is_maf){
                $outFh->print( $line . "\t" . "tier_WU\n" );
            } else {
                $outFh->print( $line . "\t" . "tier\n" );
            }
            next;
        }

        my @F = split( /\t/, $line );
        my $key;
        if( $input_is_maf ){
            $key = $F[4] . ":" . $F[5] . ":" . $F[6];
        }
        else {
            $key = $F[0] . ":" . $F[1] . ":" . $F[2];
        }
        ( defined $tierhash{$key} ) or die "$key not found in hash! Something's not right in the code.\n";
        push( @F, $tierhash{$key} );
        $outFh->print( join( "\t", @F ) . "\n" );
    }
    $inFh->close;
    $outFh->close;

    return 1;
}
