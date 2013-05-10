package Genome::Model::Tools::LiftOver;

use warnings;
use strict;
use IO::File;
use Genome;

class Genome::Model::Tools::LiftOver {
    is => 'Command',
    has_input => [
        source_file => {
            is => 'Text',
            doc => 'File containing loci to be remapped to a new genome build with the UCSC liftOver tool',
        },
        destination_file => {
            is => 'Text',
            doc => 'Where to output the remapped loci',
        },
    ],
    has_optional => [
        unmapped_file => {
            is => 'Text',
            doc => 'Where to output the unmappable loci',
        },
        allow_multiple_output_regions => {
            is => 'Boolean', default => '0',
            doc => 'Whether or not to allow multiple output regions',
        },
        file_format => {
            is => 'Text', default => 'bed',
            doc => 'The format of the source file',
            valid_values => ['bed', 'gff', 'genePred', 'sample', 'pslT'],
        },
        input_is_annoformat => {
            is => 'Boolean', default => 0,
            doc => 'Input is 1-based annotation format',
        },

        input_is_maf_format => {
            is => 'Boolean', default => 0,
            doc => 'Input is 1-based maf format',
        },

        input_is_vcf_format => {
            is => 'Boolean', default => 0,
            doc => 'Input is 1-based maf format',
        },

        chain_file => {
            is => 'Text',
            doc => 'The liftOver "chain" file that maps from source to destination build. Required if --lift-direction is unspecified',
        },
        lift_direction => {
            is => 'Text',
            doc => 'Shorthand for commonly used lift operations.',
            valid_values => ['hg18ToHg19', 'hg19ToHg18'],
            example_values => ['hg18ToHg19', 'hg19ToHg18'],
        },
    ],
    doc => "Wrapper for the UCSC liftOver with added conveniences",
};

sub execute {
    my $self = shift;
    my ( $source_file, $dest_file, $unmapped_file ) = ( $self->source_file, $self->destination_file, $self->unmapped_file );
    my ( $chain_file, $lift_direction ) = ( $self->chain_file, $self->lift_direction );
    my $file_format = ( $self->file_format eq 'bed' ? "" : " -" . $self->file_format ),
    my $multiple_regions = ( $self->allow_multiple_output_regions ? " -multiple" : "" );

    unless( defined $chain_file ) {
        if( $lift_direction eq 'hg19ToHg18' ) {
            $chain_file = "/gscmnt/gc12001/info/model_data/chain_files/hg19Tohg18/hg19ToHg18.over.chain";
        }
        elsif( $lift_direction eq 'hg18ToHg19' ) {
            $chain_file = "/gscmnt/gc12001/info/model_data/chain_files/hg18Tohg19/hg18ToHg19.over.chain";
        }
        # Invalid values of --lift-direction are already handled
    }

    # Allow unmapped loci to be an optional output, by providing liftOver a dummy file to write to
    $unmapped_file = Genome::Sys->create_temp_file_path unless( defined $unmapped_file );
    ( $unmapped_file ) or die "Unable to create temporary file $!";

    # If input is annotation format, convert to a bed for liftOver, and convert back to anno later
    my $anno_headers = "";
    #### anno ####
    if( $self->input_is_annoformat ) {
        # Create temp directory for munging, and replace the source/dest files that liftOver will use
        my $tempdir = Genome::Sys->create_temp_directory;
        ( $tempdir ) or die "Unable to create temporary directory $!";
        ( $source_file, $dest_file ) = ( "$tempdir/inbed", "$tempdir/outbed" );

        my $outFh = IO::File->new( $source_file, ">" ) or die "Can't open temp file $source_file";
        my $inFh = IO::File->new( $self->source_file ) or die "Can't open file $self->source_file";
        while( my $line = $inFh->getline ) {

            # Save header lines if any, to write to output later
            if( $line =~ m/^(#|chromosome_name)/ ) {
                $anno_headers .= $line;
                next;
            }
            chomp( $line );
            my @F = split( /\t/, $line );

            # tabbed format: 1  123  456  A  T

            # we're combining all the other fields into the name field since liftover is picky
            # about its format. I've chosen a delimiter "?|?" unlikely to be found in any real
            # files (I hope). Feel free to update this to something less hacky later on.
            my $extraFields = join( "?|?", ( @F[3..$#F] ));

            # spaces also get treated as delimiters by liftover, so we'll replace them with
            # something equally unlikely
            $extraFields =~ s/ /?_?/g;

            # liftover doesn't like insertion coordinates (start==stop), so we add one to the stop
            # to enable a liftover, and remove it on the other side. This is effectively the same
            # as leaving it unchanged. For deletions and SNVs, decrement the start locus
            --$F[1] unless( $F[3] =~ m/^[-0*]$/ );
            $outFh->print( join( "\t", ( "chr$F[0]", $F[1], $F[2], $extraFields )) . "\n" );
        }
        $inFh->close;
        $outFh->close;
    }
    #### maf ####
    elsif( $self->input_is_maf_format ) {
        # Create temp directory for munging, and replace the source/dest files that liftOver will use
        my $tempdir = Genome::Sys->create_temp_directory;
        ( $tempdir ) or die "Unable to create temporary directory $!";
        ( $source_file, $dest_file ) = ( "$tempdir/inbed", "$tempdir/outbed" );

        my $outFh = IO::File->new( $source_file, ">" ) or die "Can't open temp file $source_file";
        my $inFh = IO::File->new( $self->source_file ) or die "Can't open file $self->source_file";
        while( my $line = $inFh->getline ) {

            # Save header lines if any, to write to output later
            if( $line =~ m/^(#|Hugo_Symbol)/ ) {
                $anno_headers .= $line;
                next;
            }
            chomp( $line );
            my @F = split( /\t/, $line );

            # we're combining all the other fields into the name field since liftover is picky
            # about its format. I've chosen a delimiter "?|?" unlikely to be found in any real
            # files (I hope). Feel free to update this to something less hacky later on.
            my $extraFields = join( "?|?", ( @F[0..3],@F[7..$#F] ));

            # spaces also get treated as delimiters by liftover, so we'll replace them with
            # something equally unlikely
            $extraFields =~ s/ /?_?/g;

            # liftover doesn't like insertion coordinates (start==stop), so we add one to the stop
            # to enable a liftover, and remove it on the other side. This is effectively the same
            # as leaving it unchanged. For deletions and SNVs, decrement the start locus
            --$F[5] unless( $F[10] =~ m/^[-0*]$/ );
            $outFh->print( join( "\t", ( "chr$F[4]", $F[5], $F[6], $extraFields )) . "\n" );
        }
        $inFh->close;
        $outFh->close;
    }
    #### vcf ####
    elsif( $self->input_is_vcf_format ) {
        # Create temp directory for munging, and replace the source/dest files that liftOver will use
        my $tempdir = Genome::Sys->create_temp_directory;
        ( $tempdir ) or die "Unable to create temporary directory $!";
        ( $source_file, $dest_file ) = ( "$tempdir/inbed", "$tempdir/outbed" );

        my $outFh = IO::File->new( $source_file, ">" ) or die "Can't open temp file $source_file";
        my $inFh = IO::File->new( $self->source_file ) or die "Can't open file $self->source_file";
        while( my $line = $inFh->getline ) {

            # Save header lines if any, to write to output later
            if( $line =~ m/^#/ ) {
                #update reference to avoid confusion
                if( $line =~ /##reference=/ ){
                    $line = "##reference=Lifted with $chain_file\n";
                }
                $anno_headers .= $line;
                next;
            }
            chomp( $line );
            my @F = split( /\t/, $line );

            # we're combining all the other fields into the name field since liftover is picky
            # about its format. I've chosen a delimiter "?|?" unlikely to be found in any real
            # files (I hope). Feel free to update this to something less hacky later on.
            my $extraFields = join( "?|?", ( @F[2..$#F] ));

            # spaces also get treated as delimiters by liftover, so we'll replace them with
            # something equally unlikely
            $extraFields =~ s/ /?_?/g;

            # Since everything has one coordinate in VCF, we're going to treat everything like a
            # SNP for simplicity
            $F[1]--;
            $outFh->print( join( "\t", ( "chr$F[0]", $F[1]++, $F[1], $extraFields )) . "\n" );
        }
        $inFh->close;
        $outFh->close;
    }

    # Do the lifting over
    my $cmd = 'liftOver';
    $cmd .= "$multiple_regions$file_format ";
    $cmd .= join( ' ', ( $source_file, $chain_file, $dest_file, $unmapped_file ));
    print "RUN: $cmd\n";
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$source_file, $chain_file],
        output_files => [$dest_file],
    );

    # Convert back to annotation format if necessary
    if( $self->input_is_annoformat ) {
        my $outFh = IO::File->new( $self->destination_file, ">" ) or die "can't open output file: $!";
        $outFh->print( $anno_headers ); # Print header lines copied from the input file
        my $inFh = IO::File->new( $dest_file ) or die "can't open file\n";
        while( my $line = $inFh->getline ) {
            chomp( $line );
            my @F = split( "\t", $line );
            my @extra = split( /\?\|\?/, $F[3] );
            my $ref = $extra[0];

            # Restore trailing fields
            my $post = join( "\t", @extra );
            $post =~ s/\?_\?/ /g;

            $F[0] =~ s/^chr//g;
            # Increment start locus for deletions and SNVs, but not for insertions
            ++$F[1] unless( $ref =~ m/^[-0*]$/ );
            $outFh->print( join( "\t", ( $F[0], $F[1], $F[2], $post )) . "\n" );
        }
        $inFh->close;
        $outFh->close;

    }
    # Convert back to maf format if necessary
    elsif( $self->input_is_maf_format ) {
        my $outFh = IO::File->new( $self->destination_file, ">" ) or die "can't open output file: $!";
        $outFh->print( $anno_headers ); # Print header lines copied from the input file
        my $inFh = IO::File->new( $dest_file ) or die "can't open file\n";
        while( my $line = $inFh->getline ) {
            chomp( $line );
            my @F = split( "\t", $line );
            my @extra = split( /\?\|\?/, $F[3] );
            my $ref = $extra[7];

            # Restore non-coordinate fields
            my $pre = join( "\t", @extra[0..3] );
            $pre =~ s/\?_\?/ /g;

            my $post = join( "\t", @extra[4..$#extra] );
            $post =~ s/\?_\?/ /g;

            $F[0] =~ s/^chr//g;
            # Increment start locus for deletions and SNVs, but not for insertions
            ++$F[1] unless( $ref =~ m/^[-0*]$/ );
            $outFh->print( join( "\t", ( $pre, $F[0], $F[1], $F[2], $post )) . "\n" );
        }
        $inFh->close;
        $outFh->close;

    }
    # Convert back to vcf format if necessary
    elsif( $self->input_is_vcf_format ) {
        my $outFh = IO::File->new( $self->destination_file, ">" ) or die "can't open output file: $!";
        $outFh->print( $anno_headers ); # Print header lines copied from the input file
        my $inFh = IO::File->new( $dest_file ) or die "can't open file\n";
        while( my $line = $inFh->getline ) {
            chomp( $line );
            my @F = split( "\t", $line );
            my @extra = split( /\?\|\?/, $F[3] );
            my $post = join( "\t",@extra );
            $post =~ s/\?_\?/ /g;

            $F[0] =~ s/^chr//g;
            $outFh->print( join( "\t", ( $F[0], $F[1]+1, $post )) . "\n" );
        }
        $inFh->close;
        $outFh->close;
    }

    return 1;
}

1;
