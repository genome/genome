package Genome::Model::Tools::LiftOver;

use warnings;
use strict;
use IO::File;
use Genome;

use File::Spec;

class Genome::Model::Tools::LiftOver {
    is => 'Command::V2',
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
    has_optional_input => [
        unmapped_file => {
            is => 'Text',
            doc => 'Where to output the unmappable loci',
        },
        chain_file => {
            is => 'Text',
            doc => 'The liftOver "chain" file that maps from source to destination build. Required if --lift-direction is unspecified',
        },
    ],
    has_optional_param => [
        lift_direction => {
            is => 'Text',
            doc => 'Shorthand for commonly used lift operations.',
            valid_values => ['hg18ToHg19', 'hg19ToHg18'],
            example_values => ['hg18ToHg19', 'hg19ToHg18'],
        },
        allow_multiple_output_regions => {
            is => 'Boolean', default => '0',
            doc => 'Whether or not to allow multiple output regions',
        },

        # TODO: auto-detect this if not specified
        file_format => {
            is => 'Text', default => 'bed',
            doc => 'The format of the source file',
            valid_values => ['bed', 'gff', 'genePred', 'sample', 'pslT','bedpe'],
        },

        # TODO: make these into additional input formats
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
    ],
    doc => "wrapper for the UCSC liftOver tool with support for additional input formats, maintaining additional columns",
};

sub help_synopsis {
    return <<EOS
# liftover TGI legacy SV annotation 
gmt lift-over \
    --file-format bedpe \
    --source  /gscmnt/gc13003/info/test_suite_data//Genome-Model-Tools-LiftOver/2013-07-25-sv-bedpe/inputs/svs.bedpe \
    --dest /tmp/lifted \
    --unmapped /tmp/unlifted \
    --lift-direction hg18ToHg19
EOS
}

sub help_detail {
    return <<EOS
The UCSC liftOver tool translates coordinates from one reference genome to another.
This wrapper supports a wider variety of file formats, and also handles maintaining additional columns present in the source file.
EOS
}

our %FORMAT_TRANSLATES_TO_LIFTOVER_PARAM = map { $_ => 1 } qw/gff genePred sample pslT/;

sub execute {
    my $self = shift;
    my ( $source_file, $dest_file, $unmapped_file ) = ( $self->source_file, $self->destination_file, $self->unmapped_file );
    my ( $chain_file, $lift_direction ) = ( $self->chain_file, $self->lift_direction );

    # TODO: some formats are determined by those is_* flags instead of this field.
    # Change that.
    my $file_format = $self->file_format;

    # Some formats translate to a liftOver parameter, others we handle manually here. 
    my $file_format_param = ($FORMAT_TRANSLATES_TO_LIFTOVER_PARAM{$file_format} ? "" : " -" . $file_format );

    my $multiple_regions = ( $self->allow_multiple_output_regions ? " -multiple" : "" );

    unless( defined $chain_file ) {
        if( $lift_direction eq 'hg19ToHg18' ) {
            my $alloc = Genome::Disk::Allocation->get(allocation_path => 'model_data/chain_files/hg19Tohg18');
            $self->fatal_message("Couldn't find chain file for hg19Tohg18") unless $alloc;

            $chain_file = File::Spec->join($alloc->absolute_path, "hg19ToHg18.over.chain");
        }
        elsif( $lift_direction eq 'hg18ToHg19' ) {
            my $alloc = Genome::Disk::Allocation->get(allocation_path => 'model_data/chain_files/hg18Tohg19');
            $self->fatal_message("Couldn't find chain file for hg19Tohg18") unless $alloc;

            $chain_file = File::Spec->join($alloc->absolute_path, "hg18ToHg19.over.chain");
        }
        else {
            die "no chain file specified, and not lift direction specified from which we can derive a chain file!";
        }
        # Invalid values of --lift-direction are already handled
    }

    # Allow unmapped loci to be an optional output, by providing liftOver a dummy file to write to
    $unmapped_file = Genome::Sys->create_temp_file_path unless( defined $unmapped_file );
    ( $unmapped_file ) or die "Unable to create temporary file $!";

    # this was cut-and-pasted below multiple times
    my $tempdir;
    my $inFh;
    my $outFh;
    my $init_tmpfiles = sub {
        $tempdir = "/tmp/";
        $tempdir = Genome::Sys->create_temp_directory;
        $tempdir or die "Unable to create temporary directory $!";
        ($source_file, $dest_file) = ( "$tempdir/inbed", "$tempdir/outbed" );
        unlink $source_file if -e $source_file;
        unlink $dest_file if -e $dest_file;     
        $outFh = Genome::Sys->open_file_for_writing($source_file);
        $inFh = Genome::Sys->open_file_for_reading($self->source_file);
    };

    # If input is annotation format, and we convert to a bed for liftOver then back to anno later, preserve the header
    my $anno_headers = "";

    if ($file_format eq 'bedpe') {
        #### svpair
        $init_tmpfiles->(); 
        my $n = 0;
        while (my $line = $inFh->getline) {
            chomp $line;
            $n++;
            my @F = split(/\t/, $line );
            # liftover doesn't like insertion coordinates (start==stop), so we add one to the stop
            # to enable a liftover, and remove it on the other side. 
            $outFh->print( join( "\t", ( "chr$F[0]", $F[1], $F[2]+1, "$n.1")) . "\n" );
            $outFh->print( join( "\t", ( "chr$F[3]", $F[4], $F[5]+1, "$n.2")) . "\n" );
        }
        $inFh->close();
        $outFh->close();
    }
    elsif( $self->input_is_annoformat ) {
        #### anno ####
        # Create temp directory for munging, and replace the source/dest files that liftOver will use
        $init_tmpfiles->();
        while( my $line = $inFh->getline ) {

            # Save header lines if any, to write to output later
            if( $line =~ m/^(#|chromosome_name)/ ) {
                $anno_headers .= $line;
                next;
            }
            chomp( $line );
            my @F = split( /\t/, $line );

            # tabbed format: 1  123  456  A  T

            #Some older annotation files may have start/stop coordinates in wrong order. If so, reverse:
            if ($F[1]>$F[2]){
              @F[1,2]=@F[2,1];
            }

            my $extraFields = _merge_fields(@F[3..$#F]);

            # liftover doesn't like insertion coordinates (start==stop), so we add one to the stop
            # to enable a liftover, and remove it on the other side. This is effectively the same
            # as leaving it unchanged. For deletions and SNVs, decrement the start locus
            --$F[1] unless( $F[3] =~ m/^[-0*]$/ );
            $outFh->print( join( "\t", ( "chr$F[0]", $F[1], $F[2], $extraFields )) . "\n" );
        }
        $inFh->close;
        $outFh->close;
    }
    elsif( $self->input_is_maf_format ) {
        #### maf ####
        $init_tmpfiles->(); 
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
    elsif( $self->input_is_vcf_format ) {
        #### vcf ####
        $init_tmpfiles->();
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

            my $extraFields = _merge_fields(@F[2..$#F]);

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
    $cmd .= "$multiple_regions$file_format_param ";
    $cmd .= join( ' ', ( $source_file, $chain_file, $dest_file, $unmapped_file ));
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$source_file, $chain_file],
        output_files => [$dest_file],
    );

    # Convert back to original format if we converted *into* be format just to work with liftover
    if ($file_format eq 'bedpe') {
        # This is an old SV annotation format in the form CHR1 POS1 POS1 CHR2 POS2 POS2 TYPE ORIENTATION
        # Because POS1 appears twice and POS2 appears twice (always a zero-width range), the new format
        # was trimmed-down.  TODO: add the new format.
        my $origFh = IO::File->new($self->source_file);
        my $liftedFh = IO::File->new( $dest_file ) or die "can't open file\n";
        my $outFh = IO::File->new( $self->destination_file, ">" ) or die "can't open output file: $!";
        
        # the unmapped file will be rebuilt while comparint the inputs to the outputs 
        unlink $self->unmapped_file;
        my $unmappedFh = IO::File->new( $self->unmapped_file, ">" ) or die "can't open unmapped file: $!";

        my $orig_n = 0;
        while(my $lifted_line1 = $liftedFh->getline) {
            chomp($lifted_line1);
            my @lifted_cols1 = split("\t", $lifted_line1);
            $lifted_cols1[2]--;
            my ($lifted_from_original_row1, $colset1) = split('\.',$lifted_cols1[3]);
       
            unless ($lifted_from_original_row1 =~ /\d/ and $colset1 =~ /^(1|2)/) {
                $DB::single=1;
                die "bad row $lifted_line1";
            }

            # get the line in the original file that goes with this lifted line
            # if there are rows in the original file that come before the lifted one 
            # they will go to the unmapped as we pass over them
            my $orig_line;
            my @orig_cols;
            for (1) {
                $orig_line = $origFh->getline;
                unless ($orig_line) {
                    die "out of rows in the original file while content is still in the liftover file???!: $lifted_line1\n";
                }
                chomp $orig_line;
                @orig_cols = split("\t", $orig_line);
                $orig_n++;

                if ($lifted_from_original_row1 < $orig_n) {
                    die "lifted file not in order?";
                }
                elsif ($orig_n < $lifted_from_original_row1) {
                    # this line in the original file was not in the lifted file
                    $unmappedFh->print($orig_line,"\n");
                    redo;
                }
            }

            if ($colset1 eq '2') {
                # we expect the first row to have a value of '1'
                # if this is 2, the first row did not map
                # consider the whole thing unmapped
                $unmappedFh->print($orig_line,"\n");
                next;
            }
            elsif ($colset1 ne '1') {
                die "expected the field to have n.1 or n.2!";
            }

            my $lifted_line2 = $liftedFh->getline;
            chomp($lifted_line2);
            my @lifted_cols2 = split("\t", $lifted_line2);
            $lifted_cols2[2]--;
            my ($lifted_from_original_row2, $colset2) = split('\.',$lifted_cols2[3]);

            if ($lifted_from_original_row2 != $lifted_from_original_row1) {
                # there was no line 2 in the file for this original row
                # we just went on to the liftover data for another row
                # consider this pair unmapped
                $unmappedFh->print($orig_line,"\n");
                # and consider the lifted data again from scratch
                $lifted_line1 = $lifted_line2;
                redo;
            }
            elsif ($colset2 ne '2') {
                die "expected the extra col to be n.1 or n.2!";
            }

            $lifted_cols1[0] =~ s/chr//;
            $lifted_cols2[0] =~ s/chr//;
            $outFh->print(join("\t",@lifted_cols1[0,1,2],@lifted_cols2[0,1,2],@orig_cols[6..$#orig_cols]),"\n");
        }
        $origFh->close;
        $liftedFh->close;
        $outFh->close;
    }
    elsif( $self->input_is_annoformat ) {
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

sub _merge_fields {
    my @F = @_;

    # we're combining all the other fields into the name field since liftover is picky
    # about its format. I've chosen a delimiter "?|?" unlikely to be found in any real
    # files (I hope). Feel free to update this to something less hacky later on.
    my $extraFields = join( "?|?", ( @F ));

    # spaces also get treated as delimiters by liftover, so we'll replace them with
    # something equally unlikely
    $extraFields =~ s/ /?_?/g;

    return $extraFields;
}

1;
