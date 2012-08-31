package Genome::Model::Tools::Fasta::ApplyDiff;

use strict;
use warnings;

use Data::Dumper;
use Genome::Utility::FastaStreamOut;
use Genome::Utility::DiffStream;
use Genome::Utility::FastaStreamIn;
use Genome;
use Command;

use IO::File;


UR::Object::Type->define(
    class_name => __PACKAGE__, is => 'Command',
    has => [
        input => { is => 'String', 
            doc => 'Reference fasta file to apply diff too' 
        },
        diff => { is => 'String', 
            doc => 'file containing diffs to apply to fasta.'  
        },
        ],
    has_optional => [
        output => { is => 'String', 
            doc => 'Output file containing input fasta file with appropriate diffs applied', is_optional => 1 
        },
        diff_flank_file => { is => 'String', 
            doc => 'Output file containing flanking sequence of diff positions with the diff applied' 
        },
        ref_flank_file => { is => 'String', 
            doc => 'Output file containing flanking sequence of diff positions original ref sequence' 
        },
        flank_size => { is => 'Int', 
            doc => 'If exporting flank files, this is the amount of flanking sequence to include around each diff position', 
            default => 0 
        },
        ignore_deletion_sequence => { is => 'Boolean', 
            doc => "Flag to skip indel sequence comparison with reference fasta sequence.  Useful if diff file doesn't contain any/complete seqence for deletions",
            default => 0
        },
        column_width => { is => 'Int',
            doc => "column width for output fasta column width, defaults to 60",
        },
    ],
);

sub help_brief {
    "applies seq inserts and deletes from a diff file to a fasta file";
}

sub help_detail {                           # This is what the user will see with --help <---
    return q/ 
    Given a reference fasta file and a diff file, applies the diffs to the reference. Can produce several outputs: The original reference fasta with the diffs applied, and fastas containing sequence regions around the diff positions with or without the diff applied.
    Usage:
    genome-model tools apply-diff-to-fasta --input reference.fasta --diff reference_indels.diff --output modified_reference.fasta --diff_flank_file flanks_with_diff.fasta --ref_flank_file original_flanks.fasta --flank_size 50 --ignore_deletion_sequence

    This command streams through both files simultaneously, so for your diffs to be applied correctly, the following criteria must be met:
    The headers in your diffs must be in the same order as the headers in the reference fasta file.
    Diffs for same header must be sorted positionally within the diff file.

    Diff file format: <fasta_header> <position> <indel_type> <deletion_length> <indel_sequence> <pre_diff_sequence> <post_diff_seqeunce>  
    The deletion length field should be the lenght of the deletion sequence.  If the indel is an insertion, this value should be zero. If the indel is a substitution, this value should contain the length of the bases removed.  The fasta_header is simply the first sequence of non-whitespace characters after the '>' on the fasta section description line.  Also see Genome::Utility::DiffStream docs

    Remember to set --flank-size if producing flanking fasta files, or you will get empty or very small regions in your output files.  --flank-size is set to a minimum length of pre- or post- diff sequences in the diff file.
    /
}


sub execute {
    my $self = shift;

#ERROR HANDLE ON FAILURE
    $self->error_messages_callback(
        sub{
            $self->status_message("Removing output files due to error");
            unlink $self->output; 
            unlink $self->diff_flank_file if $self->diff_flank_file;
            unlink $self->ref_flank_file if $self->ref_flank_file;
        }
    );

#INPUT/OUTPUT STREAMS;
    my $fasta_stream    = Genome::Utility::FastaStreamIn->new(IO::File->new("< ".$self->input));
    my $diff_stream     = Genome::Utility::DiffStream->new(IO::File->new("< ".$self->diff));

    my $diff_flank_stream;
    my $ref_flank_stream;
    my $output_stream;
    $diff_flank_stream = Genome::Utility::FastaStreamOut->new(IO::File->new("> ".$self->diff_flank_file), $self->column_width) if $self->diff_flank_file;
    $ref_flank_stream =  Genome::Utility::FastaStreamOut->new(IO::File->new("> ".$self->ref_flank_file), $self->column_width) if $self->ref_flank_file;
    $output_stream   = Genome::Utility::FastaStreamOut->new(IO::File->new("> ".$self->output ), $self->column_width) if $self->output;

#LOOP VARIABLES
    my $pre_diff_sequence;
    my $post_diff_sequence;
    my $skip_diff;

    my $ignore_del = $self->ignore_deletion_sequence;

    my $read_position = 0;
    my $buffer;
    my $write_position = 0;

    my $flank_header = 1; 
    my $ref_flank_sequence;
    my $diff_flank_sequence;

    my $flank_size = $self->flank_size;
    my $min_flank_size; #used to at least grab pre and post diff sequence if they are present
    my $left_flank_position;
    my $right_flank_position;

    my $first_part_length;
    my $first_part;

    my $current_fasta_header = '';
    my $current_fasta_header_id='';

    my $successful_diffs = 0;
    my $attempted_diffs;
    my @failed_diffs;

###########################################################
#MAIN LOOP

    $DB::single = $DB::stopper;
    while (my $diff = $diff_stream->next_diff){

        # leftover buffer from previous diff
        $output_stream->print($buffer) if $output_stream;
        $buffer = undef;

        #generate flanking positions
        do { 
            no warnings;
            $min_flank_size = length $diff->{pre_diff_sequence} <=> length $diff->{post_diff_sequence} 
                ?  length $diff->{post_diff_sequence} 
                : length $diff->{pre_diff_sequence}
        };
        $flank_size             = $min_flank_size if $min_flank_size > $flank_size;
        $left_flank_position    = $diff->{position} - $flank_size;
        $left_flank_position    = 0 if $left_flank_position < 0;
        $right_flank_position   = $diff->{position} + ( $diff->{deletion_length} ) + $flank_size;

        #ADVANCE THROUGH THE FASTA FILE UNTIL CURRENT FASTA SECTION HEADER EQ DIFF HEADER
        until ($current_fasta_header_id eq $diff->{header}) { 

            #print out entire fasta section if still not at diff->{header} in fasta
            while ($buffer = $fasta_stream->next_line){
                $output_stream->print($buffer) if $output_stream;
            }
            $buffer = undef;

            #grab next header and reset write position
            unless ( $current_fasta_header = $fasta_stream->next_header ){
                use Data::Dumper;
                $self->error_message( "Can't get next fasta header and we still have diffs to process!\n".Dumper $diff);
                return;
            }

            $write_position = 0;
            $read_position = 0;
            $output_stream->print_header($current_fasta_header) if $output_stream;
            $current_fasta_header_id = $self->parse_header($current_fasta_header);
        }        

        #ERROR CHECK
        if ($write_position > $diff->{position}) {
            $self->error_message(
                "Write position is greater than diff(header:".$diff->{header}.") position! We've missed the boat! $write_position > ".$diff->{position}
            );
            return;
        } 

        #ADVANCE FASTA TO LEFT FLANK OF CURRENT DIFF

        while( $write_position <= $left_flank_position){ 
            unless (defined $buffer){
                $buffer = $fasta_stream->next_line;
                unless ($buffer){ 
                    $self->error_message(
                        "hit the end of the section and haven't reached the current diff's(header:".$diff->{header}.") left flank position! $write_position < ".$left_flank_position
                    );
                    return;
                }
                $read_position = $read_position + length $buffer;
            }
            last if $read_position >= $left_flank_position; #current buffered line contains start of flank
            $output_stream->print($buffer) if $output_stream;
            $buffer = undef;
            $write_position = $read_position;
        }

        #ADVANCE THROUGH BUFFER UNTIL FLANK POSITION
        $first_part_length = $left_flank_position - $write_position;
        $first_part = substr($buffer, 0,  $first_part_length, ''); 

        $output_stream->print($first_part) if $output_stream;

        $write_position += $first_part_length;

        $diff_flank_sequence = ''; 
        $ref_flank_sequence = '';

        #STARTING HERE WE CAN REPEAT PROCESS UNTIL WE'VE REACHED THE END OF THE FLANK SECTION IF RIGHT FLANK OF THE CURRENT DIFF OVERLAPS THE LEFT FLANK OF THE NEXT DIFF
        while( defined $diff_flank_sequence and defined $ref_flank_sequence ){

            while( $write_position <= $diff->{position}){ #advance fasta to diff position, flank printing starts here
                unless (defined $buffer){
                    $buffer = $fasta_stream->next_line;
                    unless ($buffer){ #fail condition
                        $self->error_message(
                            "Hit the end of the section and haven't reached the current diff's (header:".$diff->{header}.") position! $write_position < ".$diff->{position}
                        );
                        return;
                    }
                    $read_position = $read_position + length $buffer;
                }

                last if $read_position >= $diff->{position}; #current buffered line contains diff pos

                $output_stream->print($buffer) if $output_stream; #otherwise,
                $diff_flank_sequence .= $buffer; 
                $ref_flank_sequence .= $buffer;
                $buffer = undef; 
                $write_position = $read_position;
            }

            $first_part_length = $diff->{position} - $write_position;
            $first_part = substr($buffer, 0,  $first_part_length, '');  #this splices out the first part from the buffer string

            $output_stream->print($first_part) if $output_stream;
            $write_position += $first_part_length;
            $diff_flank_sequence .= $first_part;
            $ref_flank_sequence .= $first_part;

            if ( $diff->{pre_diff_sequence} ) {
                $pre_diff_sequence = substr($ref_flank_sequence, (0 - length( $diff->{pre_diff_sequence} ) ) );

                unless ( length ($pre_diff_sequence) eq length ($diff->{pre_diff_sequence}) ){
                    $self->status_message("Pre diff sequences are different lengths! Diff not processed!");
                }


                unless ( $pre_diff_sequence eq $diff->{pre_diff_sequence} ) {
                    $self->status_message("pre_diff_sequence in diff(".$diff->{type}."-".$diff->{deletion_length}." header:".$diff->{header}." pos:".$diff->{position}.") does not match actual fasta sequence! fasta:$pre_diff_sequence not eq diff:".$diff->{pre_diff_sequence}."  Diff not processed!" );
                    $skip_diff = 1;
                }
            }

            if ( $diff->{post_diff_sequence} and !$skip_diff ) {
                my $length = length ( $diff->{post_diff_sequence} );
                my $del_length = ( $diff->{deletion_length} );
                while ( ($length + $del_length) > length $buffer ){
                    my $nextline = $fasta_stream->next_line;
                    unless ($nextline){
                        $self->error_message("length of post diff(header:".$diff->{header}." pos:".$diff->{position}.") sequence goes beyond end of file!");
                        return;
                    }
                    $buffer.=$nextline;
                    $read_position += length $nextline;
                }
                $post_diff_sequence = substr($buffer, $del_length, $length);
                unless ( $post_diff_sequence eq $diff->{post_diff_sequence} ){
                    $self->status_message("post_diff_sequence in diff(".$diff->{type}."-".$diff->{deletion_length}." header:".$diff->{header}." pos:".$diff->{position}.") does not match actual fasta sequence! fasta:$post_diff_sequence not eq diff:".$diff->{post_diff_sequence}."  Diff not processed!" );
                    $skip_diff = 1;
                }
            }

            if ($diff->{insert} and !$skip_diff){
                $output_stream->print($diff->{insert}) if $output_stream;
                $diff_flank_sequence .= $diff->{insert};  #diff flank gets the insert
            }

#REF FLANK GETS THE DELETE, OUTPUT_STREAM DOESNT SEQUENCE
            if ($diff->{delete} and !$skip_diff) { 
                my $to_delete = $diff->{delete};
                while ( $diff->{deletion_length} > length($buffer) ) {  
                    my $nextline = $fasta_stream->next_line;
                    unless ($nextline){ #fail condition
                        $self->error_message("deletion($to_delete) substring extends beyond end of sequence, ending ($buffer)");
                        return;
                    }
                    $buffer .= $nextline;
                    $read_position += length $nextline;
                }
                my $deletion = substr($buffer, 0, $diff->{deletion_length}, ''); #check buffer at this point
                $ref_flank_sequence.=$deletion; 

                unless ($ignore_del or $deletion eq $to_delete){ 
                    $self->status_message( "deleted seq does not equal actual sequence! $deletion != $to_delete  chrom: $current_fasta_header_id position: $write_position $first_part - $deletion - $buffer ". ($diff->{position} + 1) ."\n");
                    $skip_diff = 1;
                }
                $write_position += length $deletion unless $skip_diff;
            }

            $successful_diffs++ unless $skip_diff;
            
            my ($next_diff_position, $next_pre_diff_length) = $diff_stream->next_diff_position;
            $next_pre_diff_length||=0;
            my $next_flank_length = $flank_size >= $next_pre_diff_length ? $flank_size : $next_pre_diff_length;
            if (  $next_diff_position and $right_flank_position > $next_diff_position - $next_flank_length ){ #now check if the tail of the flank overlaps the next diff flank
                $diff = $diff_stream->next_diff;
                $min_flank_size         = length $diff->{pre_diff_sequence} <=> length $diff->{post_diff_sequence} ? 
                length $diff->{post_diff_sequence} : 
                length $diff->{pre_diff_sequence};
                $flank_size             = $min_flank_size if $min_flank_size > $flank_size;
                $left_flank_position    = $diff->{position} - $flank_size;
                $left_flank_position    = 0 if $left_flank_position < 0;
                $right_flank_position   = $diff->{position} + ( $diff->{deletion_length} ) + $flank_size;
                $skip_diff = 0;

            }else{

#ADVANCE FASTA TO DIFF POSITION, FLANK PRINTING STARTS HERE
                while( $write_position <= $right_flank_position){ 

                    unless (defined $buffer){
                        $buffer = $fasta_stream->next_line;
                        unless ($buffer){ 
                            $buffer = ''; #end of section
                            $right_flank_position = $write_position;
                        }
                        $read_position = $read_position + length $buffer;
                    }

                    last if $read_position >= $right_flank_position; #current buffered line contains diff pos

                    $output_stream->print($buffer) if $output_stream; #otherwise, add to flank and print buffer
                    $diff_flank_sequence .= $buffer; 
                    $ref_flank_sequence .= $buffer;
                    $buffer = undef;
                    $skip_diff = 0;
                    $write_position = $read_position;
                }

                $first_part_length = $right_flank_position - $write_position;
                $first_part = substr($buffer, 0,  $first_part_length, '');  #this splices out the first part from the buffer string

                $output_stream->print($first_part) if $output_stream;
                $write_position += $first_part_length;
                $diff_flank_sequence .= $first_part;
                $ref_flank_sequence .= $first_part;

                if ( length $diff_flank_sequence and $diff_flank_stream){
                    $diff_flank_stream->print_header(">$current_fasta_header_id|".$diff->{position}." (Diff cluster #$flank_header, flank size $flank_size)");
                    $diff_flank_stream->print($diff_flank_sequence);
                }
                if ( length $ref_flank_sequence and $ref_flank_stream){
                    $ref_flank_stream->print_header(">$current_fasta_header_id|".$diff->{position}." (Diff cluster #$flank_header, flank size $flank_size)");
                    $ref_flank_stream->print($ref_flank_sequence);
                }
                $diff_flank_sequence = undef;
                $ref_flank_sequence = undef;
                $flank_header++;
            }
        }
    }
#FINISH PRINTING BUFFER
    if ($output_stream){
        $output_stream->print($buffer);
        while (my $line= $fasta_stream->next_line){
            $output_stream->print($line);
        }

#FINISH PRINTING FASTA FILE
        while (my $current_header = $fasta_stream->next_header){
            $output_stream->print_header($current_header);
            while (my $line= $fasta_stream->next_line){
                $output_stream->print($line);
            }
        }
    }

    $self->status_message("Diffs processed successfully: $successful_diffs");

    $output_stream->close() if $output_stream;
    $ref_flank_stream->close() if $ref_flank_stream;
    $diff_flank_stream->close() if $diff_flank_stream;

    return 1;
}

sub parse_header{
    my $self = shift;
    my $line = shift;
    my ($header) = $line =~ />(\S+)/;
    # my $header = substr($line,1,1);
    return $header;
}

1;

#$HeadURL$
#$Id$
