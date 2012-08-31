package Genome::InstrumentData::AlignmentResult::Crossmatch::GTB::Sequencing::ReadIO;

### Adapted from code by Nancy Hansen <nhansen@mail.nih.gov>


############################################################

=head1 NAME

ReadIO.pm - A Perl module to read and write Sequencing
          reads in various formats.

=head1 DESCRIPTION

  A Perl module for inputting/outputting sequencing reads.

=head1 DATE

 March 28, 2009

=head1 AUTHORS

 Nancy Hansen <nhansen@mail.nih.gov>

=head1 PUBLIC METHODS

=cut

############################################################
use strict;
use warnings;

use Carp;
use FileHandle;

use Genome::InstrumentData::AlignmentResult::Crossmatch::GTB::Sequencing::Read;

our $VERSION  = '0.01';

use vars qw( );

###########################################################

=over 4

=item new()

  This method creates a ReadIO object.

  Input:  -file - the name of a file to read from or write
              to.  Format the path passed in "-file" exactly
              as you'd pass it to perl's "open" method, e.g.,
              ">$filename" for writing, "$filename" for 
              reading.  If the file name ends in ".gz", the
              module will open through a pipe to gzip or gunzip.
          -format - a format for read data (fasta, fastq,
             sqr, seq). Default: fasta.  If the format is 'fasta'
             or 'seq', the method will look for an associated
             quality or prb file by appending ".qual" to a 
             fasta file name, or substituting "prb.txt" for 
             "seq.txt" in a seq file name.  When writing, use
             -format 'qual' to write fasta.qual format.
          -qual_offset - the offset in quality values from the
             ASCII vaue of a character used in fastq or SQR 
             format.  E.g., If the offset is 33, the character
             representing 0 quality is '!'.  (Default 33)

  Output: New ReadIO object

=cut

###########################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $file = $params{-file};
    my $format = $params{-format} || 'fasta';
    my $qual_offset = defined($params{-qual_offset}) ? $params{-qual_offset} : 33;

    # see if we need to unzip:

    my $file_to_open = $file;
    if ($file_to_open =~ /\.gz$/)
    {
        $file_to_open = "gunzip -c $file_to_open |";
    }

    # open file to be sure things work:

    my $fh = FileHandle->new($file_to_open)
        or croak "Couldn\'t create filehandle for $file_to_open: $!\n";

    my $self = {file => $file,
                filehandle => $fh,
                format => $format,
                qual_offset => $qual_offset,
                count => 0 };

    bless $self, $class;
    #if ($format eq 'fasta') # open quality file
    #{
#
    #}
    if ($format eq 'seq') # open prb file
    {
        my $prb_file = $file;
        $prb_file =~ s:seq\.txt$:prb.txt:;
        
        my $qh = FileHandle->new($prb_file)
            or croak "Couldn\'t create filehandle for $prb_file: $!\n";

        $self->qual_filehandle($qh);
    }
        
    return $self;
}

###########################################################

=item file()

  This method gets or sets the value of file for this
  ReadIO object.

  Input: Optional argument sets value.
  Output: Value of "file" (scalar string).

=cut

###########################################################
sub file {
    my $self  = shift;

    if (defined (my $new_file = shift))
    {
        $self->{file} = $new_file;
    }

    return $self->{file};

} ## end file

###########################################################

=item filehandle()

  This method gets or sets the value of filehandle for this
  ReadIO object.

  Input: Optional argument sets value.
  Output: Value of "filehandle" (FileHandle object)

=cut

###########################################################
sub filehandle {
    my $self  = shift;

    if (defined (my $new_filehandle = shift))
    {
        $self->{filehandle} = $new_filehandle;
    }

    return $self->{filehandle};

} ## end filehandle

###########################################################

=item qual_filehandle()

  This method gets or sets the value of qual_filehandle for this
  ReadIO object.

  Input: Optional argument sets value.
  Output: Value of "qual_filehandle" (FileHandle object)

=cut

###########################################################
sub qual_filehandle {
    my $self  = shift;

    if (defined (my $new_qual_filehandle = shift))
    {
        $self->{qual_filehandle} = $new_qual_filehandle;
    }

    return $self->{qual_filehandle};

} ## end qual_filehandle

###########################################################

=item format()

  This method gets or sets the value of format for this
  ReadIO object. Possible values are "fasta", "fastq",
  and "sqr" (for now).

  Input: Optional argument sets value.
  Output: Value of "format" (scalar string).

=cut

###########################################################
sub format {
    my $self  = shift;

    if (defined (my $new_format = shift))
    {
        $self->{format} = $new_format;
    }

    return $self->{format};

} ## end format

###########################################################

=item next_read()

  This method reads from the current position in the file
  and returns a GTB::Sequencing::Read object. 

  Input: None.
  Output: GTB::Sequencing::Read object.

=cut

###########################################################
sub next_read {
    my $self  = shift;

    my $format = $self->format();

    my $next_read;

    if ($format =~ /^fasta/i)
    {
        $next_read = $self->_read_fasta();
    }
    elsif ($format eq 'fastq')
    {
        $next_read = $self->_read_fastq();
    }
    elsif ($format eq 'sqr')
    {
        $next_read = $self->_read_sqr();
    }
    elsif ($format eq 'seq')
    {
        $next_read = $self->_read_seq();
    }
    elsif ($format eq 'export')
    {
        $next_read = $self->_read_export();
    }
    else
    {
        croak "Don\'t know how to read format $format in next_read!\n";
    }

    return $next_read;

} ## end next_read

###########################################################

=item write()

  This method writes a GTB::Sequencing::Read object to the
  ReadIO object's filehandle.

  Input: Read object to write to file.
  Output: 1 if successful, undef otherwise

=cut

###########################################################
sub write {
    my $self  = shift;
    my $read_obj = shift
     or croak "Must pass a valid GTB::Sequencing::Read object to write method!\n";
    my %params = @_;
    my $repeat_name = $params{-repeat_name};

    my $fh = $self->filehandle();
    my $format = $self->format();

    my $result; 

    if ($format eq 'fasta')
    {
        $result = $self->_write_fasta($read_obj);
    }
    elsif ($format eq 'fastq')
    {
        $result = $self->_write_fastq($read_obj, -repeat_name => $repeat_name);
    }
    elsif ($format eq 'sqr')
    {
        $result = $self->_write_sqr($read_obj);
    }
    elsif ($format eq 'qual')
    {
        $result = $self->_write_fasta_qual($read_obj);
    }
    else
    {
        croak "Don\'t know how to write format $format in write method!\n";
    }

    return $result;

} ## end write

###########################################################

=item _read_sqr()

  This method reads from the current position in a file of
  sqr format and returns a GTB::Sequencing::Read object. 

  Input: None.
  Output: GTB::Sequencing::Read object.

=cut

###########################################################
sub _read_sqr {
    my $self  = shift;

    my $fh = $self->filehandle();

    my $next_line = <$fh>;

    if (!defined($next_line))
    {
        $fh->close();
        return undef;
    }

    chomp $next_line;
    my ($seq, $qual, $name, $rev_qual, $rev_seq) = split /\t/, $next_line;
    $seq =~ s:\.:N:g;

    #print STDERR "$name, $seq, $qual\n";
    return Genome::InstrumentData::AlignmentResult::Crossmatch::GTB::Sequencing::Read->new(-name => $name,
                                  -seq => $seq,
                                  -qual_string => $qual);
}
 
###########################################################

=item _read_seq()

  This method reads from the current position in a file of
  seq format and returns a GTB::Sequencing::Read object. 

  Input: None.
  Output: GTB::Sequencing::Read object.

=cut

###########################################################
sub _read_seq {
    my $self  = shift;

    my $fh = $self->filehandle();

    my $next_line = <$fh>;

    if (!defined($next_line))
    {
        $fh->close();
        return undef;
    }

    chomp $next_line;
    my ($lane, $tile, $pos1, $pos2, $seq) = split /\t/, $next_line;
    $seq =~ s:\.:N:g;

    $lane = ($lane < 10) ? '00'.$lane :
             ($lane < 100) ? '0'.$lane : $lane;

    $tile = ($tile < 10) ? '00'.$tile :
             ($tile < 100) ? '0'.$tile : $tile;

    my $count = ++$self->{'count'};
    $count = ($count < 10) ? '00'.$count :
             ($count < 100) ? '0'.$count : $count;
    my $name = "$lane\_$tile\_$count";
    my $ra_quals = $self->_read_prb(-seq => $seq);

    return Genome::InstrumentData::AlignmentResult::Crossmatch::GTB::Sequencing::Read->new(-name => $name,
                                  -seq => $seq,
                                  -quals => $ra_quals);

} ## end _read_seq

###########################################################

=item _read_prb()

  This method reads from the current position in a file of
  prb format and returns a reference to an array of quality
  scores.

  Input: -seq passes the sequence of the read.  If the 
      highest Sequencing quality score if for a base other than 
      the base in the sequence, the quality score is converted
      to zero.  Also, all N's are assigned a quality zero.
  Output: Reference to an array of (phred) quality scores.

=cut

###########################################################
sub _read_prb {
    my $self  = shift;
    my %params = @_;
    my $seq = reverse $params{-seq};
    my %index = ('A' => 0, 'C' => 1, 'G' => 2, 'T' => 3);

    my $qh = $self->qual_filehandle();

    my $next_line = <$qh>;

    if (!defined($next_line))
    {
        croak "Different number of lines in prb vs. seq file!\n";
    }

    chomp $next_line;
    my @qual_strings = split "\t", $next_line;

    my $ra_phred_quals = [];
    foreach my $qual_string (@qual_strings)
    {
        my $called_base = chop $seq;
        if ($called_base eq 'N') # automatically quality 0
        {
            push @{$ra_phred_quals}, 0;
            next;
        }

        # if we get here, check that highest quality is called base
        $qual_string =~ s:^\s*(\S.*\S)\s*$:$1:;
        my @four_scores = split /\s+/, $qual_string;
        my $no_quals = @four_scores;

        my $called_q = $four_scores[$index{$called_base}];
        my $qual_pushed = 0;
        foreach my $base qw( A C G T )
        {
            my $index_no = $index{$base};
            if (!defined ($index_no))
            {
                croak "No index defined for base $base!\n";
            }
            if ($four_scores[$index{$base}] > $called_q) # uncalled base has higher score!
            {
                 push @{$ra_phred_quals}, 0;
                 $qual_pushed = 1;
                 last;
            }
        }
        if (!$qual_pushed) # quality not zero--push phred value
        {
            my $phred_score = int($called_q + 10*log(1 + 10**(-$called_q/10))/log(10));
            push @{$ra_phred_quals}, $phred_score;
        }
    }

    return $ra_phred_quals;

} ## end _read_prb

###########################################################

=item _read_fastq()

  This method reads from the current position in a file of
  fastq format and returns a GTB::Sequencing::Read object. 

  Input: None.
  Output: GTB::Sequencing::Read object.

=cut

###########################################################
sub _read_fastq {
    my $self  = shift;

    my $fh = $self->filehandle();

    my $name_line = <$fh>;

    if (!$name_line) # assume end of file:
    {
        $fh->close();
        return undef;
    }

    chomp $name_line;

    my $name = ($name_line =~ /^\@(\S+)/) ? $1 : undef;

    my $seq_line = <$fh>;
    chomp $seq_line;

    my $seq = ($seq_line =~ /^([atgcnATGCN.]+)/) ? $1 : undef;
    # change '.' char to N:

    $seq =~ tr/./N/;

    my $sec_name_line = <$fh>;
    chomp $sec_name_line;

    my $sec_name = ($sec_name_line =~ /^\+(\S+)/) ? $1 : undef;

    my $qual_line = <$fh>;
    chomp $qual_line;

    my $quals = ($qual_line =~ /^(\S+)/) ? $1 : undef;

    if ((!$seq) || (!$quals))
    {
         croak "Unable to parse fastq lines:\n$name_line\n$seq_line\n$sec_name_line\n$qual_line\n";
    }

    return Genome::InstrumentData::AlignmentResult::Crossmatch::GTB::Sequencing::Read->new(-name => $name,
                                  -seq => $seq,
                                  -qual_string => $quals);
}
 
###########################################################

=item _read_export()

  This method reads from the current position in a file of
  export format and returns a GTB::Sequencing::Read object. 

  Input: None.
  Output: GTB::Sequencing::Read object.

=cut

###########################################################
sub _read_export {
    my $self  = shift;

    my $fh = $self->filehandle();
    my $qual_offset = $self->{qual_offset};
    my $offset_char = chr($qual_offset);

    my $line = <$fh>;
    if (!defined ($line)) # assume end of file:
    {
        $fh->close();
        return undef;
    }
    chomp $line;

    my @fields = split /\t/, $line;

    my $name = join ('.', $fields[1], $fields[2], $fields[3], $fields[4], $fields[5] );

    my $seq = $fields[8];
    my $qual_string = $fields[9];

    if ((!$seq) || (!$qual_string))
    {
         croak "Unable to parse line:\n$line\n";
    }

    return Genome::InstrumentData::AlignmentResult::Crossmatch::GTB::Sequencing::Read->new(-name => $name,
                                      -seq => $seq,
                                      -qual_string => $qual_string,
                                      -qual_offset => $offset_char);
}
 
###########################################################

=item _read_fasta()

  This method reads from the current position in a file of
  fasta format and returns a GTB::Sequencing::Read object. 

  Input: None.
  Output: GTB::Sequencing::Read object.

=cut

###########################################################
sub _read_fasta {
    my $self  = shift;

    my $fh = $self->filehandle();

    my $name_line = <$fh>;

    if (!$name_line) # assume end of file:
    {
        $fh->close();
        return undef;
    }

    chomp $name_line;

    my $name = ($name_line =~ /^\@(\S+)/) ? $1 : undef;

    my $seq_line = <$fh>;
    chomp $seq_line;

    my $seq = ($seq_line =~ /^([atgcnATGCN.]+)/) ? $1 : undef;
    # change '.' char to N:

    $seq =~ tr/./N/;

    #my $sec_name_line = <$fh>;
    #chomp $sec_name_line;

    #my $sec_name = ($sec_name_line =~ /^\+(\S+)/) ? $1 : undef;

    #my $qual_line = <$fh>;
    #chomp $qual_line;

    #my $quals = ($qual_line =~ /^(\S+)/) ? $1 : undef;

    #if ((!$seq) || (!$quals))
    #{
         #croak "Unable to parse fastq lines:\n$name_line\n$seq_line\n$sec_name_line\n$qual_line\n";
    #}

    return Genome::InstrumentData::AlignmentResult::Crossmatch::GTB::Sequencing::Read->new(-name => $name,
                                  -seq => $seq);

} # end _read_fasta

###########################################################

=item _write_sqr()

  This method writes the GTB::Sequencing::Read object passed
  as an argument to the object's filehandle in sqr format.

  Input: GTB::Sequencing::Read object.
  Output: 1 if successful, undef otherwise

=cut

###########################################################
sub _write_sqr {
    my $self  = shift;
    my $read_obj = shift;

    my $fh = $self->filehandle();

    if (!defined($fh))
    {
        $fh->close();
        return undef;
    }

    my $seq = $read_obj->seq();
    my $qual_string = $read_obj->qual_string();
    my $name = $read_obj->name();

    print $fh "$seq\t$qual_string\t$name\n"
        or return undef;

    return 1;
}
 
###########################################################

=item _write_fasta()

  This method writes the GTB::Sequencing::Read object passed
  as an argument to the object's filehandle in fasta format.

  Input: GTB::Sequencing::Read object.
  Output: 1 if successful, undef otherwise

=cut

###########################################################
sub _write_fasta {
    my $self  = shift;
    my $read_obj = shift;

    my $fh = $self->filehandle();

    if (!defined($fh))
    {
        $fh->close();
        return undef;
    }

    my $seq = $read_obj->seq();
    my $name = $read_obj->name();

    print $fh ">$name\n$seq\n"
        or return undef;

    return 1;
}
 
###########################################################

=item _write_fasta_qual()

  This method writes the GTB::Sequencing::Read object passed
  as an argument to the object's filehandle in fasta.qual format.

  Input: GTB::Sequencing::Read object.
  Output: 1 if successful, undef otherwise

=cut

###########################################################
sub _write_fasta_qual {
    my $self  = shift;
    my $read_obj = shift;

    my $fh = $self->filehandle();

    if (!defined($fh))
    {
        $fh->close();
        return undef;
    }

    my $name = $read_obj->name();
    my $qualstring = $read_obj->numerical_qual_string();

    print $fh ">$name\n$qualstring\n"
        or return undef;

    return 1;
}
 
###########################################################

=item _write_fastq()

  This method writes the GTB::Sequencing::Read object passed
  as an argument to the object's filehandle in fastq format.

  Input: GTB::Sequencing::Read object.
  Output: 1 if successful, undef otherwise

=cut

###########################################################
sub _write_fastq {
    my $self  = shift;
    my $read_obj = shift;
    my %params = @_;
    my $repeat_name = $params{-repeat_name};

    my $fh = $self->filehandle();

    if (!defined($fh))
    {
        $fh->close();
        return undef;
    }

    my $seq = $read_obj->seq();
    my $name = $read_obj->name() || '';
    my $qualstring = $read_obj->qual_string() || '';

    if ($repeat_name)
    {
        print $fh "\@$name\n$seq\n\+$name\n$qualstring\n"
            or return undef;
    }
    else
    {
        print $fh "\@$name\n$seq\n\+\n$qualstring\n"
            or return undef;
    }

    return 1;
}

###########################################################

=item close_file ()

  The method closes the filehandle associated with the object.

  Input: None.
  Output: Return value of the close method on the FileHandle 
      object.

=cut

###########################################################
sub close_file {
    my $self  = shift;

    my $fh = $self->filehandle(); # returns open fh if opened

    return $fh->close();

} ## end close_file

###########################################################

1;

__END__

=back
