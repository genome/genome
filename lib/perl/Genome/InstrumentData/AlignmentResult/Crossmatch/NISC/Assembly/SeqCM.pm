
### Adapted from code by Nancy Hansen <nhansen@mail.nih.gov>


############################################################
# SeqCM.pm: Module for parsing and manipulating large
#      crossmatch output files by sequentially reading
#      the alignment strings.
#
# Author:       Nancy F. Hansen
# Version: $Id: SeqCM.pm 3256 2010-01-13 19:37:04Z nhansen $
############################################################

package Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Assembly::SeqCM;
use strict;
use Carp;
use FileHandle;
use Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Assembly::Alignment;

#######################################################################
# SeqCM constructor
#
# INPUT: parameters:
#    -file - file name of a crossmatch output file, or, if not
#            specified, will read from STDIN.
# OUTPUT: a new SeqCM object
#######################################################################

sub new {

    my $this = shift;
    my %params = @_;
    my ($file, $path);

    $file = $params{'-file'};

    if ($file) {
        unless ($file =~ m:^/:) # if no full path provided
        {
            my $dir = `pwd`;
            chomp $dir;
            $file = "$dir/$file";
        }
        $path = $file;
        $path =~ s:^(.*)/[^/]+$:$1:;
    }

    my $self = { file => $file, path => $path };
    my $class = ref($this) || $this;
    bless $self, $class;

    $self->_parse_header(); # open file and check for initial lines

    return($self);
}

###############################################################################
# Subroutine to retrieve the file name of the cm output.
#
# INPUT: SeqCM object
# OUTPUT: scalar string (filename with path)
###############################################################################

sub file_name {

    my $self = shift;

    my $file;
    if (defined ($file = $self->{file}))
    {
        return $file;
    }
    else
    {
        return "";
    }

} # end file_name

###############################################################################
# Subroutine to retrieve the path for the cm output file
#
# INPUT: SeqCM object
# OUTPUT: scalar string (path)
###############################################################################

sub path {

    my $self = shift;

    return $self->{path};

} ## end path

###############################################################################
# Subroutine to retrieve the name of the query_file
#
# INPUT: SeqCM object, optional argument sets value
# OUTPUT: scalar string (crossmatch query_file)
###############################################################################

sub query_file {

    my $self = shift;

    if (defined (my $new_query_file = shift))
    {
        $self->{query_file} = $new_query_file;
    }

    return $self->{query_file};

} ## end query_file

###############################################################################
# Subroutine to retrieve the name of the subject_file
#
# INPUT: SeqCM object, optional argument sets value
# OUTPUT: scalar string (crossmatch subject_file)
###############################################################################

sub subject_file {

    my $self = shift;

    if (defined (my $new_subject_file = shift))
    {
        $self->{subject_file} = $new_subject_file;
    }

    return $self->{subject_file};

} ## end subject_file

###############################################################################
# Subroutine to retrieve the version of crossmatch that was run.
#
# INPUT: SeqCM object, optional argument sets value
# OUTPUT: scalar string (crossmatch version)
###############################################################################

sub version {

    my $self = shift;

    if (defined (my $new_version = shift))
    {
        $self->{version} = $new_version;
    }

    return $self->{version};

} ## end version

###############################################################################
# Subroutine to retrieve the date and time this crossmatch run was done.
#
# INPUT: SeqCM object, optional argument sets value
# OUTPUT: scalar string (date and time in yymmdd:hhmmss format)
###############################################################################

sub date_time {

    my $self = shift;

    if (defined (my $new_date_time = shift))
    {
        $self->{date_time} = $new_date_time;
    }

    return $self->{date_time};

} # end date_time

###############################################################################
# Subroutine to return a list of alignment strings from the output file
#
# INPUT: SeqCM object
# OUTPUT: reference to a list of strings
###############################################################################

sub alignment_strings {

    my $self = shift;

    unless (defined ($self->{alignment_strings}))
    {
        $self->_parse_cm_file();
    }

    return $self->{alignment_strings};    

} ## end alignment_strings

###############################################################################
# Subroutine to return a list of Alignment objects from the output file
#
# INPUT: SeqCM object, optional argument sets value
# OUTPUT: reference to a list of NISC::Assembly::Alignment objects
###############################################################################

sub alignments {

    my $self = shift;

    if (defined (my $new_alignments = shift))
    {
        $self->{alignments} = $new_alignments;
    }

    unless (defined ($self->{alignments}))
    {
        $self->_parse_cm_file();
    }

    return $self->{alignments};    

} ## end alignments

###############################################################################
# Subroutine to parse the header of the crossmatch file, 
#   populating the object's query, subject, alignments_included,
#   and version fields.
#
# INPUT: SeqCM object
# OUTPUT: same objects with fields populated.
###############################################################################

sub _parse_header {

    my $self = shift;
    my $file = $self->file_name();

    if ($self->{'fh'})
    {
        print STDERR "Attempt to open file twice\n";
        return 0;
    }

    my $fh;
    if ($file) {
        my $open_file = ($file =~ /\.gz$/) ? "gunzip -c $file |" : $file;
        $fh = FileHandle->new("$open_file")
            or croak "Couldn\'t open $file for reading: $!\n";
    
    }
    else {
        $fh = FileHandle->new_from_fd(*STDIN, "r")
            or croak "Couldn\'t read from STDIN!\n";
    }

    $self->{'fh'} = $fh; # store for later use
    my $first_line = <$fh>;

    croak "Couldn\'t parse cross_match file!\n" if (!$first_line);

    if ($first_line =~ /alignments/)
    {
        $self->{alignments_included} = 1;
    }
    else
    {
        $self->{alignments_included} = 0;
    }
   
    while (<$fh>)
    { 
        if (m/Query\s+file\(s\):\s*(\S+)\s*$/)
        {
            $self->{query_file} = $1;
            print STDERR "Found query $1\n";
        }
        if (m/Subject\s+file\(s\):\s*(\S+)\s*$/)
        {
            $self->{subject_file} = $1;
            print STDERR "Found subject $1\n";
        }
        if (m/version\s+([\d\.]+)/) # found version
        {
            $self->{version} = $1;
            print STDERR "Found version $1\n";
        }
        if (m/Run\s+date:time\s+(\d{6}:\d{6})/) # found date_time
        {
            $self->{date_time} = $1;
        }
        if (m/Maximal single base matches/)
        {
            last;
        }
    }

    return;

} ## end _parse_header

###############################################################################
# Subroutine to parse the crossmatch file, populating the object's fields.
#
# INPUT: SeqCM object
# OUTPUT: same objects with fields populated.  Note that hit_start > hit_end
#     when alignment is complementary.
###############################################################################

sub _parse_cm_file {

    my $self = shift;

    my $ra_alignments = [];
    while (my $align_obj = $self->next_alignment())
    { 
        push @{$ra_alignments}, $align_obj;
    }

    $self->{alignments} = $ra_alignments;

} ## end _parse_cm_file

###############################################################################
# Subroutine to read a single alignment string out of the file
#
# INPUT: SeqCM object
# OUTPUT: same objects with fields populated.  Note that hit_start > hit_end
#     when alignment is complementary.
###############################################################################

sub next_alignment {

    my $self = shift;
    my $query_file = $self->query_file();
    my $hit_file = $self->subject_file();
    my $fh = $self->{'fh'};

    if (!$fh) # end of file
    {
        return undef;
    }

    my $align_string = ($self->{'first_alignment_line'}) ?
                        $self->{'first_alignment_line'} : '';

    if ($self->{alignments_included})
    {
         my $newline;
         if (!$align_string) # must be beginning of file
         {
             $newline = <$fh>;
             while (defined($newline) && ($newline !~ /^ALIGNMENT/))
             {
                 $newline = <$fh>;
             }
             $align_string = $newline;
         }

         while ((defined ($newline = <$fh>)) && ($newline !~ /^Transitions/))
         {
             $align_string .= $newline;
         }

         while (defined($newline) && ($newline !~ /^ALIGNMENT/))
         {
             $newline = <$fh>;
         }
         $self->{'first_alignment_line'} = $newline;
         if (!$self->{'first_alignment_line'}) # end of file
         {
             $fh->close();
             $self->{'fh'} = undef;
         }
    }
    else # this doesn't make sense to me--will need to be fixed
    {
        while (<$fh>)
        {
            if (/^ALIGNMENT/)
            {
                $align_string = $_;
                last;
            }
        }
    }

    if (!$align_string) # maybe no alignments at all?
    {
        return undef;
    }

    my ($score, $subs, $dels, $ins, $seq_name1, $start1, $end1, 
        $remaining1, $seq_name2, $start2, $end2, $remaining2, $comp);

    if ($align_string =~ m/^(ALIGNMENT\s+){0,1}(\d+)\s+          # score
                            ([\d\.]+)\s+      # substitution %
                            ([\d\.]+)\s+      # deletion %
                            ([\d\.]+)\s+      # insertion %
                            (\S+)\s+          # seq 1 name
                            (\d+)\s+          # seq 1 start
                            (\d+)\s*          # seq 1 end
                            \(\s*(\d+)\s*\)\s+ # seq 1 remaining
                            (\S+)\s+          # seq 2 name
                            (\d+)\s+          # seq 2 start
                            (\d+)\s*          # seq 2 end
                            \(\s*(\d+)\s*\)/xsg) # seq 2 remaining
    {
        ($score, $subs, $dels, $ins, $seq_name1, 
            $start1, $end1, $remaining1, $seq_name2, 
            $start2, $end2, $remaining2, $comp) = 
           ($2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, 'U');
    }
    elsif ($align_string =~ m/^(ALIGNMENT\s+){0,1}(\d+)\s+       # score
                            ([\d\.]+)\s+      # substitution %
                            ([\d\.]+)\s+      # deletion %
                            ([\d\.]+)\s+      # insertion %
                            (\S+)\s+          # seq 1 name
                            (\d+)\s+          # seq 1 start
                            (\d+)\s*          # seq 1 end
                            \(\s*(\d+)\s*\)\s+ # seq 1 remaining
                            C\s+              # comp match
                            (\S+)\s+          # seq 2 name
                            \(\s*(\d+)\s*\)\s* # seq 2 remaining
                            (\d+)\s+          # seq 2 start
                            (\d+)\s*/xsg)      # seq 2 end
    {
        ($score, $subs, $dels, $ins, $seq_name1, 
            $start1, $end1, $remaining1, $seq_name2, 
            $remaining2, $start2, $end2, $comp) = 
           ($2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, 'C');
    }

    my $align_obj = Genome::InstrumentData::AlignmentResult::Crossmatch::NISC::Assembly::Alignment->new(-score => $score, 
              -substitution_rate => $subs, -deletion_rate => $dels,
              -insertion_rate => $ins, -query_file => $query_file,
              -hit_file => $hit_file, -query => $seq_name1, 
              -query_start => $start1, -query_end => $end1, 
              -query_remaining => $remaining1, -hit => $seq_name2,
              -hit_remaining => $remaining2, -hit_start => $start2,
              -hit_end => $end2, -comp => $comp);

    if ($self->{alignments_included})
    {
        # pull alignment strings:

        my $shortened_seq_name1 = $seq_name1;
        $shortened_seq_name1 =~ s:^(.{15}).*$:$1:;
        my $shortened_seq_name2 = $seq_name2;
        $shortened_seq_name2 =~ s:^(.{15}).*$:$1:;

        my @query_strings = ($align_string =~ m/$shortened_seq_name1\s+
                                                \d+\s+([ATGCXNatgcxn-]+)
                                                \s+\d+/xsg);
        my $query_seqstring = join '', @query_strings;
        my @hit_strings = ($align_string =~ m/$shortened_seq_name2\s+
                                                \d+\s+([ATGCXNatgcxn-]+)
                                                \s+\d+/xsg);
        my $hit_seqstring = join '', @hit_strings;

        if ($shortened_seq_name1 eq $shortened_seq_name2) # need to assign separately if query and hit have same name
        {
            my @total_strings = @query_strings;
            my @new_query_strings = ();
            my @new_hit_strings = ();

            for (my $i=0; $i <= $#query_strings; $i++)
            {
                if (2*int($i/2) == $i) # even number matches are query strings
                {
                    push @new_query_strings, $total_strings[$i];
                }
                else
                {
                    push @new_hit_strings, $total_strings[$i];
                }
            }
        
            $query_seqstring = join '', @new_query_strings;
            $hit_seqstring = join '', @new_hit_strings;
        }

        $query_seqstring =~ s:-:*:g;
        $hit_seqstring =~ s:-:*:g;

        $align_obj->query_string($query_seqstring);  
        $align_obj->hit_string($hit_seqstring);  
    }

    return $align_obj;
 
} ## end next_alignment

1;

__END__

=head1 NAME

NISC::Assembly::SeqCM - Perl extension for parsing crossmatch output files.

=head1 SYNOPSIS

use NISC::Assembly::SeqCM;

my $cm_obj = NISC::Assembly::SeqCM->new( -file => $filename );

my $ra_aligns = $cm_obj->alignments();

=head1 DESCRIPTION

This module is used to pull information from a crossmatch output file.

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=head2 new

 Title   : new
 Usage   : $cm_obj = NISC::Assembly::SeqCM->new(-file => $file,
                                             -string = $string);
 Function: Returns a new SeqCM object from a file name or string. 
 Returns : a new SeqCM object
 Args	 : -file	=> a readable filename
	   -string	=> a string in cm format (in place of a filename) 

=head2 file_name

 Title   : file_name
 Usage   : $name = $cm_obj->file_name();
 
 Function: Returns the name of the file containing the cm output.
 Returns : a scalar string (file name with path)
 Args	 : none.

=head2 string

 Title   : string
 Usage   : $cm_string = $cm_obj->string();
 
 Function: Returns the crossmatch-formatted string used to construct the
           SeqCM object, or the contents of the file used to construct
           the SeqCM object.
 Returns : a scalar (string)
 Args	 : none.

=head2 _parse_cm_file

 Title   : _parse_cm_file
 Usage   : $self->_parse_cm_file();
 
 Function: Sets the values of various fields in the crossmatch file by
           parsing the object's string.
 Returns : the object with its new values set.
 Args    : None.

=head1 AUTHOR

Nancy F. Hansen <nhansen@nhgri.nih.gov>

=head1 SEE ALSO

NISC::Assembly::Alignment, NISC::Assembly::ReadInfo.

=cut

