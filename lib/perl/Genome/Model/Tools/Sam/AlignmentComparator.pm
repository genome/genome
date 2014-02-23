package Genome::Model::Tools::Sam::AlignmentComparator;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Basename;

class Genome::Model::Tools::Sam::AlignmentComparator {
    is  => 'Genome::Model::Tools::Sam',
    has => [
        files_to_compare => {
            is_many  => 1,
            shell_args_position => 1,
            doc => 'The sam/bam files to compare',
        },
        fai_file => {
            is => 'String',
            doc => 'FastA index for the reference the reads were aligned against',
        },
        output_dir => {
            is => 'String',
            default_value => '.',
            is_optional => 1,
            doc => 'Directory in which output will be placed',
        },
    ],
};

sub help_brief {
    'Tool to compare BAM or SAM files';
}

sub help_detail {
    return <<EOS
    Tool to compare BAM or SAM files.
EOS
}

sub read_line {
    my $entry = shift;
    my $chr2idx = shift;
    my $fh = $entry->{'in'};
    my $line = scalar(<$fh>);
    return 0 if (not $line) || (substr($line,0,1) eq '@');
    my @fds = split(' ',$line);
    $entry->{'line'} = $line;
    $entry->{'name'} = $fds[0];
    $entry->{'chr'}  = $chr2idx->{$fds[2]};
    $entry->{'pos'}  = $fds[3];
    return 1;
}

sub alignments_agree {
    my ($entry,$data) = @_;
    for my $e (@$data){
        return 0 unless $e->{'pos'} == $entry->{'pos'} && $e->{'name'} eq $entry->{'name'};
    }
    return 1;
}

sub execute {
    my $self = shift;

    my @files = $self->files_to_compare;
    my $fai_file = $self->fai_file;
    my $out_dir = $self->output_dir;
    
    # input checking
    unless (-d $out_dir) {
       $self->error_message("The specified output directory does not exist: $out_dir");
       return;
    }
    unless (-s $fai_file) {
        $self->error_message("The FAI file doesn't exist, or is empty: $fai_file");
        return;
    }
    if (@files < 2) {
        $self->error_message("Not enough files to compare."); 
        return;
    }
    
    # set up the chromosome position hash for easier sorting
    my $fai_fh = IO::File->new("cut -f1 $fai_file |");
    my @chr_lst = map { substr($_,0,-1) } <$fai_fh>;
    $fai_fh->close;
    my %chr2idx;
    @chr2idx{@chr_lst} = (1..scalar(@chr_lst)); # last chr (*) maps to undef
    
    $self->debug_message("Comparing files: ". join(",",@files) );

    my @data = map {
        my ($name,undef,$ext) = fileparse($_,'\.[^\.]*');
        my $in_file = $_;
        if ( $ext =~ /\.bam$/i ){
           die('64-bit OS needed to analyze BAM files') unless (Genome::Config->arch_os() =~ /64/);
           $in_file = "samtools view $_ |"; 
        } elsif ( not $ext =~ /\.sam$/i ){
            $self->error_message("Invalid filetype: $_. Must be SAM or BAM file.");
            return;
        }
        {'in' => IO::File->new($in_file) || die("Error opening $in_file"),
         'out' => IO::File->new(">$out_dir/$name.diff.sam") || die("Error creating >$out_dir/$name.diff.sam"),
        };
    } @files;
    
    for my $entry (@data) {
        die("Encountered a file with no alignment data") if eof $entry->{'in'};
        redo unless read_line($entry,\%chr2idx); # read until the first real line of data
    }
    
    my $old_chr = 0;
    my $old_pos = 0;
    MAIN: while (1) {
        # sort in position order
        @data = sort { $a->{'chr'} <=> $b->{'chr'} || $a->{'pos'} <=> $b->{'pos'} } @data;
        my $entry = shift @data;
        # check for a unanimous alignment, else write out to the 'out' fh
        if (alignments_agree($entry,\@data)) {
            for my $e (@data){ read_line($e,\%chr2idx); }
            for my $e (@data){ 
                unless ($e->{'chr'}){ 
                    push @data, $entry;
                    last MAIN;
                }
            }
        } else {
            my $outfh = $entry->{'out'};
            $outfh->print($entry->{'line'});
        }
        # reset position trackers
        $old_chr = $entry->{'chr'};
        $old_pos = $entry->{'pos'};
        # end when one file runs out
        unless (read_line($entry,\%chr2idx) && $entry->{'chr'}){
            push @data, $entry;
            last;
        }
        # ensure sorted order
        if ($entry->{'chr'} < $old_chr || ($entry->{'pos'} < $old_pos && $entry->{'chr'} == $old_chr)){
            print %$entry,"\n";
            die("Input files must be sorted in position order. Try 'samtools sort'.");
        }
        # re-add the entry to our list
        push @data, $entry;
    }

    # dump any remaining reads and close the handles
    for my $entry (@data) {
        my $infh = $entry->{'in'};
        my $outfh = $entry->{'out'};
        while (<$infh>){
            $outfh->print($_);
        }
        $infh->close;
        $outfh->close;
    }
    $self->debug_message("Comparison finished. Outputs placed in '$out_dir'");

    return 1;
}

1;
