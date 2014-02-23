package Genome::Model::Tools::Fasta::Chunk;

use strict;
use warnings;

use Genome;
use IO::File;
use Bio::SeqIO;
use File::Temp;
use File::Basename;

class Genome::Model::Tools::Fasta::Chunk {
    is  => 'Genome::Model::Tools::Fasta',
    has => [
        chunk_size => {
            is  => 'Integer',
            doc => 'Number of query fasta for each chunk',
        },
    ],
    has_optional => [
        chunk_dir => {
            is  => 'String',
            doc => 'Root directory of temp dir to hold chunk fasta(qual) files, default is dir of fasta_file',
        },
        show_list => {
            is  => 'boolean',
            doc => 'print file path of fasta/qual chunk file fof to std out',
            default => 0,
        },
    ],
};


sub help_brief {
    'Divide fasta into chunk by chunk_size' 
}


sub help_detail {  
    return <<EOS
    Divide the fasta(qaul) file of multi-fasta(qual) into chunk by given chunk_size. --show-list option will show the file path of chunk file list.
EOS
}


sub create {
    my $class = shift;
    
    my $self = $class->SUPER::create(@_);
    return unless $self;
    my $root_dir = $self->chunk_dir || dirname $self->fasta_file;
    
    my $chunk_dir = File::Temp::tempdir(
        "FastaChunkDir_XXXXXX", 
        DIR => $root_dir,
    );
    $self->chunk_dir($chunk_dir);

    return $self;
}
        
    
sub execute {
    my $self = shift;
    $self->dump_status_messages($self->show_list);

    my $fa_in_io = $self->get_fasta_reader($self->fasta_file);
    my $qa_in_io = $self->get_qual_reader($self->qual_file) if $self->have_qual_file;
    
    my $seq_ct  = 0;
    my $file_ct = 0;

    my ($fa_out_io, $qa_out_io, $chunk_fa_file, $chunk_qa_file);
    my (@chunk_fa_files, @chunk_qa_files);

    my $chunk_dir = $self->chunk_dir;

    while (my $seq = $fa_in_io->next_seq) {
        my $qual = $qa_in_io ->next_seq if $self->have_qual_file;
        $seq_ct++;
        
        if ($seq_ct > $self->chunk_size || !defined $chunk_fa_file) {
            $seq_ct = 1;
            $file_ct++;
            
            $chunk_fa_file = $chunk_dir."/Chunk$file_ct.fasta";
            $fa_out_io = $self->get_fasta_writer($chunk_fa_file);
            push @chunk_fa_files, $chunk_fa_file;

            if ($self->have_qual_file) {
                my ($f_id, $q_id) = ($seq->id, $qual->id);
                
                unless ($f_id eq $q_id) {
                    $self->error_message("id of fasta and quality not equal: $f_id <=> $q_id");
                    return;
                }
                unless ($seq->length == $qual->length) {
                    $self->error_message("length of fasta and quality not equal: $f_id <=> $q_id");
                    return;
                }
                
                my $chunk_qa_file = $chunk_fa_file.'.qual';
                $qa_out_io = $self->get_qual_writer($chunk_qa_file);
                push @chunk_qa_files, $chunk_qa_file;
            }
        }
        $fa_out_io->write_seq($seq);
        $qa_out_io->write_seq($qual) if $self->have_qual_file;
    }
    
    my $rv = $self->_write_to_file($chunk_dir, 'fasta', \@chunk_fa_files);
    return unless $rv;
    $rv = $self->_write_to_file($chunk_dir, 'qual', \@chunk_qa_files) if $self->have_qual_file;
    return unless $rv;
    
    my @out = (\@chunk_fa_files);
    push @out, \@chunk_qa_files if $self->have_qual_file;

    $self->_set_chunk_dir_permissions();
    return \@out;
}

# File::Temp sets the mode of the created directory to 700.
# This changes it to match the user's umask
sub _set_chunk_dir_permissions {
    my $self = shift;

    my $umask = umask();
    my $mode = 0777 & ( ~ $umask);
    chmod $mode, $self->chunk_dir();
}


sub _write_to_file {
    my ($self, $dir, $type, $files) = @_;
    
    my $fof_file = $dir.'/chunk_'.$type.'_file.fof';
    my $fh = IO::File->new(">$fof_file") or
        ($self->error_message("can't write to $fof_file") and return);
    map{$fh->print($_."\n")}@$files;
    $fh->close;
    $self->debug_message("List of chunk $type files:\n$fof_file");

    return 1;
}
    
    
1;

#$HeadURL$
#$Id$
