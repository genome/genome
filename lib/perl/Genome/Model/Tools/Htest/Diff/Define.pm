package Genome::Model::Tools::Htest::Diff::Define;

use strict;
use warnings;

use Genome;
use Command;
use Genome::Model;

use Bio::Seq;
use Bio::SeqIO;

class Genome::Model::Tools::Htest::Diff::Define {
    is => 'Command',
    has => [
	from_path => { is => 'String', is_optional => 0, doc => 'Path to the original fasta refseq' },
        to_path   => { is => 'String', is_optional => 0, doc => 'Path to where the new, patched reference in bfa format will be located' },
        changes   => { is => 'String', is_optional => 0, doc => 'Path to the formatted diff file' },
    ],
};

sub help_brief {
    "Load a reference diff into the database";
}

#sub help_synopsis {
#    return <<"EOS"
#    genome-model htest diff define --from-path /my/old/consensus --lst-diff /my/file --to-path /my/new/consensus
#EOS
#}

sub help_detail {                           
    return <<EOS 
Load a reference diff into the database.  See the sub-commands for what formats are acceptable.
EOS
}


sub load_changes_file {
    my $self = shift;
    if (ref($self) eq __PACKAGE__) {
        $self->error_message("This command shouldn't be called directly, use one of its sub-commands");
    } else {
        $self->error_message("Class ",ref($self)," didn't define load_changes_file()");
    }
    return;
}


sub execute {
    my $self = shift;

$DB::single = $DB::stopper;
    my $diff_obj = $self->load_changes_file();  # subclasses for each diff type should define this
    return unless $diff_obj;

    my @diffs = Genome::SequenceDiffPart->get(diff_id => $diff_obj->diff_id);
    #unless (@diffs) {  # no diffs means just copy the original
    #    $self->error_message('No diff parts found for diff_id ',$diff_obj->diff_id,'?!');
    #    return;
    #}

    # Group them up by reference file
    my %diffs_by_refseq;
    foreach my $diff_part ( @diffs ) {
        push (@{$diffs_by_refseq{$diff_part->refseq_path}}, $diff_part);
    }
    @diffs = ();  # save some memory

    my $orig_fasta = Bio::SeqIO->new(-file => $self->from_path, -format => 'fasta');
    my %orig_by_refseq;
    my @orig_refseq_order = ();
    # FIXME When working with the human reference, this is going to be _huge_!
    while (my $seq = $orig_fasta->next_seq) {
        push @orig_refseq_order, $seq->primary_id;
        $orig_by_refseq{$seq->primary_id} = $seq;  # These primary_ids should also be the keys for %diffs_by_refseq
    }

    my $patched_ref_path = $self->to_path . '.fasta';
    my $patched_fasta = Bio::SeqIO->new(-file => ">$patched_ref_path", format => 'fasta');

    foreach my $refseq ( @orig_refseq_order ) {
        my @diff_parts = sort { $a->orig_position <=> $b->orig_position } @{$diffs_by_refseq{$refseq}};
        my $orig_seq = $orig_by_refseq{$refseq};

        if (my @bad = grep { $_->orig_position > $orig_seq->length or $_->orig_position < 1} @diff_parts) {
            $self->error_message(sprintf("Attempt to patch off the end.  patch position %d of refseq %s.  Sequence length %d\n",
                                         $bad[0]->orig_position, $refseq, $bad[0]->length));
            return;
        }
        

        unless (@diff_parts) {  # No diffs for this sequence, just make a copy
            $patched_fasta->write_seq($orig_seq);
            next;
        }

        my $curr_pos = 1;   # Current position in the original file
        my $patched_sequence = '';

        foreach my $diff_part ( @diff_parts ) {
            if ($curr_pos < $diff_part->orig_position) {
                # There's unpatched sequence between the current position and the next diff
                $patched_sequence .= $orig_seq->subseq($curr_pos, $diff_part->orig_position - 1);
                $curr_pos = $diff_part->orig_position;
            }
            
            if ($diff_part->orig_length) {  # This is a deletion
                $curr_pos += $diff_part->orig_length;
            }
            if ($diff_part->patched_length) { # This is an insertion
                $patched_sequence .= $diff_part->patched_sequence;
            }
        }

        # That's all the diffs, write out the rest of the original sequence
        $patched_sequence .= $orig_seq->subseq($curr_pos, $orig_seq->length());

    
        my $seq_obj = Bio::Seq->new(-alphabet => $orig_seq->alphabet,
                                    -desc => $orig_seq->desc . ' ' . $diff_obj->description,
                                    -display_id => $orig_seq->display_id,
                                    -primary_id => $orig_seq->primary_id,
                                    -seq => $patched_sequence);
        $patched_fasta->write_seq($seq_obj);
    }

    $patched_fasta->close();

    # maq needs everything in binary fasta (bfa) format
    my $cmdline = sprintf("maq fasta2bfa %s %s", $patched_ref_path, $self->to_path);
    if (my $retval = system($cmdline)) {
        $retval = $retval >> 8 unless ($retval < 0);
        $self->error_message("running '$cmdline' failed, exit code $retval");
        return;
    }

    # We don't need the new reference in fasta format hanging around any more
    # FIXME maybe make this a File::Temp file after we know this works right 
    # so it'll get cleaned up automaticly.
    unlink($patched_ref_path);

    $self->status_message("diff_id is ".$diff_obj->diff_id);

    return 1;
}

1;

