package Genome::Model::Tools::Htest::Diff::Define::Lst;

use strict;
use warnings;

use Genome;
use Command;
use Genome::Model;

use IO::File;

class Genome::Model::Tools::Htest::Diff::Define::Lst {
    is => 'Genome::Model::Tools::Htest::Diff::Define',
};

sub help_brief {
    "Load a reference diff in .lst format into the database and create a new, patched reference from it ";
}

sub help_synopsis {
    return <<"EOS"
    genome-model htest diff define --from-path /my/old/consensus --changes /my/file --to-path /my/new/consensus
EOS
}

sub help_detail {                           
    return <<EOS 
Load a reference diff into the database.  See the sub-commands for what formats are acceptable.
EOS
}

# split the changes file into buckets by refseq_path
sub _bucketize_changes_by_refseq {
    my($self,$changes_filename) = @_;

    my $fh = IO::File->new($changes_filename);
    unless ($fh) {
        $self->error_message("Can't open changes file $changes_filename: $!");
        return;
    }

    my $changes_by_refseq = {};
    while(<$fh>) {
        chomp;
        my($refseq_name, $position, $repl_seq, $orig_seq, $confidence) = split;

        unless ($refseq_name && $position && $repl_seq && $orig_seq && $confidence) {
            $self->error_message("Couldn't parse line ",$fh->input_line_number," of .lst file ",$self->changes);
            return;
        }

        if ($confidence eq '+') {
            $confidence = 99;
        } elsif ($confidence eq '-') {
            $confidence = 1;
        }

        $orig_seq = '' if ($orig_seq eq '-');
        $repl_seq = '' if ($repl_seq eq '-');

        my $change_type;
        if ($orig_seq and $repl_seq) {
            $change_type = 'snp';
        } elsif ($orig_seq or $repl_seq) {
            $change_type = 'indel';
        } else {
            $self->error_message("orig_seq and repl_seq are both blank on line  ",$fh->input_line_number," of .lst file ",$self->changes);
            return
        }

        push @{$changes_by_refseq->{$refseq_name}}, { refseq_name => $refseq_name,
                                                      position => $position,
                                                      orig_seq => uc $orig_seq,
                                                      repl_seq => uc $repl_seq,
                                                      confidence => $confidence,
                                                      change_type => $change_type,
                                                    };
    }

    # Sort the changes by position
    foreach my $refseq_name ( keys %$changes_by_refseq ) {
        my @sorted = sort {$a->{'position'} <=> $b->{'position'}} @{$changes_by_refseq->{$refseq_name}};
        $changes_by_refseq->{$refseq_name} = \@sorted;
    }

    return $changes_by_refseq;
}




sub load_changes_file {
    my $self = shift;

    my $changes_by_refseq = $self->_bucketize_changes_by_refseq($self->changes);

    my $diff_obj = Genome::SequenceDiff->create(from_path => $self->from_path, to_path => $self->to_path, description => 'imported from .lst format');
    

    # group contiguous changes together
    # Do we need to worry about an insertion and a deletion overlapping?
    foreach my $refseq_name ( keys %$changes_by_refseq ) {
        my $patched_offset = 0;  # How different the orig and patched positions are
        my $last_position = 0;
        my %accumulate;  # Holds the accumulated, contiguous indel information

        foreach my $change ( @{$changes_by_refseq->{$refseq_name}} ) {
            no warnings 'uninitialized';

            if ($change->{'change_type'} ne 'snp' and 
                 ($last_position == $change->{'position'}  # part of some kind of overlapping indel
                  or $last_position + 1 == $change->{'position'})) {  # or otherwise contiguous

                $accumulate{'orig_seq'} .= $change->{'orig_seq'};
                $accumulate{'repl_seq'} .= $change->{'repl_seq'};
                $accumulate{'position'} = $change->{'position'} unless defined($accumulate{'position'});
                $accumulate{'confidence'} = $change->{'confidence'} unless defined($accumulate{'confidence'});

            } else {
                # Not contiguous, save what we've got and start a new contiguous region
                if ($accumulate{'position'}) {
                    my $diff_part = $self->_combine_and_store_accumulations(refseq_name => $refseq_name,
                                                                            diff_id => $diff_obj->id,
                                                                            patched_offset => $patched_offset,
                                                                            %accumulate);
                    $patched_offset += $diff_part->orig_length - $diff_part->patched_length;

                }

                if ($change->{'change_type'} eq 'snp') {
                    # we have to handle these by themselves, I think
                    $self->_combine_and_store_accumulations(refseq_name => $refseq_name,
                                                            diff_id => $diff_obj->diff_id,
                                                            patched_offset => $patched_offset,
                                                            %$change);
                    %accumulate = ();
                } else {
                    %accumulate = ( position => $change->{'position'},
                                    confidence => $change->{'confidence'},
                                    orig_seq => $change->{'orig_seq'},
                                    repl_seq => $change->{'repl_seq'});
                }
            }

            $last_position = $change->{'position'};
        }
        # Ran through the whole list, save the last contiguous part
        $self->_combine_and_store_accumulations(refseq_name => $refseq_name,
                                                diff_id => $diff_obj->id,
                                                patched_offset => $patched_offset,
                                                %accumulate) if ($accumulate{'position'});
    }

                                                                
    return $diff_obj;
}

sub _combine_and_store_accumulations {
    my($self,%args) = @_;

    # FIXME what do we do with the confidence for sequences longer than 1?
    my $diff_part = Genome::SequenceDiffPart->create(diff_id => $args{'diff_id'},
                                                     refseq_path => $args{'refseq_name'},
                                                     orig_position => $args{'position'} + $args{'patched_offset'},
                                                     orig_length => length($args{'orig_seq'}),
                                                     orig_sequence => $args{'orig_seq'},
                                                     patched_position => $args{'position'},
                                                     patched_length => length($args{'repl_seq'}),
                                                     patched_sequence => $args{'repl_seq'},
                                                     confidence_value => $args{'confidence'},
                                                   );
    return $diff_part;
}

1;

