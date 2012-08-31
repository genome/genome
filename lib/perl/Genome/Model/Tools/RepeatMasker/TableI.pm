package Genome::Model::Tools::RepeatMasker::TableI;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RepeatMasker::TableI {
    is => ['Command'],
    is_abstract => 1,
    has => [
        output_table => {
            is => 'Text',
            doc => 'This is the file location where the output table from this command will be written.  This is NOT a RepeatMasker(.tbl) file!',
        },
        _total_count => {
            is => 'Integer',
            is_optional => 1,
        },
        _total_bp => {
            is => 'Integer',
            is_optional => 1,
        },
        _aligned => {
            is => 'Integer',
            is_optional => 1,
        },
        _repeat_aligned => {
            is => 'Integer',
            is_optional => 1,
        },
    ],
};

sub print_table_from_hash_ref {
    my $self = shift;
    my $repeats_hash_ref = shift;

    my $masked_bp = 0;
    my $string = '';

    my %repeats = %{$repeats_hash_ref};
    for my $family (sort keys %repeats) {
        my $family_bp = delete($repeats{$family}{base_pair});
        $masked_bp += $family_bp;
        my $family_elements = delete($repeats{$family}{elements});
        my $family_pc = sprintf("%.02f",(($family_bp / $self->_total_bp ) * 100)) .'%';
        $string .= $family .":\t". $family_elements ."\t". $family_bp ."\t". $family_pc."\n";
        for my $class (sort keys %{$repeats{$family}}) {
            my $class_elements = $repeats{$family}{$class}{elements};
            my $class_bp = $repeats{$family}{$class}{base_pair};
            my $class_pc = sprintf("%.02f",(($class_bp / $self->_total_bp ) * 100)) .'%';
            if ($class) {
                $string .= "\t". $class .":\t". $class_elements ."\t". $class_bp ."\t". $class_pc."\n";
            }
        }
        $string .= "\n";
    }

    my $table_fh = Genome::Sys->open_file_for_writing($self->output_table);
    unless ($table_fh) {
        die('Failed to open table file for output: '. $self->output_table);
    }
    print $table_fh "sequences:\t". $self->_total_count ."\n";
    print $table_fh "total length:\t". $self->_total_bp ."\n";
    if ($self->_aligned) {
       print $table_fh "aligned:\t". $self->_aligned ."\n";
    }
    if ($self->_repeat_aligned) {
        print $table_fh "repeat aligned:\t". $self->_repeat_aligned ."\n";
    }
    print $table_fh "masked:\t". $masked_bp ." bp ( ". sprintf("%02f",(($masked_bp / $self->_total_bp ) * 100)) ." %) \n";
    print $table_fh $string ."\n";
    $table_fh->close;
    return 1;
}

sub print_samples_summary_from_hash_ref {
    my $self = shift;
    my $samples_ref = shift;

    my %samples = %{$samples_ref};

    my %ignore = (
        'sequences' => 1,
        'base_pair' => 1,
        'masked' => 1,
        'aligned' => 1,
        'repeat_aligned' => 1,
    );

    my %families;
    for my $sample (sort keys %samples) {
        for my $family (grep { !defined($ignore{$_}) } keys %{$samples{$sample}}) {
            $families{$family} = {};
            for my $class (grep { !defined($ignore{$_}) } keys %{$samples{$sample}{$family}}) {
                $families{$family}{$class} = 1;
            }
        }
    }

    my %data;
    for my $sample (sort keys %samples) {
        push @{$data{'samples'}}, $sample;
        push @{$data{'sequences'}}, delete($samples{$sample}{'sequences'});
        push @{$data{'base_pair'}}, delete($samples{$sample}{'base_pair'});
        push @{$data{'masked'}}, delete($samples{$sample}{'masked'});
        push @{$data{'aligned'}}, delete($samples{$sample}{'aligned'}) if $samples{$sample}{'aligned'};
        push @{$data{'repeat_aligned'}}, delete($samples{$sample}{'repeat_aligned'}) if $samples{$sample}{'repeat_aligned'};
        for my $family (keys %families) {
            if ($samples{$sample}{$family}) {
                push @{$data{$family}{base_pair}}, delete($samples{$sample}{$family}{base_pair});
            } else {
                push @{$data{$family}{base_pair}},0;
            }
            for my $class (keys %{$families{$family}}) {
                if ($samples{$sample}{$family}{$class}) {
                    push @{$data{$family}{$class}{base_pair}}, delete($samples{$sample}{$family}{$class}{base_pair});
                } else {
                    push @{$data{$family}{$class}{base_pair}}, 0;
                }
            }
        }
    }


    my $table_fh = Genome::Sys->open_file_for_writing($self->output_table);
    unless ($table_fh) {
        die('Failed to open table file for output: '. $self->output_table);
    }
    print $table_fh "samples:\t". join("\t",@{$data{'samples'}}) ."\n";
    print $table_fh "sequences:\t". join("\t",@{$data{'sequences'}}) ."\n";
    if ($data{'aligned'}) {
        print $table_fh "aligned:\t". join("\t",@{$data{'aligned'}}) ."\n";
    }
    if ($data{'repeat_aligned'}) {
        print $table_fh "repeat aligned:\t". join("\t",@{$data{'repeat_aligned'}}) ."\n";
    }
    my @base_pair = @{$data{'base_pair'}};
    print $table_fh "total length:\t". join("\t",@base_pair) ."\n";

    my $masked_array_ref = $data{'masked'};
    print $table_fh 'masked:';
    for (my $i = 0; $i < scalar @{$masked_array_ref}; $i++) {
        my $total_bp = $base_pair[$i];
        my $masked_bp = $masked_array_ref->[$i];
        print $table_fh "\t". $masked_bp .'('. sprintf("%.02f",(($masked_bp/$total_bp) * 100)) .'%)';
    }
    print $table_fh "\n\n";
    for my $family (sort keys %families) {
        my $family_bp_array_ref = delete($data{$family}{base_pair});
        print $table_fh $family .':';
        for (my $i = 0; $i < scalar @{$family_bp_array_ref}; $i++) {
            my $total_bp = $base_pair[$i];
            my $family_bp = $family_bp_array_ref->[$i];
            print $table_fh "\t". $family_bp .'('. sprintf("%.02f",(($family_bp/$total_bp) * 100)) .'%)';
        }
        print $table_fh "\n";
        for my $class (sort keys %{$families{$family}}) {
            my $class_bp_array_ref = delete($data{$family}{$class}{base_pair});
            print $table_fh $class .':';
            for (my $i = 0; $i < scalar @{$class_bp_array_ref}; $i++) {
                my $total_bp = $base_pair[$i];
                my $class_bp = $class_bp_array_ref->[$i];
                print $table_fh "\t". $class_bp .'('. sprintf("%.02f",(($class_bp/$total_bp) * 100)) .'%)';
            }
            print $table_fh "\n";
        }
        print $table_fh "\n";
    }
    $table_fh->close;
    return 1;
}

1;
