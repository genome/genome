package Genome::Model::Tools::Snp::Sort;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Sort::Naturally qw| nsort |;

class Genome::Model::Tools::Snp::Sort {
    is => 'Command',
    has => [
        snp_file => {
            type => 'String',
            is_optional => 0,
            doc => "maq cns2snp output",
            shell_args_position => 1,
        },
        ],
    has_optional => [
        output_file => {
            type => 'Text',
            doc => 'optional output file',
        },

        force => {
            type => 'Boolean',
            doc => 'force overwriting of output file, even if it already exists',
            default => 0,
        },


        ],
};

sub help_brief {
    "Sorts a SNP file using Sort::Naturally to sort the chromosomes";
}

sub help_detail {
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless $self;

    unless (Genome::Sys->validate_file_for_reading($self->snp_file)) {
        $self->error_message('Failed to validate snp file '. $self->snp_file .'  for reading.');
        return;
    }
    unless ($self->force){
        if ($self->output_file) {
            unless (Genome::Sys->validate_file_for_writing($self->output_file)) {
                $self->error_message('Failed to validate output file '. $self->output_file .' for writing.');
                return;
            }
        }
    }
    return $self;
}

sub execute {
    my $self=shift;

    my $snp_fh = Genome::Sys->open_file_for_reading($self->snp_file);
    my $output_fh;
    if ($self->output_file) {
        if($self->force){
            $output_fh =  Genome::Sys->open_file_for_overwriting($self->output_file);
        } else {
            $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
        }
    } else {
        $output_fh = IO::Handle->new;
        $output_fh->fdopen(fileno(STDOUT),'w');
    }

    my %snp_at;
    while(my $line = $snp_fh->getline) {
        my ($chr, $pos,) = split /\t/, $line;
        push @{$snp_at{$chr}{$pos}}, $line;  #in case sometimes got multiple lines with the same chr, pos
    }
    $snp_fh->close;
    for my $chr (nsort keys %snp_at) {
        for my $pos (sort { $a <=> $b } keys %{$snp_at{$chr}}) {
            map{$output_fh->print($_)}@{$snp_at{$chr}{$pos}};
        }
    }
    $output_fh->close;
    return 1;
}


1;




