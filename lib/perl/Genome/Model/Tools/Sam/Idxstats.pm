package Genome::Model::Tools::Sam::Idxstats;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sam::Idxstats {
    is => 'Genome::Model::Tools::Sam',
    has => [
        bam_file    => { },
        output_file => { },
        include_stderr => { is => 'Boolean', is_optional => 1, default_value => 0, doc => 'Include any error output from flagstat in the output file.'}
    ],
};

sub execute {
    my $self = shift;

    my $version = $self->use_version;
    ($version) = $version =~ /(\d+)/;

    if ($version < 613) {
        $self->warning_message('samtools version: '.$self->use_version .' do not have idxstats option. Instead r783 will be used');
        $self->use_version('r783');
    }

    my $stderr_redirector = $self->include_stderr ? ' 2>&1 ' : '';
    my $cmd = $self->samtools_path .' idxstats '. $self->bam_file .' > '. $self->output_file . $stderr_redirector;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files =>  [$self->bam_file],
        output_files => [$self->output_file],
    );
    return 1;
}

sub parse_file_into_hashref {
    my ($proto, $stat_file) = @_;

    unless ($stat_file and -s $stat_file) {
        warn "Bam idxstats file: $stat_file is not valid";
        return;
    }

    my $stat_fh = Genome::Sys->open_file_for_reading($stat_file);
    unless($stat_fh) {
        warn 'Fail to open ' . $stat_file . ' for reading';
        return;
    }
    
    my $data;
    my @lines = <$stat_fh>;
    
    for (@lines){
        chomp $_;
        next if /^\*\s+/;
        my ($chr, $rd_len, $map_ct, $unmap_ct) = split /\s+/, $_;
        $data->{$chr} = {
            read_length   => $rd_len,
            map_read_ct   => $map_ct,
            unmap_read_ct => $unmap_ct,
        };
    }
    
    return $data;
}

sub unmap_ref_list {
    my ($proto, $stat_file) = @_;
    my $data = $proto->parse_file_into_hashref($stat_file);

    return [grep{$data->{$_}->{map_read_ct} == 0}keys%$data];
}

sub map_ref_list {
    my ($proto, $stat_file) = @_;
    my $data = $proto->parse_file_into_hashref($stat_file);

    return [grep{$data->{$_}->{map_read_ct} != 0}keys%$data];
}

1;
