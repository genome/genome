package Genome::Model::Tools::Sam::VarFilter;

use strict;
use warnings;

use Genome;
use File::Temp;

my %common_options = (
    min_map_qual_snp => {
        is  => 'Integer',
        doc => '-Q minimum RMS mapping quality for SNPs',
        default => 10,
    },
    min_read_depth => {
        is  => 'Integer',
        doc => '-d minimum read depth',
        default => 2,
    },
    max_read_depth => {
        is  => 'Integer',
        doc => '-D maximum read depth',
        default => 10000000,
    },
    gap_win_size => {
        is  => 'Integer',
        doc => '-W window size for filtering adjacent gaps',
        default => 10,
    },
    gap_nearby_size => {
        is  => 'Integer',
        doc => '-w SNP within INT bp around a gap to be filtered',
        default => 3,
    },
);

my %vcf_options = (
    min_alt_bases  => {
        is  => 'Integer',
        doc => '-a minimum number of alternate bases',
        default => 2,
    },
    min_strand_bias => {
        is  => 'Float',
        doc => '-1 min P-value for strand bias (given PV4)',
        default => 1e-4,
    },
    min_baseQ_bias  => {
        is  => 'Float',
        doc => '-2 min P-value for baseQ bias',
        default => 1e-100,
    },
    min_mapQ_bias   => {
        is  => 'Float',
        doc => '-3 min P-value for mapQ bias',
        default => 0,
    },
    min_end_dis_bias=> {
        is  => 'Float',
        doc => '-4 min P-value for end distance bias',
        default => 1e-4,
    },
    #samtools 0.1.16 vcfutils.pl does not have this -e option
    #min_hwe     => {
    #    is  => 'Float',
    #    doc => '-e min P-value for HWE (plus F<0)',
    #    default => 1e-4,
    #},
    max_gq      => {
        is  => 'Integer',
        doc => '-G maximum number of GQ, comparing MXGQ',
        default => 0,
    },
    max_sp      => {
        is  => 'Integer',
        doc => '-S maximum number of SP, comparing MXSP',
        default => 1000,
    },
);

my %pileup_options = (
    min_map_qual_gap => {
        is  => 'Integer',
        doc => '-q minimum RMS mapping quality for gaps, default 10',
        default => 10,
    },
    snp_win_size => {
        is  => 'Integer',
        doc => '-W window size for filtering dense SNPs, default 10',
        default => 10,
    },
    max_snp_per_win => {
        is  => 'Integer',
        doc => '-N maximum number of SNPs in a sized window',
        default => 2,
    },
    min_indel_score => {
        is  => 'Integer',
        doc => '-G minimum indel score for nearby SNP filtering, default is 25',
        default => 25,
    },
);

my %other_options = (
    snv_out_file   => {
        is  => 'String',
        doc => 'snv output file after filter',
    },
    filtered_snv_out_file => {
        is  => 'String',
        doc => 'filtered snv output file',
        default => 'snv.varfilter.filtered',
    },
    indel_out_file => {
        is  => 'String',
        doc => 'indel output file after filter',
    },
    filtered_indel_out_file => {
        is  => 'String',
        doc => 'filtered indel output file',
        default => 'indel.varfilter.filtered',
    },
    input_var_file   => {
        is  => 'FilePath',
        doc => 'input var file, vcf or pileup',
    },
    input_var_format => {
        is  => 'String',
        doc => 'input var file format either vcf or pileup',
        default      => 'bcf',
        valid_values => ['bcf', 'vcf', 'pileup'],
    },
);

#for the common options, vcf uses diff param from pileup, 
my %opt_pu_conv = (
    gap_win_size => '-l',
);


class Genome::Model::Tools::Sam::VarFilter {
    is  => 'Genome::Model::Tools::Sam',
    has_optional => [%common_options, %vcf_options, %pileup_options, %other_options],
};


sub help_brief {
    'Filter samtools mpileup/pileup snp indel output.';
}

sub help_detail {
    return <<EOS
    Filter samtools-mpileup snp indel output.
EOS
}


sub execute {
    my $self = shift;
    my $input_file         = $self->input_var_file;
    my $input_var_format   = $self->input_var_format;
    my $snv_out_file       = $self->snv_out_file;
    my $indel_out_file     = $self->indel_out_file;
    my $flt_snv_out_file   = $self->filtered_snv_out_file;
    my $flt_indel_out_file = $self->filtered_indel_out_file;

    unless (-e $input_file) {
        $self->error_message("Input var file: $input_file not existing");
        return;
    }
    
    unless ($snv_out_file or $indel_out_file) {
        $self->error_message('Either snv_out_file or indel_out_file is needed');
        return;
    }

    my ($tool_path, @fhs);

    if ($input_var_format eq 'bcf') {
        $tool_path = $self->path_for_bcftools($self->use_version);
        $tool_path .= " view $input_file | " . $self->path_for_vcfutils($self->use_version);
    }
    elsif ($input_var_format eq 'vcf') {
        $tool_path = $self->path_for_vcfutils($self->use_version);
    }
    elsif ($input_var_format eq 'pileup') {
        $tool_path = $self->samtools_pl_path;
    }
    else {
        $self->error_message("Invalid input var format $input_var_format");
        return;
    }
    
    my $cmd = $tool_path . ' varFilter';

    if ($input_var_format =~ /cf$/) {
        $cmd .= $self->_get_cmd_opts(%common_options, %vcf_options);
    }
    else {
        $cmd .= $self->_get_cmd_opts(%common_options, %pileup_options);
    }

    my $tmp_file     = Genome::Sys->create_temp_file_path('tmp_var.'.$input_var_format);
    my $flt_tmp_file = Genome::Sys->create_temp_file_path('tmp_var_filtered.'.$input_var_format);
   
    if ($input_var_format eq 'bcf') { #stdout for good vars, stderr for filtered vars
        $cmd .= " - 1> $tmp_file 2> $flt_tmp_file";
    }
    else {
        $cmd .= " $input_file 1> $tmp_file 2> $flt_tmp_file";
    }

    my $rv = Genome::Sys->shellcmd(
        cmd => $cmd,
        output_files => [$tmp_file, $flt_tmp_file],
        skip_if_output_is_present => 0,
    );
        
    unless ($rv == 1) {
        $self->error_message("Failed to run command: $cmd");
        return;
    }

    my ($tmp_snv_fh, $tmp_snv_file, $tmp_flt_snv_fh, $tmp_flt_snv_file, $tmp_indel_fh, $tmp_indel_file, $tmp_flt_indel_fh, $tmp_flt_indel_file);

    if ($snv_out_file) {
        ($tmp_snv_fh, $tmp_snv_file)         = Genome::Sys->create_temp_file('tmp_snv.'.$input_var_format);
        ($tmp_flt_snv_fh, $tmp_flt_snv_file) = Genome::Sys->create_temp_file('tmp_snv_filtered.'.$input_var_format);
    }

    if ($indel_out_file) {
        ($tmp_indel_fh, $tmp_indel_file)         = Genome::Sys->create_temp_file('tmp_indel.'.$input_var_format);
        ($tmp_flt_indel_fh, $tmp_flt_indel_file) = Genome::Sys->create_temp_file('tmp_indel_filtered.'.$input_var_format);
    }

    my $header;
    my $tmp_fh = Genome::Sys->open_file_for_reading($tmp_file) or return;
    while (my $line = $tmp_fh->getline) {
        if ($line =~ /^#/) { # store the header info, vcf has header while pileup does not
            $header .= $line;
            next;
        }
        $self->_print_var($tmp_snv_fh, $tmp_indel_fh, $line);
    }
    $tmp_fh->close;
    map{$_->close if $_}($tmp_snv_fh, $tmp_indel_fh);

    my $flt_tmp_fh = Genome::Sys->open_file_for_reading($flt_tmp_file) or return;
    while (my $flt_line = $flt_tmp_fh->getline) {
        $flt_line =~ s/^\S+\s+//;   #The filtered var file has the first column as a letter
        $self->_print_var($tmp_flt_snv_fh, $tmp_flt_indel_fh, $flt_line);
    }
    $flt_tmp_fh->close;
    map{$_->close if $_}($tmp_flt_snv_fh, $tmp_flt_indel_fh);

    if ($header) {
        my ($header_fh, $header_file) = Genome::Sys->create_temp_file('tmp_header_'.$input_var_format);
        $header_fh->print($header);
        $header_fh->close;

        my @touch_files;

        if ($snv_out_file) {
            Genome::Sys->cat(input_files => [$header_file, $tmp_snv_file],     output_file => $snv_out_file);
            Genome::Sys->cat(input_files => [$header_file, $tmp_flt_snv_file], output_file => $flt_snv_out_file);
        }
        if ($indel_out_file) { #sometimes for testing data, pass_indel_filter or fail_indel_filter file could be empty
            if (-z $tmp_indel_file) {
                $self->warning_message('No indel calls passed the filter');
                push @touch_files, $indel_out_file;
            }
            else {
                Genome::Sys->cat(input_files => [$header_file, $tmp_indel_file], output_file => $indel_out_file);
            }
            if (-z $tmp_flt_indel_file) {
                $self->warning_message('No indel calls failed the filter');
                push @touch_files, $flt_indel_out_file;
            }
            else {
                Genome::Sys->cat(input_files => [$header_file, $tmp_flt_indel_file], output_file => $flt_indel_out_file);
            }
        }
        `touch @touch_files` if @touch_files;
    }
    else {
        if ($snv_out_file) {
            Genome::Sys->copy_file($tmp_snv_file, $snv_out_file);
            Genome::Sys->copy_file($tmp_flt_snv_file, $flt_snv_out_file);
        }
        if ($indel_out_file) {
            Genome::Sys->copy_file($tmp_indel_file, $indel_out_file);
            Genome::Sys->copy_file($tmp_flt_indel_file, $flt_indel_out_file);
        }
    }
    return 1;
}
    

sub _get_cmd_opts {
    my ($self, %options) = @_;
    my $cmd_opt;

    for my $option (keys %options) {
        if (defined $self->$option) {
            my ($opt) = $options{$option}->{doc} =~ /^(\-\S)\s/;
            $opt = $opt_pu_conv{$option} if $self->input_var_format eq 'pileup' and $opt_pu_conv{$option};
            $cmd_opt .= " $opt " . $self->$option;
        }
    }
    return $cmd_opt . ' -p'; #output filtered list
}


sub _print_var {
    my ($self, $snv_fh, $indel_fh, $line) = @_;
    my ($id) = $line =~ /^\S+\s+\S+\s+(\S+)\s+/ if $self->input_var_format eq 'pileup';

    if ($id and $id eq '*') {
        $indel_fh->print($line) if $indel_fh;
    }
    elsif ($line =~ /\sINDEL;/) {
        $indel_fh->print($line) if $indel_fh;
    }
    else {
        $snv_fh->print($line) if $snv_fh;
    }

    return 1;
}
1;
