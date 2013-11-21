# review tmooney
# @chrs is set to parameter value, but later overwritten by (1..22, 'X')?  Then later (1..22, 'X') hardcoded into loop.

package Genome::Model::Tools::Somatic::BroadBamToCna;

## This script analyzes a pair of tumor and normal bam files
## and outputs chromosomal based copy number alteration

use strict;
use warnings;
use Genome;
use Statistics::Descriptive;
use Statistics::R;
require Genome::Sys;

my $DEFAULT_VERSION = '0.1';
my $BAMWINDOW_COMMAND = 'bam-window';

class Genome::Model::Tools::Somatic::BroadBamToCna {
    is => 'Command',
    has => [
    use_version => {
        is => 'Version',
        is_optional => 1,
        default_value => $DEFAULT_VERSION,
        doc => "Version of bam-window to use, default is $DEFAULT_VERSION"
    },
    tumor_bam_file => {
        type => 'String',
        is_input => 1,
        is_optional => 0,
        doc => 'Location of tumor bam file.'
    },
    normal_bam_file => {
        type => 'String',
        is_input => 1,
        is_optional => 0,
        doc => 'Location of normal bam file.'
    },
    output_file => {
        type => 'String',
        is_input => 1,
        is_output => 1,
        is_optional => 0,
        doc => 'Copy number analysis output file (full path).'
    },
    window_size => {
        type => 'Number',
        is_optional => 1,
        default => 10000,
        doc => 'Window size (bp) for counting reads contributing to copy number in that window (resolution, default = 10000 bp).'
    },
    maq_quality_cutoff => {
        type => 'Number',
        is_optional => 1,
        default => 35,
        doc => 'MAQ mapping quality cutoff for contributing reads (default = 35).'
    },
    ratio => {
        type => 'Number',
        is_optional => 1,
        default => 0.25,
        doc => 'Ratio diverged from median, used to find copy number neutral region (default = 0.25).'
    },
    chromosome_list => {
        type => 'String',
        is_optional => 1,
        default => '',
        doc => 'List of chromosomes (comma separated) to use for calculation of median coverage (default = all).'
    },
    tumor_downsample_percentage => {
        type => 'Number',
        is_optional => 1,
        default => 1,
        doc => 'Percent of reads (value x 100%) to use in calculations (max & default = 1).'
    },
    normal_downsample_percentage => {
        type => 'Number',
        is_optional => 1,
        default => 1,
        doc => 'Percent of reads (value x 100%) to use in calculations (max,default = 1).'
    },
    plot => {
        type => 'Boolean',
        is_optional => 1,
        doc => "whether or not to run R plot command at end to create .png image of data. Use --noplot to skip plot. Default is to make a plot.",
        default => 1,
    },
    skip_if_output_present => {
        is => 'Boolean',
        is_optional => 1,
        is_input => 1,
        default => 0,
        doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
    },
    # Make workflow choose 64 bit blades
    lsf_resource => {
        is_param => 1,
        default_value => 'rusage[mem=4000] select[type==LINUX64] span[hosts=1]',
    },
    lsf_queue => {
        is_param => 1,
        default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
    },
    ]
};

my %BAMWINDOW_VERSIONS = (
    '0.1' => $ENV{GENOME_SW} . '/bamwindow/bamwindow-v0.1/' . $BAMWINDOW_COMMAND,
);

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt somatic bam-to-cna --tumor-bam-file tumor.bam --normal-bam-file normal.bam --output-file copy_number.out
gmt somatic bam-to-cna --tumor-bam-file tumor.bam --normal-bam-file normal.bam --output-file copy_number.out --window-size 100000 --ratio 0.3 --maq-quality-cutoff 40 
EOS
}

sub help_brief {
    "This tool analyzes a tumor and normal bam file and outputs chromosomal-based copy number alteration."
}

sub help_detail {
    return<<EOS
This tool analyzes a tumor and normal bam file and outputs chromosomal-based copy number alteration in the form of an output file which displays \"chromosome, position, tumor cn, normal cn, difference\". Also, the tool plots a grid of copy number graphs, one for each chromosome, on a single .png image the size of one full page using an embedded R script.
EOS
}

sub execute {
    my $self = shift;

    my %maps = (tumor => $self->tumor_bam_file, normal => $self->normal_bam_file);
    my @samples = ("tumor","normal");
    my %downratios = (tumor => $self->tumor_downsample_percentage, normal => $self->normal_downsample_percentage);
    my $outfile = $self->output_file;

    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    #test architecture to make sure bam-window program can run (req. 64-bit)
    unless (`uname -a` =~ /x86_64/) {
        $self->error_message("Must run on a 64 bit machine");
        die;
    }

    ####################### Compute read counts in sliding windows ########################
    my %data;
    my %statistics;

    for my $sample (@samples){
        my $cmd = sprintf("%s -w %d -q %d -s -p -d %f %s |", $self->bamwindow_path, $self->window_size, $self->maq_quality_cutoff, $downratios{$sample}, $maps{$sample});
        open(MAP,$cmd) || die "unable to open $maps{$sample}\n";
        $statistics{$sample} = Statistics::Descriptive::Sparse->new();

        my ($previous_chr,$window)=(0,0);
        while(<MAP>){
            chomp;
            my ($chr,$pos,$nread) = split /\t/;
            $window = 0 if($chr ne $previous_chr);
            $data{$sample}{$chr}[$window++]=$nread;
            $previous_chr = $chr;
            $statistics{$sample}->add_data($nread);
        }
        close(MAP);
    }

    my @chrs;
    unless($self->chromosome_list) {
        for my $chr (1..22) {
            my $broad_chr = "chr" . $chr;
            push @chrs,$broad_chr;
        }
        push @chrs,"chrX";
    }
    #@chrs=(chr1..22,'X');
    else {
        @chrs = split /,/, $self->chromosome_list;
    }

    #Estimate genome-wide tumor/normal 2X read count
    my %medians;
    for my $sample (@samples){
        my $median=Statistics::Descriptive::Full->new();
        foreach my $chr(@chrs){
            next unless (defined $data{$sample}{$chr});
            my $md = $self->get_median($data{$sample}{$chr});
            $median->add_data($md);
        }
        $medians{$sample} = $median->median();
    }

    @chrs = ();
    for my $chr (1..22) {
        my $broad_chr = "chr" . $chr;
        push @chrs,$broad_chr;
    }
    push @chrs,"chrX";
    #@chrs = (1..22,'X');

    my %num_CN_neutral_pos;
    my %NReads_CN_neutral;
    foreach my $chr(@chrs){
        next unless (defined $data{tumor}{$chr} && defined $data{normal}{$chr});
        my $tumor_window_count = $#{$data{tumor}{$chr}};
        my $normal_window_count = $#{$data{normal}{$chr}};
        my $window_count = ($tumor_window_count < $normal_window_count)? $tumor_window_count : $normal_window_count; #mininum

        for my $window (0..$window_count){
            my $f2x=1;
            next unless (defined $data{tumor}{$chr}[$window] && defined $data{normal}{$chr}[$window]);
            for my $sample (@samples){
                $f2x=0 if($data{$sample}{$chr}[$window]<$medians{$sample}*(1-$self->ratio) || $data{$sample}{$chr}[$window]>$medians{$sample}*(1+$self->ratio));
            }
            next if(! $f2x);
            $num_CN_neutral_pos{$chr}++;
            $num_CN_neutral_pos{'allchr'}++;
            for my $sample (@samples){
                $NReads_CN_neutral{$chr}{$sample}+=$data{$sample}{$chr}[$window];
                $NReads_CN_neutral{'allchr'}{$sample}+=$data{$sample}{$chr}[$window];
            }
        }
    }

    #subtract the normal from the tumor
    open(OUT, ">$outfile") || die "Unable to open output file $outfile: $!";
    my %depth2x;
    foreach my $chr(@chrs,'allchr'){
        printf OUT "#Chr%s median read count",$chr;
        for my $sample (@samples){
            if($num_CN_neutral_pos{$chr} && $num_CN_neutral_pos{$chr} > 10){
                $depth2x{$chr}{$sample}=$NReads_CN_neutral{$chr}{$sample}/$num_CN_neutral_pos{$chr};
            }
            else{  # Sample < 10, backoff to genome-wide estimation
                $depth2x{$chr}{$sample}=$NReads_CN_neutral{'allchr'}{$sample}/$num_CN_neutral_pos{'allchr'};
            }
            printf OUT "\t%s\:%d",$sample,$depth2x{$chr}{$sample};
        }
        print OUT "\n";
    }
    print OUT "CHR\tPOS\tTUMOR\tNORMAL\tDIFF\n";

    my @included_chrs = ();
    #for my $chr (1..22,'X'){
    for my $chr (@chrs){
        next unless (defined $data{tumor}{$chr} && defined $data{normal}{$chr});
        if($NReads_CN_neutral{$chr}{tumor} and $NReads_CN_neutral{$chr}{normal}) {
            push @included_chrs, $chr;
        } else {
            $self->warning_message('No reads within specified ratio of median for chromosome ' . $chr . '. Skipping. ');
            next;
        }
        my $cov_ratio=$NReads_CN_neutral{$chr}{tumor}/$NReads_CN_neutral{$chr}{normal};

        my $tumor_window_count = $#{$data{tumor}{$chr}};
        my $normal_window_count = $#{$data{normal}{$chr}};
        my $window_count = ($tumor_window_count < $normal_window_count)? $tumor_window_count : $normal_window_count; #mininum

        for my $window (0..$window_count){
            next unless (defined $data{tumor}{$chr}[$window] && defined $data{normal}{$chr}[$window]);
            my $cna_unadjusted=($data{tumor}{$chr}[$window]-$cov_ratio*$data{normal}{$chr}[$window])*2/$depth2x{$chr}{tumor};
            my $poschr=$window*$self->window_size;
            #printf OUT "%s\t%d\t%d\t%d\t%.6f\n",$chr,$poschr,${$data{tumor}{$chr}}[$window]*2/$depth2x{$chr}{normal},${$data{normal}{$chr}}[$window]*2/$depth2x{$chr}{normal},$diff_copy;
            printf OUT "%s\t%d\t%d\t%d\t%.6f\n",$chr,$poschr,$data{tumor}{$chr}[$window],$data{normal}{$chr}[$window],$cna_unadjusted;
        }
    }
    close(OUT);

    #clear some memory
    undef %data;
    undef %statistics;
    undef %medians;
    undef %num_CN_neutral_pos;
    undef %NReads_CN_neutral;
    undef %depth2x;

    #plot output
    if ($self->plot) { 
        $self->plot_output($outfile, \@included_chrs);
    }

    return 1;
}

sub get_median {
    my $self = shift;
    my $rpole = shift;
    my @pole = @$rpole;
    my $ret;

    @pole= sort {$a<=>$b} @pole;
    if( (@pole % 2) == 1 ) {
        $ret = $pole[((@pole+1) / 2)-1];
    } else {
        $ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
    }
    return $ret;
}

sub plot_output {
    use Cwd qw(abs_path cwd);
    my $self = shift;
    my $datafile = shift;
    my $chr_array = shift;

    $datafile = abs_path($datafile);
    my $Routfile = $datafile.".png";
    my $tempdir = Genome::Sys->create_temp_directory();
    my $chr_list = join(',', map("'$_'", @$chr_array));

    #R automatically sets the working directory to its tmp_dir, which prevents Genome::Sys from cleaning it up...
    #So save the original beforehand and restore it after we're done
    my $cwd = cwd();

    my $R = Statistics::R->new(tmp_dir => $tempdir);
    $R->startR();
    $R->send(qq{
        bitmap('$Routfile', height = 8.5, width=11, res=200);
        par(mfrow=c(4,6));
        x=read.table('$datafile',comment.char='#',header=TRUE);
        for (i in c($chr_list)) { y=subset(x,CHR==i); plot( y\$POS, y\$DIFF, main=paste(i), xlab='mb', ylab='cn', type='p', col=rgb(0,0,0), pch='.', ylim=c(-4,4) ) };
        dev.off();
    });
    $R->stopR();

    chdir $cwd; 
}

sub bamwindow_path {
    my $self = $_[0];
    return $self->path_for_bamwindow_version($self->use_version);
}

sub available_bamwindow_versions {
    my $self = shift;
    return keys %BAMWINDOW_VERSIONS;
}

sub path_for_bamwindow_version {
    my $class = shift;
    my $version = shift;

    if (defined $BAMWINDOW_VERSIONS{$version}) {
        return $BAMWINDOW_VERSIONS{$version};
    }
    die('No path for bam-window version '. $version);
}

sub default_bamwindow_version {
    die "default bam-window version: $DEFAULT_VERSION is not valid" unless $BAMWINDOW_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

1;
