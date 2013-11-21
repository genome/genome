package Genome::Model::Tools::Somatic::HighConfidence;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Readonly;
use Genome::Info::IUB;

my $DEFAULT_VERSION = '0.2';
my $READCOUNT_COMMAND = 'bam-readcount';

class Genome::Model::Tools::Somatic::HighConfidence {
    is => 'Command',
    has => [
        bam_readcount_version => {
            is => 'Version',
            is_input=>1, 
            is_optional => 1,
            default_value => $DEFAULT_VERSION,
            doc => "Version of bam-readcount to use, default is $DEFAULT_VERSION"
        },
        bam_readcount_params => {
            is => 'String',
            is_optional => 1,
            is_input=>1, 
            default_value => "-q 1",
            doc => "Parameters to pass to bam-readcount"
        },
       'sniper_file' => {
           type => 'String',
           is_input => 1,
           doc => 'somatic sniper or pileup output',
       },
       'output_file' => {
           type => 'String',
           is_input => 1,
           is_output => 1,
           doc => 'File name in which to write output',
       },
       'tumor_bam_file' => {
            type => 'String',
            doc => 'Tumor bam file in which to examine reads',
            is_input => 1,
       },
       'min_mapping_quality' => {
            type => 'String',
            default => '70',
            is_optional => 1,
            is_input => 1,
            doc => 'minimum average mapping quality threshold for high confidence call',
       },
       'min_somatic_quality' => {
            type => 'String',
            default => '40',
            is_optional => 1,
            is_input => 1,
            doc => 'minimum somatic quality threshold for high confidence call',
       },
       prepend_chr => {
           is => 'Boolean',
           default => '0',
           is_optional => 1,
           is_input => 1,
           doc => 'prepend the string "chr" to chromosome names. This is primarily used for external/imported bam files.',
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
       skip => {
           is => 'Boolean',
           default => '0',
           is_input => 1,
           is_optional => 1,
           doc => "If set to true... this will do nothing! Fairly useless, except this is necessary for workflow.",
       },
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
        },
    ]
};

my %READCOUNT_VERSIONS = (
    '0.2' => $ENV{GENOME_SW} . '/samtools/readcount/readcount-v0.2/' . $READCOUNT_COMMAND,
);

sub help_brief {
    return "This module takes in somatic sniper output and filters it to high confidence variants";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    gmt somatic high-confidence --sniper-file sniper,out --tumor-bam somefile.bam 
EOS
}

sub help_detail {                           
    return <<EOS 
This module takes in somatic sniper output and filters it to high confidence variants
EOS
}

sub execute {
    my $self = shift;

    if ($self->skip) {
        $self->status_message("Skipping execution: Skip flag set");
        return 1;
    }
    
    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }

    #test architecture to make sure we can run read count program
    unless (`uname -a` =~ /x86_64/) {
       $self->error_message("Must run on a 64 bit machine");
       die;
    }

    #check on BAM file
    unless(-e $self->tumor_bam_file) {
        $self->error_message("Tumor bam file: " . $self->tumor_bam_file . " does not exist");
        die;
    }

    unless(-e $self->tumor_bam_file . ".bai") {
        $self->error_message("Tumor bam must be indexed");
        die;
    }
    
    my $fh = IO::File->new($self->sniper_file, "r");
    unless($fh) {
        $self->error_message("Unable to open " . $self->sniper_file . ". $!");
        die;
    }
    
    my $ofh = IO::File->new($self->output_file, "w");
    unless($ofh) {
        $self->error_message("Unable to open " . $self->output_file . " for writing. $!");
        die;
    }

    my ($tfh,$temp_path) = Genome::Sys->create_temp_file;
    unless($tfh) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }
    $temp_path =~ s/\:/\\\:/g;
   
    #read sniper and skip indels
    my $somatic_threshold = $self->min_somatic_quality;
    
    my @sniper_lines;
    while(my $line = $fh->getline) {
        chomp $line;
        my ($chr, $start, $stop, $ref, $iub, $type, $somatic_score, @annotation_columns) = split /\t/, $line;
        if ($self->prepend_chr) {
            $chr = "chr$chr";
            $chr =~ s/MT$/M/;
        };
        next if $ref eq "*";
        if($somatic_score >= $somatic_threshold) {
            print $tfh "$chr\t$start\t$stop\n";
            push @sniper_lines, $line
        }
    }
    $tfh->close;
    
    #Nothing to do if no lines were not "*"
    return 1 unless scalar @sniper_lines;

    #Run readcount program 
    my $readcount_command = sprintf("%s %s -l %s %s |",$self->readcount_path, $self->bam_readcount_params, $temp_path, $self->tumor_bam_file);
    $self->status_message("Running: $readcount_command");
    my $readcounts = IO::File->new($readcount_command);

    while(my $count_line = $readcounts->getline) {
        chomp $count_line;
        my ($chr, $pos, $ref, $depth, @base_stats) = split /\t/, $count_line;
        
        my $current_variant = shift @sniper_lines;
        last unless $current_variant;
        my ($vchr, $vstart, $vstop, $vref, $viub) = split /\t/, $current_variant;

        #check if the sniper line was present in the readcount output
        while($vchr ne $chr && $vstart != $pos && @sniper_lines) {
            $self->status_message("Skipped $current_variant");
            
            $current_variant = shift @sniper_lines;
            last unless $current_variant;
            ($vchr, $vstart, $vstop, $vref, $viub) = split /\t/, $current_variant;
        }
        last unless $current_variant;
        
        my %bases;
        for my $base_stat (@base_stats) {
            my ($base,$reads,$avg_mq, $avg_bq) = split /:/, $base_stat;
            #Leaving bases coded as '=' unhandled
            next if($base eq "=");
            $bases{$base} = $avg_mq;
        }

        my @vars = Genome::Info::IUB->variant_alleles_for_iub($vref,$viub);
        foreach my $var (@vars) {
            if(exists($bases{$var}) && $bases{$var} >= $self->min_mapping_quality) {
                print $ofh $current_variant, "\n";
                last;
            }
        }
    }

    unless($readcounts->close()) {
        $self->error_message("Error running " . $self->readcount_path);
        die;
    }
    
    return 1;
}

sub readcount_path {
    my $self = $_[0];
    return $self->path_for_readcount_version($self->bam_readcount_version);
}

sub available_readcount_versions {
    my $self = shift;
    return keys %READCOUNT_VERSIONS;
}

sub path_for_readcount_version {
    my $class = shift;
    my $version = shift;

    if (defined $READCOUNT_VERSIONS{$version}) {
        return $READCOUNT_VERSIONS{$version};
    }
    die('No path for bam-readcount version '. $version);
}

sub default_readcount_version {
    die "default bam-readcount version: $DEFAULT_VERSION is not valid" unless $READCOUNT_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

1;
