package Genome::Model::Tools::Sv::Pairoscope;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Which;

class Genome::Model::Tools::Sv::Pairoscope {
    is => 'Command',
    has => [
        input_file => 
        { 
            type => 'String',
            is_optional => 0,
            doc => "Input file of svs in primer design input format",
        },
        output_dir =>
        {
            type => 'String',
            is_optional => 0,
            doc => "Output directory name for placement of directories",
        },        
        tumor_bam =>
        {
            type => 'String',
            is_optional => 0,
            doc => "bam file location for tumor",
        },
        normal_bam =>
        {
            type => 'String',
            is_optional => 0,
            doc => "bam file location for normal",
        },
        types => {
            type => 'String',
            is_optional => 1,
            doc => "Comma separated string of types to graph",
            default => "INV,INS,DEL,ITX,CTX",
        },
        possible_BD_type => {
            type => 'hashref',
            doc => "hashref of possible BreakDancer SV types",
            is_optional => 1,
            default => {INV => 1,INS => 1,DEL => 1,ITX => 1,CTX => 1,},
        },
        pairoscope_program => {
            type => "String",
            default => "pairoscope0.3",
            doc => "executable of pairoscope to use", 
            is_optional => 1,
        },
        exon_bam => {
            type => "String",
            example_values => ["/gscmnt/sata135/info/medseq/dlarson/hg18.NCBI-human.combined-annotation.54_36p_v2.sorted.bam"],
            doc => "bam file of exons to use for displaying gene models", 
            is_optional => 1,
        },
        buffer_size => {
            type => "Integer",
            default => 500,
            doc => "number of bases to include on either side of the predicted breakpoint(s)",
            is_optional => 1,
        },
        pairoscope_options => {
            type => "String",
            default => "",
            doc => "option string to pass through to pairoscope for experimental options etc",
            is_optional => 1,
        },
        output_prefix => {
            type => "String",
            default => "",
            doc => "String that will be prepended to all output file names. Cannot contain an underscore or period",
            is_optional => 1,
        },
        simplified_plot => {
            type => "Boolean",
            default => 0,
            doc => "If specified, produces a simplified plot showing only the read connections (middle section)",
            is_optional => 1,
        },

        simplified_plot_title => {
            type => "String",
            doc => "If simplified plot is being produced, adds this title",
            default => "",
            is_optional => 1,
        },


        # output_width => {
        #     type => "Integer",
        #     default => 1024,
        #     doc => "width in pixels of the output document"
        #     is_optional => 1,
        # }
        # output_height => {
        #     type => "Integer",
        #     default => 768,
        #     doc => "height in pixels of the output document"
        #     is_optional => 1,
        # }

    ],
};


sub execute {
    my $self=shift;

    #test architecture to make sure we can run pairoscope program
    #copied from G::M::T::Maq""Align.t 
    unless (`uname -a` =~ /x86_64/) {
        $self->error_message(`uname -a`); #FIXME remove
        $self->error_message("Must run on a 64 bit machine");
        die;
    }
    #Not allowed to store hash in UR?

    my @types = map { uc $_ } split /,/, $self->types;
    my $allowed_types = $self->possible_BD_type;
    foreach my $type (@types) {
        unless(exists($allowed_types->{$type})) {
            $self->error_message("$type type is not a valid BreakDancer SV type");
            return;
        }
    }
    my %types = map {$_ => 1} @types; #create types hash


    unless(-f $self->input_file) {
        $self->error_message("primer design file is not a file: " . $self->input_file);
        return;
    }

    my $indel_fh = IO::File->new($self->input_file);
    unless($indel_fh) {
        $self->error_message("Failed to open filehandle for: " .  $self->input_file );
        return;
    }

    #TODO These should all get checked somehow
    my $tumor_bam = $self->tumor_bam;
    unless(-e $tumor_bam) {
        $self->error_message("$tumor_bam does not exist");
        return;
    }

    my $normal_bam = $self->normal_bam;
    unless(-e $normal_bam) {
        $self->error_message("$normal_bam does not exist");
        return;
    }

    my $output_dir = $self->output_dir;
    unless(-d $output_dir) {
        $self->error_message("$output_dir does not exist");
        return;
    }

    my $grapher = which($self->pairoscope_program);
    unless(defined $grapher) {
        #perhaps the user gave us an absolute path
        $grapher = $self->pairoscope_program;
    }
    unless(-e $grapher && -x $grapher) {
        $self->error_message("$grapher does not exists or is not an executable");
        return;
    }

    my $simple_plot = $self->simplified_plot;
    my $simple_plot_title = $self->simplified_plot_title;

    my $additional_opts = $self->pairoscope_options;

    my $exon_file = $self->exon_bam;
    if($exon_file) {
        unless(-e $exon_file) {
            $self->error_message("$exon_file does not exist");
            return;
        }
        $additional_opts = $additional_opts ? $additional_opts . " -g $exon_file" : "-g $exon_file";
    }
    
    my $prefix = $self->output_prefix;
    if($prefix =~ /[_.]/) {
        $self->error_message("Output prefix cannot contain an underscore or period");
        return;
    }
    $prefix .= '.' if $prefix ne q{};

    $self->debug_message("Using option string $additional_opts");
    my $buffer = $self->buffer_size;

    my $count = 0;
    #assuming we are reasonably sorted
    while ( my $line = $indel_fh->getline) {
        chomp $line;
        $line =~ s/"//g; #kill any quotes that may have snuck in
        my @fields = split /\s+/, $line; #assuming tab separated
        my ($chr1,
            $chr1_pos,
            $chr2,
            $chr2_pos,
            $type,
        ); 
        if($fields[0] =~ /\./) {
            #probably is HQfiltered input
            $self->debug_message("First column contains a period. Assuming Ken Chen's HQFiltered file format.");
            ($chr1,
                $chr1_pos,
                $chr2,
                $chr2_pos,
                $type,
            ) = @fields[1,2,4,6,7];
        }
        else {
            ($chr1,
                $chr1_pos,
                $chr2,
                $chr2_pos,
                $type,
            ) = @fields[0,1,3,4,6];
        }


        #skip headers
        next if $line =~ /^#|START|TYPE/i;
        #validate columns
        unless($chr1 =~ /^[0-9XYNMT_]+$/i) {
            $self->error_message("First chromosome name $chr1 is invalid at line " . $indel_fh->input_line_number);
            return;
        }
        unless($chr2 =~ /^[0-9XYNMT_]+$/i) {
            $self->error_message("Second chromosome name $chr2 is invalid at line " . $indel_fh->input_line_number);
            return;
        }
        unless($chr1_pos =~ /^\d+$/) {
            $self->error_message("First breakpoint coordinate contains non-digit characters $chr1_pos at line " . $indel_fh->input_line_number);
            return;
        }
        unless($chr2_pos =~ /^\d+$/) {
            $self->error_message("Second breakpoint coordinate contains non-digit characters $chr2_pos at line " . $indel_fh->input_line_number);
            return;
        }
        if(exists($types{$type})) {
            $count++;
            #then we should graph it
            #Doing this based on chromosomes in case types ever change
            if($chr1 eq $chr2) {
                my $name = "$output_dir/${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Tumor_${type}.q1.png";
                my $cmd = "$grapher -q 1 -b $buffer -o $name $additional_opts $tumor_bam $chr1 $chr1_pos $chr2_pos";
                $self->debug_message("Running: $cmd");
                system($cmd);
                $name = "$output_dir/${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Normal_${type}.q1.png";
                $cmd = "$grapher -q 1 -b $buffer  -o $name $additional_opts $normal_bam $chr1 $chr1_pos $chr2_pos";
                $self->debug_message("Running: $cmd");
                system($cmd);
                $name = "$output_dir/${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Tumor_${type}.q0.png";
                $cmd = "$grapher -q 0 -b $buffer  -o $name $additional_opts $tumor_bam $chr1 $chr1_pos $chr2_pos";
                $self->debug_message("Running: $cmd");
                system($cmd);
                $name = "$output_dir/${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Normal_${type}.q0.png";
                $cmd = "$grapher -q 0 -b $buffer  -o $name $additional_opts $normal_bam $chr1 $chr1_pos $chr2_pos";
                $self->debug_message("Running: $cmd");
                system($cmd);
            }
            else {
                my $name = "$output_dir/${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Tumor_${type}.q1.png";
                my $cmd = "$grapher -q 1 -b $buffer -o $name $additional_opts $tumor_bam $chr1 $chr1_pos $chr1_pos $tumor_bam $chr2 $chr2_pos $chr2_pos";
                $self->debug_message("Running: $cmd");
                system($cmd);
                $name = "$output_dir/${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Normal_${type}.q1.png";
                $cmd = "$grapher -q 1 -b $buffer  -o $name $additional_opts $normal_bam $chr1 $chr1_pos $chr1_pos $normal_bam $chr2 $chr2_pos $chr2_pos";
                $self->debug_message("Running: $cmd");
                system($cmd);
                $name = "$output_dir/${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Tumor_${type}.q0.png";
                $cmd = "$grapher -q 0 -b $buffer  -o $name $additional_opts $tumor_bam $chr1 $chr1_pos $chr1_pos $tumor_bam $chr2 $chr2_pos $chr2_pos";
                $self->debug_message("Running: $cmd");
                system($cmd);
                $name = "$output_dir/${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Normal_${type}.q0.png";
                $cmd = "$grapher -q 0 -b $buffer  -o $name $additional_opts $normal_bam $chr1 $chr1_pos $chr1_pos $normal_bam $chr2 $chr2_pos $chr2_pos";
                $self->debug_message("Running: $cmd");
                system($cmd);
            }
        }
        unless(exists($allowed_types->{$type})) {
            $self->error_message("Type $type invalid");
            $self->error_message("Valid types are " . join("\t",keys %{$allowed_types}));
            return;
        }

        if($simple_plot){
            #two temp files, one for each plot
            my ($tfh,$newfile1) = Genome::Sys->create_temp_file;
            unless($tfh) {
                $self->error_message("Unable to create temporary file $!");
                die;
            }
            my ($tfh2,$newfile2) = Genome::Sys->create_temp_file;
            unless($tfh2) {
                $self->error_message("Unable to create temporary file $!");
                die;
            }
            my $cmd;
            #use imagemagick to create the plot
            $cmd = "convert -crop 1024x260+0+250 ${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Normal_${type}.q1.png $newfile1"; 
            system($cmd);
            $cmd = "convert -crop 1024x260+0+250 ${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Tumor_${type}.q1.png $newfile2";
            system($cmd);

            # fonts won't work unless you run this script and save the output in ~/.magick/type.xml
            # http://www.imagemagick.org/Usage/scripts/imagick_type_gen
            if($simple_plot_title eq ""){                
                $cmd = "montage -tile 1x2 -geometry 1024x260 -shadow $newfile1 $newfile2 ${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Both_${type}.q1.simple.png";
                system($cmd);
            } else {
                $cmd = "montage -tile 1x2 -geometry 1024x260 -font Kayrawan -title " . $simple_plot_title . " -shadow $newfile1 $newfile2 ${prefix}${chr1}_${chr1_pos}_${chr2}_${chr2_pos}_Both_${type}.q1.simple.png";
                system($cmd);
            }
        }


    }

    $indel_fh->close; 

    return 1;
}

1;

sub help_detail {
    my $help = <<HELP;
Ken Chen's BreakDancer predicts large structural variations by examining read pairs. This module uses the pairoscope program to graph read pairs for a given set of regions. pairoscope operates by scanning a maq map file for reads in the regions and matches up pairs across those regions. The output consists of a set of tracks for each region. One track is the read depth across the region (excluding gapped reads) the other is a so called barcode output. For multiple regions, the regions are displayed in order listed in the filename. Read depth tracks first, then the barcode graphs. Reads are represented as lines and pairs are joined by arcs. These are color coded by abnormal read pair type as follows:

Mapping status                                      Color
Forward-Reverse, abnormal insert size               magenta
Forward-Forward                                     red
Reverse-Reverse                                     blue
Reverse-Forward                                     green
One read unmapped                                   yellow
One read mapped to a different chromosome           cyan

Pairoscope.pm generates 4 PNG images for each predicted SV, 2 for tumor and 2 for normal. There is a q0 file showing reads of all mapping qualities and a q1 file showing reads of mapping quality 1 or more. A maq mapping quality of zero indicates a repeat region that mapped multiple places in the genome equally well.

The naming convention of the files produced is as follows:
(prefix.)chr_pos_chr_pos_tumor/normal_type.q#.png

The input file must be formatted as follows:
chr1	position1	?	chr2	position2	?	type

or

ID	CHR1	OUTER_START	INNER_START	CHR2	INNER_END	OUTER_END	TYPE

Please note that the second format is identified by the presence of a period in the id. If this assumption is wrong for some reason then this script may break.

Problems or questions? Email dlarson\@genome.wustl.edu

HELP

}

sub help_brief {
    return "This module takes a sv primer design file and uses the rudimentary graphical tool pairoscope to graph the read pairs.";
}


