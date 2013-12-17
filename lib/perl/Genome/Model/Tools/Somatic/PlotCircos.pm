#review tmooney
#What was added for BRC1 only?

package Genome::Model::Tools::Somatic::PlotCircos;

use warnings;
use strict;

use Genome;
use Carp;
use IO::File;
use Genome::Sys;
use Cwd qw( abs_path getcwd );
use File::Basename;

class Genome::Model::Tools::Somatic::PlotCircos{
    is => 'Command',
    has_optional => [
        cna_file => {
            is  => 'String',
            is_input  => 1,
            doc => 'Windowed copy number from sequence data (e.g., as output from gmt somatic bam-to-cna)',
        },
        tier1_lc_file => {
            is  => 'String',
            doc => 'The list of tier1 variants',
        },
        tier1_hc_file => {
            is  => 'String',
            is_input  => 1,
            doc => 'The list of tier1 variants in SOMATIC format if you not need labels',
        },
         tier1_hclabel_file => {
            is  => 'String',
            is_input => 1,
            doc => 'The list of tier1 variants in ANNOTATION format if you want them labelled. somatic format BAD',
        },
        tier1_rchclabel_file => {
            is  => 'String',
            doc => 'The list of tier1 recurrent variants in ANNOTATION format if you want them labelled. somatic format BAD',
        },
        sv_file => {
            is  => 'String',
            is_input  => 1,
            doc => 'The list of sv in BREAKDANCER output format... define this OR all of ctx, itx, del, ins, inv files... but not both',
        },
        ctx_file => {
            is  => 'String',
            doc => 'The file containing sv - trans-chromosomal translocations',
        },
        itx_file => {
            is  => 'String',
            doc => 'The file containing sv - intra-chromosomal translocations',
        },
        del_file => {
            is  => 'String',
            doc => 'The file containing sv - deletions, specify this OR sv_file but not both',
        },
        ins_file => {
            is  => 'String',
            doc => 'The file containing sv - insertions, specify this OR sv_file but not both',
        },
        inv_file => {
            is  => 'String',
            doc => 'The file containing sv - inversions, specify this OR sv_file but not both',
        },
        output_file => {
            is => 'String',
            is_input => 1,
            is_output => 1,
            doc=> 'The output png. Provide a base file name, .png will be appended',
        },
        skip_if_output_present => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
        },
        minimum_score_graphed => {
            is => 'Integer',
            is_optional => 1,
            default => 45,
            doc => 'Minimum score to graph... all below this score are ignored. Default is 45.',
        },
        config_file => {
            is => 'String',
            doc => 'circos config file',
        },
        output_dir => {
            is => 'String',
            doc => 'output to dump circos files instead of temp',
        },
        #private variables
        _cna_circos_file => {
            is => 'String',
        },
        _ideogram_file => {
            is => 'String',
        },
        _colors_file => {
            is => 'String',
        },
        _circos_config_dir => {
            is => 'String',
        },
        _tier1_hc_circos_file => {
            is  => 'String',
        },
         _tier1_hclabel_circos_file => {
            is  => 'String',
        },
         _tier1_rchclabel_circos_file => {
            is  => 'String',
        },
        _tier1_lc_circos_file => {
            is  => 'String',
        },
        _ctx_circos_file => {
            is  => 'String',
        },
        _itx_circos_file => {
            is  => 'String',
        },
        _ins_circos_file => {
            is  => 'String',
        },
        _inv_circos_file => {
            is  => 'String',
        },
        _del_circos_file => {
            is  => 'String',
        },
        lsf_resource => {
            is_param => 1,
            default_value => 'rusage[mem=8000, tmp=2000] select[type==LINUX64 && mem > 8000 && tmp > 2000] span[hosts=1] -M 8000000'
        }, 
        lsf_queue => {
            is_param => 1,
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
    ],
};

sub help_brief {
    "make circos graph",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt somatic plot-circos...    
EOS
}

sub help_detail {                           
    return <<EOS 
makes a circos graph 
EOS
}

sub execute {
    my $self = shift;
        
    #test files 

    # Check/setup and split the 1 sv_file into multiple sv files by category if necessary
    if( (defined $self->sv_file) &&
           ( (defined $self->ctx_file) || (defined $self->itx_file) || (defined $self->ins_file) || (defined $self->del_file) || (defined $self->inv_file) ) ) {
        $self->error_message("Either only sv_file OR some/all of (ctx_file, itx_file, ins_file, del_file, inv_file) must be provided... but not both");
        die;
    }

    # Append .png to the output name... grab the base first
    my $circos_output_base = $self->output_file;
    $self->output_file($self->output_file . ".png");

    if (($self->skip_if_output_present)&&(-s $self->output_file)) {
        $self->status_message("Skipping execution: Output is already present and skip_if_output_present is set to true");
        return 1;
    }
    if (defined $self->sv_file) {
        $self->split_sv_input;
    }
    
    my $cna_fh;
    
    if($self->cna_file) {
        $cna_fh = IO::File->new($self->cna_file,"r");

        unless($cna_fh) {
            $self->error_message("Couldn't open " . $self->cna_file);
            return;
        }
    }

    #copynumber information
    my ($cna_circos_fh, $cna_temp_file) = Genome::Sys->create_temp_file();
    $self->_cna_circos_file($cna_temp_file);
    #First convert everything to Circos format.

    #write out converted cna file
    $self->convert_cna_file($cna_fh,$cna_circos_fh);
    $cna_fh->close if($cna_fh);
    $cna_circos_fh->close;

    #convert and write out the tier1 hc sites
    my $hc_fh;
    if($self->tier1_hc_file) {
        $hc_fh = IO::File->new($self->tier1_hc_file,"r");
        unless($hc_fh) {
            $self->error_message("Couldn't open " . $self->tier1_hc_file);
            return;
        }
    }
    my ($hc_circos_fh, $hc_temp_file) = Genome::Sys->create_temp_file();
    $self->_tier1_hc_circos_file($hc_temp_file);
    
    $self->convert_sniper_file($hc_fh, $hc_circos_fh);
    $hc_fh->close if $hc_fh;
    $hc_circos_fh->close;
    
    my $hclabel_fh;
    if($self->tier1_hclabel_file) {
        $hclabel_fh = IO::File->new($self->tier1_hclabel_file,"r");
        unless($hclabel_fh) {
            $self->error_message("Couldn't open " . $self->tier1_hclabel_file);
            return;
        }
    }
    my ($hclabel_circos_fh, $hclabel_temp_file) = Genome::Sys->create_temp_file();
    $self->_tier1_hclabel_circos_file($hclabel_temp_file);
    
    $self->convert_anno_file($hclabel_fh, $hclabel_circos_fh);
    $hclabel_fh->close if $hclabel_fh;
    $hclabel_circos_fh->close;

    #for labeling recurrent Tier1s special!!
    my $rchclabel_fh;
    if($self->tier1_rchclabel_file) {
        $rchclabel_fh = IO::File->new($self->tier1_rchclabel_file,"r");
        unless($rchclabel_fh) {
            $self->error_message("Couldn't open " . $self->tier1_rchclabel_file);
            return;
        }
    }
    my ($rchclabel_circos_fh, $rchclabel_temp_file) = Genome::Sys->create_temp_file();
    $self->_tier1_rchclabel_circos_file($rchclabel_temp_file);
    
    $self->convert_anno_file($rchclabel_fh, $rchclabel_circos_fh);
    $rchclabel_fh->close if $rchclabel_fh;
    $rchclabel_circos_fh->close;
    #now do lc tier1 mutations
    my $lc_fh;
    if($self->tier1_lc_file) {
        $lc_fh = IO::File->new($self->tier1_lc_file,"r");
        unless($lc_fh) {
            $self->error_message("Couldn't open " . $self->tier1_lc_file);
            return;
        }
    }
    my ($lc_circos_fh, $lc_temp_file) = Genome::Sys->create_temp_file();
    $self->_tier1_lc_circos_file($lc_temp_file);
    
    $self->convert_sniper_file($lc_fh, $lc_circos_fh);
    $lc_fh->close if $lc_fh;
    $lc_circos_fh->close;

    #now do SV
    my $ctx_fh;
    if($self->ctx_file) {
        $ctx_fh = IO::File->new($self->ctx_file,"r");
        unless($ctx_fh) {
            $self->error_message("Couldn't open " . $self->ctx_file);
            return;
        }
    }
    my ($ctx_circos_fh, $ctx_temp_file) = Genome::Sys->create_temp_file();
    $self->_ctx_circos_file($ctx_temp_file);
    
    $self->convert_breakdancer_file($ctx_fh, $ctx_circos_fh, "dgreen2");
    $ctx_fh->close if $ctx_fh;
    $ctx_circos_fh->close;
    my $itx_fh;
    if($self->itx_file) {
        $itx_fh = IO::File->new($self->itx_file,"r");
        unless($itx_fh) {
            $self->error_message("Couldn't open " . $self->itx_file);
            return;
        }
    }
    my ($itx_circos_fh, $itx_temp_file) = Genome::Sys->create_temp_file();
    $self->_itx_circos_file($itx_temp_file);
    
    $self->convert_breakdancer_file($itx_fh, $itx_circos_fh, "lgreen2");
    $itx_fh->close if $itx_fh;
    $itx_circos_fh->close;
    
    my $ins_fh;
    if($self->ins_file) {
        $ins_fh = IO::File->new($self->ins_file,"r");
        unless($ins_fh) {
            $self->error_message("Couldn't open " . $self->ins_file);
            return;
        }
    }
    my ($ins_circos_fh, $ins_temp_file) = Genome::Sys->create_temp_file();
    $self->_ins_circos_file($ins_temp_file);
    
    $self->convert_breakdancer_file($ins_fh, $ins_circos_fh, "orange2");
    $ins_fh->close if $ins_fh;
    $ins_circos_fh->close;
    
    my $inv_fh;
    if($self->inv_file) {
        $inv_fh = IO::File->new($self->inv_file,"r");
        unless($inv_fh) {
            $self->error_message("Couldn't open " . $self->inv_file);
            return;
        }
    }
    my ($inv_circos_fh, $inv_temp_file) = Genome::Sys->create_temp_file();
    $self->_inv_circos_file($inv_temp_file);
    
    $self->convert_breakdancer_file($inv_fh, $inv_circos_fh, "purple2");
    $inv_fh->close if $inv_fh;
    $inv_circos_fh->close;

    my $del_fh;
    if($self->del_file) {
        $del_fh = IO::File->new($self->del_file,"r");
        unless($del_fh) {
            $self->error_message("Couldn't open " . $self->del_file);
            return;
        }
    }
    my ($del_circos_fh, $del_temp_file) = Genome::Sys->create_temp_file();
    $self->_del_circos_file($del_temp_file);
    
    $self->convert_breakdancer_file($del_fh, $del_circos_fh, 'dblue2');
    $del_fh->close if $del_fh;
    $del_circos_fh->close;
    
    ####ADDED FOR BRC1 ONLY####
    
    #write out config files etc 
    my ($ideogram_fh, $ideogram_path) = Genome::Sys->create_temp_file($self->ideogram_file_name);
    print $ideogram_fh $self->ideogram_file_contents;
    $ideogram_fh->close;
    $self->_ideogram_file($ideogram_path);

    my ($color_fh, $color_path) = Genome::Sys->create_temp_file($self->color_file_name);
    print $color_fh $self->color_file_contents;
    $color_fh->close;
    $self->_colors_file($color_path);

    my ($config_fh, $config_path) = Genome::Sys->create_temp_file($self->config_file_name);
    print $config_fh $self->config_file_contents;
    $config_fh->close;

    
    print `circos -conf $config_path`;
   
    my $circos_big = $self->output_file;
    my $circos_smallest = "$circos_output_base.920x920.png";
    my $circos_small = "$circos_output_base.3000x3000.png";
    `convert $circos_big -resize 3000x3000 -interpolate bicubic -quality 1000 $circos_small`;
    `convert $circos_small -resize 920x920 -interpolate bicubic -quality 1000 $circos_smallest`;
    
    #Then graph. Done!
    return 1;
}

sub convert_sniper_file {
    my ($self,$sniper_fh, $output_fh) = @_;
    if($sniper_fh) { 
        unless($sniper_fh->opened && $output_fh->opened) {
            return;
        }

        my $label = 0; #this will give each SNP a unique label
        while(my $line = $sniper_fh->getline) {
            $label++;
            chomp $line;
            my ($chr,$start, @rest) = split /\t/, $line;

            print $output_fh "hs$chr $start $start\n";
        }
    }
    return 1;
}

sub convert_anno_file {
    my ($self,$sniper_fh, $output_fh) = @_;
    if($sniper_fh) { 
        unless($sniper_fh->opened && $output_fh->opened) {
            return;
        }

        my $label = 0; #this will give each SNP a unique label
        while(my $line = $sniper_fh->getline) {
            next if $line =~ /gene_name/;
            $label++;
            chomp $line;
            my ($chr,$start, $stop, $ref, $var, $mut, $gene, $transcript, $species, $source, $version, $strand, $dunno, $type, $c_position, $amino_acid_change, @rest) = split /\t/, $line;
            my $label = '';
            $amino_acid_change =~ s/p\.//;
            $label = "$gene\[$amino_acid_change\]" if $gene;
            print $output_fh "hs$chr $start $start $label\n";
        }
    }
    return 1;
}

sub split_sv_input {
    my $self = shift;

    my $sv_fh = IO::File->new($self->sv_file,"r");
    unless($sv_fh) {
        $self->error_message("Couldn't open " . $self->sv_file);
        return;
    }

    # For each of the 5 SV types calculate a name for that file and assign it to the object level accessor... 
    # Also make filehandles for each so we can split input into each of these files (saves us from doing all this 5 times)
    my %fh_by_type;
    for my $type ('ctx', 'itx', 'ins', 'del', 'inv') {
        my $property_name = $type . "_file";
        my $file = $self->sv_file . ".$type";
        $self->$property_name($file);
        my $fh = IO::File->new($self->$property_name,"w");
        unless($fh) {
            $self->error_message("Couldn't open $file");
            return;
        }
        $fh_by_type{$type} = $fh;
    }

    while (my $line = $sv_fh->getline) {
        # Filter out the top 3 "comments" / headers in breakdancer files which have # at the beginning
        next if $line =~ m/^#/;
        my ($chr1, $pos1, $orientation1, $chr2, $pos2, $orientation2, $type, $size, $score, $num_reads, $num_reads_lib, $allele_frequency, $version, $run_param) = split("\t", $line);

        my $fh = $fh_by_type{lc($type)};
        unless (defined $fh) {
            $self->error_message("Error in split_sv_input: SV type $type is unrecognized... should be CTX, ITX, INS, DEL, INS");
            die;
        }

        $fh->print($line);
    }

    $sv_fh->close;
    for my $type ('ctx', 'itx', 'ins', 'del', 'inv') {
        my $fh = $fh_by_type{$type};
        $fh->close;
    }
    
    return 1;
}

sub convert_breakdancer_file {
    my ($self, $breakdancer_fh, $output_fh, $color) = @_;
    if($breakdancer_fh) { 
        unless($breakdancer_fh->opened && $output_fh->opened) {
            return;
        }
        my $color_label;
        my $label = 0; #this will give each SV a unique label
        while(my $line = $breakdancer_fh->getline) {
            $label++;
            chomp $line;
            my ($chr1,$breakpoint1,$orientation1,$chr2,$breakpoint2,$orientation2,$type,$size, $score, $reads,)= split /\t/, $line;
            unless($score > $self->minimum_score_graphed) {
                next;
            }
            # (score range - cutoff) / 100 = divisor
            my $bin;

            print  "$type$label\ths$chr1\t$breakpoint1\t$breakpoint1\tcolor=$color\t$score\n";

            print $output_fh "$type$label\ths$chr1\t$breakpoint1\t$breakpoint1\tcolor=$color\n";
            print $output_fh "$type$label\ths$chr2\t$breakpoint2\t$breakpoint2\tcolor=$color\n";
        }
    }
    return 1;
}

sub convert_cna_file {
    my ($self, $map2cna_fh, $output_fh) = @_;
    if($map2cna_fh) {
        unless($map2cna_fh->opened && $output_fh->opened) {
            return;
        }

        while(my $line = $map2cna_fh->getline) {
            chomp $line;
            next if $line =~ /^#|^CHR/i; #ignore comment lines and header line
            my @fields = split /\t/, $line;
            printf $output_fh "hs%s\t%d\t%d\t%f\n",$fields[0], $fields[1],$fields[1],$fields[-1]; #chr\tstart\tend\tdifference between tumor and normal\n
        }
    }
    return 1;
}

sub ideogram_file_name {
    return "pipeline_ideogram.conf";
}

sub ideogram_file_contents {
    #This just stores the contents of the pipeline's ideogram configuration file with the source. It could be excised and stored in an appropriate
    #central location
    return <<IDEOGRAM;
###
### ideogram.conf
###

<ideogram>

<spacing>

default = 10u

</spacing>

# thickness (px) of chromosome ideogram
thickness        = 300p
stroke_thickness = 2
# ideogram border color
stroke_color     = black2
fill             = yes
# the default chromosome color is set here and any value
# defined in the karyotype file overrides it
fill_color       = black

# fractional radius position of chromosome ideogram within image
radius         = 0.8r
show_label     = yes
label_font     = condensedbold
label_center= yes
label_radius = (dims(ideogram,radius_outer)+dims(ideogram,radius_inner))/2
label_size     = 80 
label_color     = black2
# cytogenetic bands
band_stroke_thickness = 2

# show_bands determines whether the outline of cytogenetic bands
# will be seen
show_bands            = no 
# in order to fill the bands with the color defined in the karyotype
# file you must set fill_bands
fill_bands            = no

band_transparency     = 1

</ideogram>
IDEOGRAM
    
}


sub config_file_name {
    my $self = shift;
    return "pipeline.conf";
}

sub config_file_contents {
    my $self = shift;
    #set up filenames here
    my $ideogram = $self->_ideogram_file;
    my $colors = $self->_colors_file;
    my $cna = $self->_cna_circos_file;
    my $hc = $self->_tier1_hc_circos_file;
    my $lc = $self->_tier1_lc_circos_file;
    my $ctx = $self->_ctx_circos_file;
    my $itx = $self->_itx_circos_file;
    my $ins = $self->_ins_circos_file;
    my $inv = $self->_inv_circos_file;
    my $del = $self->_del_circos_file;
    my $hc_label=$self->_tier1_hclabel_circos_file;
    my $rchc_label=$self->_tier1_rchclabel_circos_file;
    
    my $file = $self->output_file;
    my $dir;
    #Store Cwd if necessary
    if($file) {
        ($file,$dir) = fileparse($file); #circos will change the path of this file
        $dir = abs_path($dir);
    }
    else {
        $file = 'test.png';
        $dir = getcwd();
    }
        

    return <<CONF;
###
### gmt somatic plot-circos configuration file
###

<colors>
<<include $colors>>
</colors>

<fonts>
<<include $ENV{GENOME_SW}/circos/installed/etc/fonts.conf>>
</fonts>

<<include $ideogram>>

karyotype = $ENV{GENOME_SW}/circos/installed/data/7/karyotype.human.colorbychr.txt

<image>

background = white
dir = $dir 
file  = $file 
24bit = yes
#auto_alpha_colors = yes
#auto_alpha_steps  = 10 
png = yes
#svg = yes
# radius of inscribed circle in image
radius         = 5400p
# by default angle=0 is at 3 o'clock position
angle_offset   = -90
#angle_orientation = counterclockwise
</image>

chromosomes_units           = 1000000

chromosomes_display_default = yes

anglestep       = 0.5
minslicestep    = 10
beziersamples   = 40
debug           = no
warnings        = no
imagemap        = no

units_ok = bupr
units_nounit = n

<plots>
<plot>
type             = text
color            = black2
file       = $hc_label

r0 = 1r
r1 = 1.2r

label_snuggle             = yes
# shift label up to its height in pixels in the angular direction
max_snuggle_distance      = 50r
snuggle_sampling          = 2
snuggle_tolerance         = 1r
snuggle_link_overlap_test = yes
snuggle_link_overlap_tolerance = 1p
snuggle_refine            = yes


show_links     = yes
link_dims      = 6p,6p,70p,6p,6p
link_thickness = 9p
link_color     = black2

label_size   = 90p
label_font   = condensed

padding  = 12p
rpadding = 12p

</plot>

<plot>
type             = text
color            = red
file       = $rchc_label

r0 = 1r
r1 = 1.2r

label_snuggle             = yes
# shift label up to its height in pixels in the angular direction
max_snuggle_distance      = 50r
snuggle_sampling          = 2
snuggle_tolerance         = 1r
snuggle_link_overlap_test = yes
snuggle_link_overlap_tolerance = 1p
snuggle_refine            = yes


show_links     = yes
link_dims      = 6p,6p,70p,6p,6p
link_thickness = 9p
link_color     = red

label_size   = 90p
label_font   = condensed

padding  = 12p
rpadding = 12p

</plot>


<plot>

 show  = yes
 type  = scatter
 file  = $cna 
 glyph = circle
 glyph_size = 8
 fill_color = black2
 stroke_color = black2
 stroke_thickness = 1
 min   = -4 
 max   = 4 
 r0    = 0.74r
 r1    = .94r
 
 axis           = yes
 axis_color     = lgrey
 axis_thickness = 2
 axis_spacing   = 0.001
</plot>



</plots>


<highlights>

<highlight>
fill_color = dred 
stroke_color = dred
stroke_thickness = 15
file       = $hc 
r0         = .95r
r1         = .99r
</highlight>

<highlight>
fill_color = lred 
stroke_color = lred
stroke_thickness = 15
file       = $lc 
r0         = .95r
r1         = .99r
</highlight>

</highlights>

<links>

z      = 10
radius = 0.7r
bezier_radius = 0.1r

<link deletions>
color = dblue2
thickness = 15
radius = 0.73r
bezier_radius = .7r
crest = 10 
file = $del 
</link>

<link inversions>
color = purple2
thickness = 15
radius = 0.73r
bezier_radius = .7r
crest = 10
file = $inv 
record_limit = 10000
</link>

<link insertions>
color = orange2 
thickness = 15
radius = 0.73r
bezier_radius = .7r
crest = 10
file = $ins 
record_limit = 100000
</link>

<link translocations>
thickness = 15
radius = 0.73r
bezier_radius = .7r
crest = 1
file = $ctx 
</link>

<link intratranslocations>
color = lgreen2 
thickness = 5
radius = 0.73r
bezier_radius = .7r
crest = 2
file = $itx 
</link>

</links>
CONF
}

sub color_file_name {
    my $self=shift;
    return "colors.conf";
}


sub color_file_contents {
    return <<COLORS;
optblue =  55,133,221
optblue_a1 =  55,133,221,12
optblue_a2 =  55,133,221,24
optblue_a3 =  55,133,221,36
optblue_a4 =  55,133,221,48
optblue_a5 =  55,133,221,60
optblue_a6 =  55,133,221,78
optblue_a7 =  55,133,221,90
optblue_a8 =  55,133,221,102
optblue_a9 =  55,133,221,114
optgreen =  55,221,125
optgreen_a1 =  55,221,125,12
optgreen_a2 =  55,221,125,24
optgreen_a3 =  55,221,125,36
optgreen_a4 =  55,221,125,48
optgreen_a5 =  55,221,125,60
optgreen_a6 =  55,221,125,78
optgreen_a7 =  55,221,125,90
optgreen_a8 =  55,221,125,102
optgreen_a9 =  55,221,125,114
optyellow =  221,215,55
optyellow_a1 =  221,215,55,12
optyellow_a2 =  221,215,55,24
optyellow_a3 =  221,215,55,36
optyellow_a4 =  221,215,55,48
optyellow_a5 =  221,215,55,60
optyellow_a6 =  221,215,55,78
optyellow_a7 =  221,215,55,90
optyellow_a8 =  221,215,55,102
optyellow_a9 =  221,215,55,114
optorange =  221,164,55
optorange_a1 =  221,164,55,12
optorange_a2 =  221,164,55,24
optorange_a3 =  221,164,55,36
optorange_a4 =  221,164,55,48
optorange_a5 =  221,164,55,60
optorange_a6 =  221,164,55,78
optorange_a7 =  221,164,55,90
optorange_a8 =  221,164,55,102
optorange_a9 =  221,164,55,114
optred =  221,55,55
optred_a1 =  221,55,55,12
optred_a2 =  221,55,55,24
optred_a3 =  221,55,55,36
optred_a4 =  221,55,55,48
optred_a5 =  221,55,55,60
optred_a6 =  221,55,55,78
optred_a7 =  221,55,55,90
optred_a8 =  221,55,55,102
optred_a9 =  221,55,55,114
optviolet =  145,55,221
optviolet_a1 =  145,55,221,12
optviolet_a2 =  145,55,221,24
optviolet_a3 =  145,55,221,36
optviolet_a4 =  145,55,221,48
optviolet_a5 =  145,55,221,60
optviolet_a6 =  145,55,221,78
optviolet_a7 =  145,55,221,90
optviolet_a8 =  145,55,221,102
optviolet_a9 =  145,55,221,114
optpurple =  219,55,221,0
optpurple_a1 =  219,55,221,12
optpurple_a2 =  219,55,221,24
optpurple_a3 =  219,55,221,36
optpurple_a4 =  219,55,221,48
optpurple_a5 =  219,55,221,60
optpurple_a6 =  219,55,221,78
optpurple_a7 =  219,55,221,90
optpurple_a8 =  219,55,221,102
optpurple_a9 =  219,55,221,114
white =  255,255,255
white_a1 =  255,255,255,12
white_a2 =  255,255,255,24
white_a3 =  255,255,255,36
white_a4 =  255,255,255,48
white_a5 =  255,255,255,60
white_a6 =  255,255,255,78
white_a7 =  255,255,255,90
white_a8 =  255,255,255,102
white_a9 =  255,255,255,114
vvvvlgrey =  250,250,250
vvvvlgrey_a1 =  250,250,250,12
vvvvlgrey_a2 =  250,250,250,24
vvvvlgrey_a3 =  250,250,250,36
vvvvlgrey_a4 =  250,250,250,48
vvvvlgrey_a5 =  250,250,250,60
vvvvlgrey_a6 =  250,250,250,78
vvvvlgrey_a7 =  250,250,250,90
vvvvlgrey_a8 =  250,250,250,102
vvvvlgrey_a9 =  250,250,250,114
vvvlgrey =  240,240,240
vvvlgrey_a1 =  240,240,240,12
vvvlgrey_a2 =  240,240,240,24
vvvlgrey_a3 =  240,240,240,36
vvvlgrey_a4 =  240,240,240,48
vvvlgrey_a5 =  240,240,240,60
vvvlgrey_a6 =  240,240,240,78
vvvlgrey_a7 =  240,240,240,90
vvvlgrey_a8 =  240,240,240,102
vvvlgrey_a9 =  240,240,240,114
vvlgrey =  230,230,230
vvlgrey_a1 =  230,230,230,12
vvlgrey_a2 =  230,230,230,24
vvlgrey_a3 =  230,230,230,36
vvlgrey_a4 =  230,230,230,48
vvlgrey_a5 =  230,230,230,60
vvlgrey_a6 =  230,230,230,78
vvlgrey_a7 =  230,230,230,90
vvlgrey_a8 =  230,230,230,102
vvlgrey_a9 =  230,230,230,114
vlgrey =  220,220,220
vlgrey_a1 =  220,220,220,12
vlgrey_a2 =  220,220,220,24
vlgrey_a3 =  220,220,220,36
vlgrey_a4 =  220,220,220,48
vlgrey_a5 =  220,220,220,60
vlgrey_a6 =  220,220,220,78
vlgrey_a7 =  220,220,220,90
vlgrey_a8 =  220,220,220,102
vlgrey_a9 =  220,220,220,114
lgrey =  210,210,210
lgrey_a1 =  210,210,210,12
lgrey_a2 =  210,210,210,24
lgrey_a3 =  210,210,210,36
lgrey_a4 =  210,210,210,48
lgrey_a5 =  210,210,210,60
lgrey_a6 =  210,210,210,78
lgrey_a7 =  210,210,210,90
lgrey_a8 =  210,210,210,102
lgrey_a9 =  210,210,210,114
grey =  200,200,200
grey_a1 =  200,200,200,12
grey_a2 =  200,200,200,24
grey_a3 =  200,200,200,36
grey_a4 =  200,200,200,48
grey_a5 =  200,200,200,60
grey_a6 =  200,200,200,78
grey_a7 =  200,200,200,90
grey_a8 =  200,200,200,102
grey_a9 =  200,200,200,114
dgrey =  170,170,170
dgrey_a1 =  170,170,170,12
dgrey_a2 =  170,170,170,24
dgrey_a3 =  170,170,170,36
dgrey_a4 =  170,170,170,48
dgrey_a5 =  170,170,170,60
dgrey_a6 =  170,170,170,78
dgrey_a7 =  170,170,170,90
dgrey_a8 =  170,170,170,102
dgrey_a9 =  170,170,170,114
vdgrey =  140,140,140
vdgrey_a1 =  140,140,140,12
vdgrey_a2 =  140,140,140,24
vdgrey_a3 =  140,140,140,36
vdgrey_a4 =  140,140,140,48
vdgrey_a5 =  140,140,140,60
vdgrey_a6 =  140,140,140,78
vdgrey_a7 =  140,140,140,90
vdgrey_a8 =  140,140,140,102
vdgrey_a9 =  140,140,140,114
vvdgrey =  100,100,100
vvdgrey_a1 =  100,100,100,12
vvdgrey_a2 =  100,100,100,24
vvdgrey_a3 =  100,100,100,36
vvdgrey_a4 =  100,100,100,48
vvdgrey_a5 =  100,100,100,60
vvdgrey_a6 =  100,100,100,78
vvdgrey_a7 =  100,100,100,90
vvdgrey_a8 =  100,100,100,102
vvdgrey_a9 =  100,100,100,114
vvvdgrey =  70,70,70
vvvdgrey_a1 =  70,70,70,12
vvvdgrey_a2 =  70,70,70,24
vvvdgrey_a3 =  70,70,70,36
vvvdgrey_a4 =  70,70,70,48
vvvdgrey_a5 =  70,70,70,60
vvvdgrey_a6 =  70,70,70,78
vvvdgrey_a7 =  70,70,70,90
vvvdgrey_a8 =  70,70,70,102
vvvdgrey_a9 =  70,70,70,114
vvvvdgrey =  40,40,40
vvvvdgrey_a1 =  40,40,40,12
vvvvdgrey_a2 =  40,40,40,24
vvvvdgrey_a3 =  40,40,40,36
vvvvdgrey_a4 =  40,40,40,48
vvvvdgrey_a5 =  40,40,40,60
vvvvdgrey_a6 =  40,40,40,78
vvvvdgrey_a7 =  40,40,40,90
vvvvdgrey_a8 =  40,40,40,102
vvvvdgrey_a9 =  40,40,40,114
black =  0,0,0,0
black2 = 1,1,1,0
black_a1 =  0,0,0,12
black_a2 =  0,0,0,24
black_a3 =  0,0,0,36
black_a4 =  0,0,0,48
black_a5 =  0,0,0,60
black_a6 =  0,0,0,78
black_a7 =  0,0,0,90
black_a8 =  0,0,0,102
black_a9 =  0,0,0,114
vlred =  255,193,200
vlred_a1 =  255,193,200,12
vlred_a2 =  255,193,200,24
vlred_a3 =  255,193,200,36
vlred_a4 =  255,193,200,48
vlred_a5 =  255,193,200,60
vlred_a6 =  255,193,200,78
vlred_a7 =  255,193,200,90
vlred_a8 =  255,193,200,102
vlred_a9 =  255,193,200,114
lred =  255,122,137
lred_a1 =  255,122,137,12
lred_a2 =  255,122,137,24
lred_a3 =  255,122,137,36
lred_a4 =  255,122,137,48
lred_a5 =  255,122,137,60
lred_a6 =  255,122,137,78
lred_a7 =  255,122,137,90
lred_a8 =  255,122,137,102
lred_a9 =  255,122,137,114
red =  247,42,66
red_a1 =  247,42,66,12
red_a2 =  247,42,66,24
red_a3 =  247,42,66,36
red_a4 =  247,42,66,48
red_a5 =  247,42,66,60
red_a6 =  247,42,66,78
red_a7 =  247,42,66,90
red_a8 =  247,42,66,102
red_a9 =  247,42,66,114
dred =  205,51,69
dred_a1 =  205,51,69,12
dred_a2 =  205,51,69,24
dred_a3 =  205,51,69,36
dred_a4 =  205,51,69,48
dred_a5 =  205,51,69,60
dred_a6 =  205,51,69,78
dred_a7 =  205,51,69,90
dred_a8 =  205,51,69,102
dred_a9 =  205,51,69,114
vlgreen =  204,255,218,0
vlgreen_a1 =  204,255,218,12
vlgreen_a2 =  204,255,218,24
vlgreen_a3 =  204,255,218,36
vlgreen_a4 =  204,255,218,48
vlgreen_a5 =  204,255,218,60
vlgreen_a6 =  204,255,218,78
vlgreen_a7 =  204,255,218,90
vlgreen_a8 =  204,255,218,102
vlgreen_a9 =  204,255,218,114
lgreen =  128,255,164,0
lgreen2 = 127,254,163,0
lgreen_a1 =  128,255,164,12
lgreen_a2 =  128,255,164,24
lgreen_a3 =  128,255,164,36
lgreen_a4 =  128,255,164,48
lgreen_a5 =  128,255,164,60
lgreen_a6 =  128,255,164,78
lgreen_a7 =  128,255,164,90
lgreen_a8 =  128,255,164,102
lgreen_a9 =  128,255,164,114
green =  51,204,94,0
green_a1 =  51,204,94,12
green_a2 =  51,204,94,24
green_a3 =  51,204,94,36
green_a4 =  51,204,94,48
green_a5 =  51,204,94,60
green_a6 =  51,204,94,78
green_a7 =  51,204,94,90
green_a8 =  51,204,94,102
green_a9 =  51,204,94,114
dgreen =  38,153,71,0
dgreen2 =  37,152,70,0
dgreen_a1 =  38,153,71,12
dgreen_a2 =  38,153,71,24
dgreen_a3 =  38,153,71,36
dgreen_a4 =  38,153,71,48
dgreen_a5 =  38,153,71,60
dgreen_a6 =  38,153,71,78
dgreen_a7 =  38,153,71,90
dgreen_a8 =  38,153,71,102
dgreen_a9 =  38,153,71,114
vlblue =  128,176,255
vlblue_a1 =  128,176,255,12
vlblue_a2 =  128,176,255,24
vlblue_a3 =  128,176,255,36
vlblue_a4 =  128,176,255,48
vlblue_a5 =  128,176,255,60
vlblue_a6 =  128,176,255,78
vlblue_a7 =  128,176,255,90
vlblue_a8 =  128,176,255,102
vlblue_a9 =  128,176,255,114
lblue =  64,137,255
lblue_a1 =  64,137,255,12
lblue_a2 =  64,137,255,24
lblue_a3 =  64,137,255,36
lblue_a4 =  64,137,255,48
lblue_a5 =  64,137,255,60
lblue_a6 =  64,137,255,78
lblue_a7 =  64,137,255,90
lblue_a8 =  64,137,255,102
lblue_a9 =  64,137,255,114
blue =  54,116,217,0
blue_a1 =  54,116,217,12
blue_a2 =  54,116,217,24
blue_a3 =  54,116,217,36
blue_a4 =  54,116,217,48
blue_a5 =  54,116,217,60
blue_a6 =  54,116,217,78
blue_a7 =  54,116,217,90
blue_a8 =  54,116,217,102
blue_a9 =  54,116,217,114
dblue =  38,82,153,0
dblue2 =  39,83,154,0
dblue_a1 =  38,82,153,12
dblue_a2 =  38,82,153,24
dblue_a3 =  38,82,153,36
dblue_a4 =  38,82,153,48
dblue_a5 =  38,82,153,60
dblue_a6 =  38,82,153,78
dblue_a7 =  38,82,153,90
dblue_a8 =  38,82,153,102
dblue_a9 =  38,82,153,114
lpurple =  236,64,255
lpurple_a1 =  236,64,255,12
lpurple_a2 =  236,64,255,24
lpurple_a3 =  236,64,255,36
lpurple_a4 =  236,64,255,48
lpurple_a5 =  236,64,255,60
lpurple_a6 =  236,64,255,78
lpurple_a7 =  236,64,255,90
lpurple_a8 =  236,64,255,102
lpurple_a9 =  236,64,255,114
purple =  189,51,204
purple2 =  190,52,205,0
purple_a1 =  189,51,204,12
purple_a2 =  189,51,204,24
purple_a3 =  189,51,204,36
purple_a4 =  189,51,204,48
purple_a5 =  189,51,204,60
purple_a6 =  189,51,204,78
purple_a7 =  189,51,204,90
purple_a8 =  189,51,204,102
purple_a9 =  189,51,204,114
dpurple =  118,32,128
dpurple_a1 =  118,32,128,12
dpurple_a2 =  118,32,128,24
dpurple_a3 =  118,32,128,36
dpurple_a4 =  118,32,128,48
dpurple_a5 =  118,32,128,60
dpurple_a6 =  118,32,128,78
dpurple_a7 =  118,32,128,90
dpurple_a8 =  118,32,128,102
dpurple_a9 =  118,32,128,114
vlyellow =  255,253,202
vlyellow_a1 =  255,253,202,12
vlyellow_a2 =  255,253,202,24
vlyellow_a3 =  255,253,202,36
vlyellow_a4 =  255,253,202,48
vlyellow_a5 =  255,253,202,60
vlyellow_a6 =  255,253,202,78
vlyellow_a7 =  255,253,202,90
vlyellow_a8 =  255,253,202,102
vlyellow_a9 =  255,253,202,114
lyellow =  255,252,150
lyellow_a1 =  255,252,150,12
lyellow_a2 =  255,252,150,24
lyellow_a3 =  255,252,150,36
lyellow_a4 =  255,252,150,48
lyellow_a5 =  255,252,150,60
lyellow_a6 =  255,252,150,78
lyellow_a7 =  255,252,150,90
lyellow_a8 =  255,252,150,102
lyellow_a9 =  255,252,150,114
yellow =  255,255,0
yellow_a1 =  255,255,0,12
yellow_a2 =  255,255,0,24
yellow_a3 =  255,255,0,36
yellow_a4 =  255,255,0,48
yellow_a5 =  255,255,0,60
yellow_a6 =  255,255,0,78
yellow_a7 =  255,255,0,90
yellow_a8 =  255,255,0,102
yellow_a9 =  255,255,0,114
dyellow =  191,186,48
dyellow_a1 =  191,186,48,12
dyellow_a2 =  191,186,48,24
dyellow_a3 =  191,186,48,36
dyellow_a4 =  191,186,48,48
dyellow_a5 =  191,186,48,60
dyellow_a6 =  191,186,48,78
dyellow_a7 =  191,186,48,90
dyellow_a8 =  191,186,48,102
dyellow_a9 =  191,186,48,114
lime =  186,255,0
lime_a1 =  186,255,0,12
lime_a2 =  186,255,0,24
lime_a3 =  186,255,0,36
lime_a4 =  186,255,0,48
lime_a5 =  186,255,0,60
lime_a6 =  186,255,0,78
lime_a7 =  186,255,0,90
lime_a8 =  186,255,0,102
lime_a9 =  186,255,0,114
vlorange =  255,228,193
vlorange_a1 =  255,228,193,12
vlorange_a2 =  255,228,193,24
vlorange_a3 =  255,228,193,36
vlorange_a4 =  255,228,193,48
vlorange_a5 =  255,228,193,60
vlorange_a6 =  255,228,193,78
vlorange_a7 =  255,228,193,90
vlorange_a8 =  255,228,193,102
vlorange_a9 =  255,228,193,114
lorange =  255,187,110
lorange_a1 =  255,187,110,12
lorange_a2 =  255,187,110,24
lorange_a3 =  255,187,110,36
lorange_a4 =  255,187,110,48
lorange_a5 =  255,187,110,60
lorange_a6 =  255,187,110,78
lorange_a7 =  255,187,110,90
lorange_a8 =  255,187,110,102
lorange_a9 =  255,187,110,114
orange =  255,136,0,0
orange2 = 254,135,0,0
orange_a1 =  255,136,0,12
orange_a2 =  255,136,0,24
orange_a3 =  255,136,0,36
orange_a4 =  255,136,0,48
orange_a5 =  255,136,0,60
orange_a6 =  255,136,0,78
orange_a7 =  255,136,0,90
orange_a8 =  255,136,0,102
orange_a9 =  255,136,0,114
dorange =  221,143,55
dorange_a1 =  221,143,55,12
dorange_a2 =  221,143,55,24
dorange_a3 =  221,143,55,36
dorange_a4 =  221,143,55,48
dorange_a5 =  221,143,55,60
dorange_a6 =  221,143,55,78
dorange_a7 =  221,143,55,90
dorange_a8 =  221,143,55,102
dorange_a9 =  221,143,55,114
gpos100 =  0,0,0,0
gpos100_a1 =  0,0,0,12
gpos100_a2 =  0,0,0,24
gpos100_a3 =  0,0,0,36
gpos100_a4 =  0,0,0,48
gpos100_a5 =  0,0,0,60
gpos100_a6 =  0,0,0,78
gpos100_a7 =  0,0,0,90
gpos100_a8 =  0,0,0,102
gpos100_a9 =  0,0,0,114
gpos =  0,0,0,0
gpos_a1 =  0,0,0,12
gpos_a2 =  0,0,0,24
gpos_a3 =  0,0,0,36
gpos_a4 =  0,0,0,48
gpos_a5 =  0,0,0,60
gpos_a6 =  0,0,0,78
gpos_a7 =  0,0,0,90
gpos_a8 =  0,0,0,102
gpos_a9 =  0,0,0,114
gpos75 =  130,130,130
gpos75_a1 =  130,130,130,12
gpos75_a2 =  130,130,130,24
gpos75_a3 =  130,130,130,36
gpos75_a4 =  130,130,130,48
gpos75_a5 =  130,130,130,60
gpos75_a6 =  130,130,130,78
gpos75_a7 =  130,130,130,90
gpos75_a8 =  130,130,130,102
gpos75_a9 =  130,130,130,114
gpos66 =  160,160,160
gpos66_a1 =  160,160,160,12
gpos66_a2 =  160,160,160,24
gpos66_a3 =  160,160,160,36
gpos66_a4 =  160,160,160,48
gpos66_a5 =  160,160,160,60
gpos66_a6 =  160,160,160,78
gpos66_a7 =  160,160,160,90
gpos66_a8 =  160,160,160,102
gpos66_a9 =  160,160,160,114
gpos50 =  200,200,200
gpos50_a1 =  200,200,200,12
gpos50_a2 =  200,200,200,24
gpos50_a3 =  200,200,200,36
gpos50_a4 =  200,200,200,48
gpos50_a5 =  200,200,200,60
gpos50_a6 =  200,200,200,78
gpos50_a7 =  200,200,200,90
gpos50_a8 =  200,200,200,102
gpos50_a9 =  200,200,200,114
gpos33 =  210,210,210
gpos33_a1 =  210,210,210,12
gpos33_a2 =  210,210,210,24
gpos33_a3 =  210,210,210,36
gpos33_a4 =  210,210,210,48
gpos33_a5 =  210,210,210,60
gpos33_a6 =  210,210,210,78
gpos33_a7 =  210,210,210,90
gpos33_a8 =  210,210,210,102
gpos33_a9 =  210,210,210,114
gpos25 =  200,200,200
gpos25_a1 =  200,200,200,12
gpos25_a2 =  200,200,200,24
gpos25_a3 =  200,200,200,36
gpos25_a4 =  200,200,200,48
gpos25_a5 =  200,200,200,60
gpos25_a6 =  200,200,200,78
gpos25_a7 =  200,200,200,90
gpos25_a8 =  200,200,200,102
gpos25_a9 =  200,200,200,114
gvar =  220,220,220
gvar_a1 =  220,220,220,12
gvar_a2 =  220,220,220,24
gvar_a3 =  220,220,220,36
gvar_a4 =  220,220,220,48
gvar_a5 =  220,220,220,60
gvar_a6 =  220,220,220,78
gvar_a7 =  220,220,220,90
gvar_a8 =  220,220,220,102
gvar_a9 =  220,220,220,114
gneg =  255,255,255
gneg_a1 =  255,255,255,12
gneg_a2 =  255,255,255,24
gneg_a3 =  255,255,255,36
gneg_a4 =  255,255,255,48
gneg_a5 =  255,255,255,60
gneg_a6 =  255,255,255,78
gneg_a7 =  255,255,255,90
gneg_a8 =  255,255,255,102
gneg_a9 =  255,255,255,114
acen =  217,47,39
acen_a1 =  217,47,39,12
acen_a2 =  217,47,39,24
acen_a3 =  217,47,39,36
acen_a4 =  217,47,39,48
acen_a5 =  217,47,39,60
acen_a6 =  217,47,39,78
acen_a7 =  217,47,39,90
acen_a8 =  217,47,39,102
acen_a9 =  217,47,39,114
stalk =  100,127,164
stalk_a1 =  100,127,164,12
stalk_a2 =  100,127,164,24
stalk_a3 =  100,127,164,36
stalk_a4 =  100,127,164,48
stalk_a5 =  100,127,164,60
stalk_a6 =  100,127,164,78
stalk_a7 =  100,127,164,90
stalk_a8 =  100,127,164,102
stalk_a9 =  100,127,164,114
select =  135,177,255
select_a1 =  135,177,255,12
select_a2 =  135,177,255,24
select_a3 =  135,177,255,36
select_a4 =  135,177,255,48
select_a5 =  135,177,255,60
select_a6 =  135,177,255,78
select_a7 =  135,177,255,90
select_a8 =  135,177,255,102
select_a9 =  135,177,255,114
nyt_blue =  104,152,178
nyt_blue_a1 =  104,152,178,12
nyt_blue_a2 =  104,152,178,24
nyt_blue_a3 =  104,152,178,36
nyt_blue_a4 =  104,152,178,48
nyt_blue_a5 =  104,152,178,60
nyt_blue_a6 =  104,152,178,78
nyt_blue_a7 =  104,152,178,90
nyt_blue_a8 =  104,152,178,102
nyt_blue_a9 =  104,152,178,114
nyt_green =  137,129,96
nyt_green_a1 =  137,129,96,12
nyt_green_a2 =  137,129,96,24
nyt_green_a3 =  137,129,96,36
nyt_green_a4 =  137,129,96,48
nyt_green_a5 =  137,129,96,60
nyt_green_a6 =  137,129,96,78
nyt_green_a7 =  137,129,96,90
nyt_green_a8 =  137,129,96,102
nyt_green_a9 =  137,129,96,114
nyt_yellow =  241,221,117
nyt_yellow_a1 =  241,221,117,12
nyt_yellow_a2 =  241,221,117,24
nyt_yellow_a3 =  241,221,117,36
nyt_yellow_a4 =  241,221,117,48
nyt_yellow_a5 =  241,221,117,60
nyt_yellow_a6 =  241,221,117,78
nyt_yellow_a7 =  241,221,117,90
nyt_yellow_a8 =  241,221,117,102
nyt_yellow_a9 =  241,221,117,114
nyt_orange =  230,146,57
nyt_orange_a1 =  230,146,57,12
nyt_orange_a2 =  230,146,57,24
nyt_orange_a3 =  230,146,57,36
nyt_orange_a4 =  230,146,57,48
nyt_orange_a5 =  230,146,57,60
nyt_orange_a6 =  230,146,57,78
nyt_orange_a7 =  230,146,57,90
nyt_orange_a8 =  230,146,57,102
nyt_orange_a9 =  230,146,57,114
nyt_red =  217,47,39
nyt_red_a1 =  217,47,39,12
nyt_red_a2 =  217,47,39,24
nyt_red_a3 =  217,47,39,36
nyt_red_a4 =  217,47,39,48
nyt_red_a5 =  217,47,39,60
nyt_red_a6 =  217,47,39,78
nyt_red_a7 =  217,47,39,90
nyt_red_a8 =  217,47,39,102
nyt_red_a9 =  217,47,39,114
chr1 =  153,102,0
chr1_a1 =  153,102,0,12
chr1_a2 =  153,102,0,24
chr1_a3 =  153,102,0,36
chr1_a4 =  153,102,0,48
chr1_a5 =  153,102,0,60
chr1_a6 =  153,102,0,78
chr1_a7 =  153,102,0,90
chr1_a8 =  153,102,0,102
chr1_a9 =  153,102,0,114
chr2 =  102,102,0
chr2_a1 =  102,102,0,12
chr2_a2 =  102,102,0,24
chr2_a3 =  102,102,0,36
chr2_a4 =  102,102,0,48
chr2_a5 =  102,102,0,60
chr2_a6 =  102,102,0,78
chr2_a7 =  102,102,0,90
chr2_a8 =  102,102,0,102
chr2_a9 =  102,102,0,114
chr3 =  153,153,30
chr3_a1 =  153,153,30,12
chr3_a2 =  153,153,30,24
chr3_a3 =  153,153,30,36
chr3_a4 =  153,153,30,48
chr3_a5 =  153,153,30,60
chr3_a6 =  153,153,30,78
chr3_a7 =  153,153,30,90
chr3_a8 =  153,153,30,102
chr3_a9 =  153,153,30,114
chr4 =  204,0,0
chr4_a1 =  204,0,0,12
chr4_a2 =  204,0,0,24
chr4_a3 =  204,0,0,36
chr4_a4 =  204,0,0,48
chr4_a5 =  204,0,0,60
chr4_a6 =  204,0,0,78
chr4_a7 =  204,0,0,90
chr4_a8 =  204,0,0,102
chr4_a9 =  204,0,0,114
chr5 =  255,0,0
chr5_a1 =  255,0,0,12
chr5_a2 =  255,0,0,24
chr5_a3 =  255,0,0,36
chr5_a4 =  255,0,0,48
chr5_a5 =  255,0,0,60
chr5_a6 =  255,0,0,78
chr5_a7 =  255,0,0,90
chr5_a8 =  255,0,0,102
chr5_a9 =  255,0,0,114
chr6 =  255,0,204
chr6_a1 =  255,0,204,12
chr6_a2 =  255,0,204,24
chr6_a3 =  255,0,204,36
chr6_a4 =  255,0,204,48
chr6_a5 =  255,0,204,60
chr6_a6 =  255,0,204,78
chr6_a7 =  255,0,204,90
chr6_a8 =  255,0,204,102
chr6_a9 =  255,0,204,114
chr7 =  255,204,204
chr7_a1 =  255,204,204,12
chr7_a2 =  255,204,204,24
chr7_a3 =  255,204,204,36
chr7_a4 =  255,204,204,48
chr7_a5 =  255,204,204,60
chr7_a6 =  255,204,204,78
chr7_a7 =  255,204,204,90
chr7_a8 =  255,204,204,102
chr7_a9 =  255,204,204,114
chr8 =  255,153,0
chr8_a1 =  255,153,0,12
chr8_a2 =  255,153,0,24
chr8_a3 =  255,153,0,36
chr8_a4 =  255,153,0,48
chr8_a5 =  255,153,0,60
chr8_a6 =  255,153,0,78
chr8_a7 =  255,153,0,90
chr8_a8 =  255,153,0,102
chr8_a9 =  255,153,0,114
chr9 =  255,204,0
chr9_a1 =  255,204,0,12
chr9_a2 =  255,204,0,24
chr9_a3 =  255,204,0,36
chr9_a4 =  255,204,0,48
chr9_a5 =  255,204,0,60
chr9_a6 =  255,204,0,78
chr9_a7 =  255,204,0,90
chr9_a8 =  255,204,0,102
chr9_a9 =  255,204,0,114
chr10 =  255,255,0
chr10_a1 =  255,255,0,12
chr10_a2 =  255,255,0,24
chr10_a3 =  255,255,0,36
chr10_a4 =  255,255,0,48
chr10_a5 =  255,255,0,60
chr10_a6 =  255,255,0,78
chr10_a7 =  255,255,0,90
chr10_a8 =  255,255,0,102
chr10_a9 =  255,255,0,114
chr11 =  204,255,0
chr11_a1 =  204,255,0,12
chr11_a2 =  204,255,0,24
chr11_a3 =  204,255,0,36
chr11_a4 =  204,255,0,48
chr11_a5 =  204,255,0,60
chr11_a6 =  204,255,0,78
chr11_a7 =  204,255,0,90
chr11_a8 =  204,255,0,102
chr11_a9 =  204,255,0,114
chr12 =  0,255,0
chr12_a1 =  0,255,0,12
chr12_a2 =  0,255,0,24
chr12_a3 =  0,255,0,36
chr12_a4 =  0,255,0,48
chr12_a5 =  0,255,0,60
chr12_a6 =  0,255,0,78
chr12_a7 =  0,255,0,90
chr12_a8 =  0,255,0,102
chr12_a9 =  0,255,0,114
chr13 =  53,128,0
chr13_a1 =  53,128,0,12
chr13_a2 =  53,128,0,24
chr13_a3 =  53,128,0,36
chr13_a4 =  53,128,0,48
chr13_a5 =  53,128,0,60
chr13_a6 =  53,128,0,78
chr13_a7 =  53,128,0,90
chr13_a8 =  53,128,0,102
chr13_a9 =  53,128,0,114
chr14 =  0,0,204
chr14_a1 =  0,0,204,12
chr14_a2 =  0,0,204,24
chr14_a3 =  0,0,204,36
chr14_a4 =  0,0,204,48
chr14_a5 =  0,0,204,60
chr14_a6 =  0,0,204,78
chr14_a7 =  0,0,204,90
chr14_a8 =  0,0,204,102
chr14_a9 =  0,0,204,114
chr15 =  102,153,255
chr15_a1 =  102,153,255,12
chr15_a2 =  102,153,255,24
chr15_a3 =  102,153,255,36
chr15_a4 =  102,153,255,48
chr15_a5 =  102,153,255,60
chr15_a6 =  102,153,255,78
chr15_a7 =  102,153,255,90
chr15_a8 =  102,153,255,102
chr15_a9 =  102,153,255,114
chr16 =  153,204,255
chr16_a1 =  153,204,255,12
chr16_a2 =  153,204,255,24
chr16_a3 =  153,204,255,36
chr16_a4 =  153,204,255,48
chr16_a5 =  153,204,255,60
chr16_a6 =  153,204,255,78
chr16_a7 =  153,204,255,90
chr16_a8 =  153,204,255,102
chr16_a9 =  153,204,255,114
chr17 =  0,255,255
chr17_a1 =  0,255,255,12
chr17_a2 =  0,255,255,24
chr17_a3 =  0,255,255,36
chr17_a4 =  0,255,255,48
chr17_a5 =  0,255,255,60
chr17_a6 =  0,255,255,78
chr17_a7 =  0,255,255,90
chr17_a8 =  0,255,255,102
chr17_a9 =  0,255,255,114
chr18 =  204,255,255
chr18_a1 =  204,255,255,12
chr18_a2 =  204,255,255,24
chr18_a3 =  204,255,255,36
chr18_a4 =  204,255,255,48
chr18_a5 =  204,255,255,60
chr18_a6 =  204,255,255,78
chr18_a7 =  204,255,255,90
chr18_a8 =  204,255,255,102
chr18_a9 =  204,255,255,114
chr19 =  153,0,204
chr19_a1 =  153,0,204,12
chr19_a2 =  153,0,204,24
chr19_a3 =  153,0,204,36
chr19_a4 =  153,0,204,48
chr19_a5 =  153,0,204,60
chr19_a6 =  153,0,204,78
chr19_a7 =  153,0,204,90
chr19_a8 =  153,0,204,102
chr19_a9 =  153,0,204,114
chr20 =  204,51,255
chr20_a1 =  204,51,255,12
chr20_a2 =  204,51,255,24
chr20_a3 =  204,51,255,36
chr20_a4 =  204,51,255,48
chr20_a5 =  204,51,255,60
chr20_a6 =  204,51,255,78
chr20_a7 =  204,51,255,90
chr20_a8 =  204,51,255,102
chr20_a9 =  204,51,255,114
chr21 =  204,153,255
chr21_a1 =  204,153,255,12
chr21_a2 =  204,153,255,24
chr21_a3 =  204,153,255,36
chr21_a4 =  204,153,255,48
chr21_a5 =  204,153,255,60
chr21_a6 =  204,153,255,78
chr21_a7 =  204,153,255,90
chr21_a8 =  204,153,255,102
chr21_a9 =  204,153,255,114
chr22 =  102,102,102
chr22_a1 =  102,102,102,12
chr22_a2 =  102,102,102,24
chr22_a3 =  102,102,102,36
chr22_a4 =  102,102,102,48
chr22_a5 =  102,102,102,60
chr22_a6 =  102,102,102,78
chr22_a7 =  102,102,102,90
chr22_a8 =  102,102,102,102
chr22_a9 =  102,102,102,114
chr23 =  153,153,153
chr23_a1 =  153,153,153,12
chr23_a2 =  153,153,153,24
chr23_a3 =  153,153,153,36
chr23_a4 =  153,153,153,48
chr23_a5 =  153,153,153,60
chr23_a6 =  153,153,153,78
chr23_a7 =  153,153,153,90
chr23_a8 =  153,153,153,102
chr23_a9 =  153,153,153,114
chrX =  153,153,153
chrX_a1 =  153,153,153,12
chrX_a2 =  153,153,153,24
chrX_a3 =  153,153,153,36
chrX_a4 =  153,153,153,48
chrX_a5 =  153,153,153,60
chrX_a6 =  153,153,153,78
chrX_a7 =  153,153,153,90
chrX_a8 =  153,153,153,102
chrX_a9 =  153,153,153,114
chr24 =  204,204,204
chr24_a1 =  204,204,204,12
chr24_a2 =  204,204,204,24
chr24_a3 =  204,204,204,36
chr24_a4 =  204,204,204,48
chr24_a5 =  204,204,204,60
chr24_a6 =  204,204,204,78
chr24_a7 =  204,204,204,90
chr24_a8 =  204,204,204,102
chr24_a9 =  204,204,204,114
chrY =  204,204,204
chrY_a1 =  204,204,204,12
chrY_a2 =  204,204,204,24
chrY_a3 =  204,204,204,36
chrY_a4 =  204,204,204,48
chrY_a5 =  204,204,204,60
chrY_a6 =  204,204,204,78
chrY_a7 =  204,204,204,90
chrY_a8 =  204,204,204,102
chrY_a9 =  204,204,204,114
chrM =  204,204,153
chrM_a1 =  204,204,153,12
chrM_a2 =  204,204,153,24
chrM_a3 =  204,204,153,36
chrM_a4 =  204,204,153,48
chrM_a5 =  204,204,153,60
chrM_a6 =  204,204,153,78
chrM_a7 =  204,204,153,90
chrM_a8 =  204,204,153,102
chrM_a9 =  204,204,153,114
chr0 =  204,204,153
chr0_a1 =  204,204,153,12
chr0_a2 =  204,204,153,24
chr0_a3 =  204,204,153,36
chr0_a4 =  204,204,153,48
chr0_a5 =  204,204,153,60
chr0_a6 =  204,204,153,78
chr0_a7 =  204,204,153,90
chr0_a8 =  204,204,153,102
chr0_a9 =  204,204,153,114
chrUn =  121,204,61
chrUn_a1 =  121,204,61,12
chrUn_a2 =  121,204,61,24
chrUn_a3 =  121,204,61,36
chrUn_a4 =  121,204,61,48
chrUn_a5 =  121,204,61,60
chrUn_a6 =  121,204,61,78
chrUn_a7 =  121,204,61,90
chrUn_a8 =  121,204,61,102
chrUn_a9 =  121,204,61,114
chrNA =  255,255,255
chrNA_a1 =  255,255,255,12
chrNA_a2 =  255,255,255,24
chrNA_a3 =  255,255,255,36
chrNA_a4 =  255,255,255,48
chrNA_a5 =  255,255,255,60
chrNA_a6 =  255,255,255,78
chrNA_a7 =  255,255,255,90
chrNA_a8 =  255,255,255,102
chrNA_a9 =  255,255,255,114
COLORS
}
