#extract reads from a set of bams in breakdancer predicted regions
package Genome::Model::Tools::Sv::AssemblyValidation;


use strict;
use warnings;
use Genome;
use Bio::Seq;
use Bio::SeqIO;
use File::Temp;
use File::Basename;

=cut
my %opts = (l=>500,p=>1000,s=>0,q=>0,n=>0,m=>10,x=>3,P=>-10,G=>-10,S=>0.02,A=>500,Q=>0);
getopts('l:d:c:p:r:s:Q:q:t:n:a:b:f:km:MRzv:hD:x:i:P:G:I:A:S:L:',\%opts);
die("
Usage:   AssemblyValidation.pl <SV file, default in BreakDancer format> <bam files ... >
Options:
         -d DIR     Directory where all intermediate files are saved
         -I DIR     Read intermediate files from DIR instead of creating them
         -f FILE    Save Breakpoint sequences in file
         -r FILE    Save relevant cross_match alignment results to file
         -z         Customized SV format, interpret column from header that must start with # and contain chr1, start, chr2, end, type, size
         -v FILE    Save unconfirmed predictions in FILE
         -h         Make homo variants het by adding same amount of randomly selected wiletype reads
         -l INT     Flanking size [$opts{l}]
         -A INT     Esimated maximal insert size [$opts{A}]
         -q INT     Only assemble reads with mapping quality > [$opts{q}]
         -m INT     Minimal size (bp) for an assembled SV to be called confirmed [$opts{m}]
         -i INT     invalidate indels are -i bp bigger or smaller than the predicted size, usually 1 std insert size
         -a INT     Get reads with start position bp into the left breakpoint, default [50,100,150]
         -b INT     Get reads with start position bp into the right breakpoint, default [50,100,150]
         -S FLOAT   Maximally allowed polymorphism rate in the flanking region [$opts{S}];
         -P INT     Substitution penalty in cross_match alignment [$opts{P}]
         -G INT     Gap Initialization penalty in cross_match alignment [$opts{G}]
         -M         Prefix reference name with \'chr\' when fetching reads using samtools view
         -R         Assemble Mouse calls, NCBI reference build 37

Filtering:
         -p INT     Ignore cases that have average read depth greater than [$opts{p}]
         -c STRING  Specify a single chromosome
         -s INT     Minimal size of the region to analysis [$opts{s}]
         -Q INT     minimal BreakDancer score required for analysis [$opts{Q}]
         -t STRING  type of SV
         -n INT     minimal number of supporting reads [$opts{n}]
         -L STRING  Ingore calls supported by libraries that contains (comma separated) STRING
         -k         Attach assembly results as additional columns in the input file
         -D DIR     A directory that contains a set of supplementary reads (when they are missing from the main data stream)
         -x INT     Duplicate supplementary reads [$opts{x}] times
\n") unless ($#ARGV>=1);
=cut


class Genome::Model::Tools::Sv::AssemblyValidation {
    is  => 'Genome::Model::Tools::Sv',
    has => [
        sv_file => {
            type => 'String',
            doc  => 'input SV file, default in BreakDancer format',
        },
        bam_files => {
            type => 'String',
            doc  => 'bam files, comma separated',
        },
        output_file  => {
            type => 'String',
            doc  => 'output file of assembly validation',
        },
    ],
    has_optional => [
        tigra_sv_output_description_file => {
            is => "String",
            doc => "Optional path to save description of tigra-sv output",
        },
        intermediate_save_dir => {
            type => 'String',
            doc  => 'Directory where all intermediate files are saved',
        },
        intermediate_read_dir => {
            type => 'String',
            doc  => 'Read intermediate files from DIR instead of creating them',
        },
        breakpoint_seq_file => {
            type => 'String',
            doc  => 'Save Breakpoint sequences in file',
        },
        specify_chr => {
            type => 'String',
            doc => 'Specify a single chromosome',
        },
        cm_aln_file => {
            type => 'String',
            doc  => 'Save relevant cross_match alignment results to file',
        },
        custom_sv_format => {
            type => 'Boolean',
            doc  => 'Customized SV format, interpret column from header that must start with # and contain chr1, start, chr2, end, type, size',
        },
        asm_high_coverage => {
            type => 'Boolean',
            doc  => 'Assemble high coverage data (such as capture)',
        },
        unconfirm_predict_file => {
            type => 'String',
            doc  => 'Save unconfirmed predictions in file',
        },
        pad_local_ref => {
            type => 'Integer',
            doc  => 'Pad local reference by additional bp on both ends',
            default_value => 200,
        },
        mismatch_limit => {
            type => 'Integer',
            doc  => 'Number of mismatches required to be tagged as poorly mapped',
            default_value => 5,
        },
        flank_size => {
            type => 'Integer',
            doc  => 'Flanking size',
            default_value => 500,
        },
        est_max_ins_size => {
            type => 'Integer',
            doc  => 'Esimated maximal insert size',
            default_value => 500,
        },
        map_qual_to_asm => {
            type => 'Integer',
            doc  => 'Only assemble reads with mapping quality >',
            default_value => 0,
        },
        min_size_of_confirm_asm_sv => {
            type => 'Integer',
            doc  => 'Minimal size (bp) for an assembled SV to be called confirmed',
            default_value => 3,
        },
        invalid_indel_range => {
            type => 'Integer',
            doc  => 'invalidate indels are -i bp bigger or smaller than the predicted size, usually 1 std insert size',
        },
        max_polymorphism_rate => {
            type => 'Number',
            doc  => 'Maximally allowed polymorphism rate in the flanking region',
            default_value => 0.02,
        },
        cm_sub_penalty => {
            type => 'Number',
            doc  => 'Substitution penalty in cross_match alignment',
            default_value => -10,
        },
        cm_gap_init_penalty => {
            type => 'Number',
            doc  => 'Gap Initialization penalty in cross_match alignment',
            default_value => -10,
        },
        assemble_mouse => {
            type => 'Boolean',
            doc  => 'Assemble Mouse calls, NCBI reference build 37',
        },
        avg_read_depth_limit => {
            type => 'Integer',
            doc  => 'Ignore cases that have average read depth greater than',
            default_value => 1000,
        },
        min_breakdancer_score => {
            type => 'Number',
            doc  => 'minimal BreakDancer score required for analysis',
            default_value => 0,
        },
        skip_libraries => {
            type => 'String',
            doc  => 'Ingore calls supported by libraries that contains (comma separated)',
        },
        skip_call => {
            type => 'Integer',
            doc => 'Ignore tigra_sv for the calls before this number',
        },
        reference_file => {
            type => 'String',
            doc => 'Specify reference in particular',
        },
        _unc_pred_fh => {
            type => 'SCALAR',
        },
        _data_dir => {
            type => 'SCALAR',
        },
        _maxSV => {
            type => 'HASH',
        },
        _N50size => {
            type => 'Number',
        },
        _WeightAvgSize => {
            type => 'Number',
        },
        _tigra_sv_output_map => {
            is => "HASH",
            is_transient => 1,
        }
    ],
};


       
sub execute {
    my $self    = shift;
    my $sv_file = $self->sv_file;
    my $out_file= $self->output_file;

    my $bp_file = $self->breakpoint_seq_file;
    my $bp_io;
    my $cm_aln_fh = new IO::File;
    my $unc_pred_fh;

    if ($bp_file) {
        if (-s $bp_file) {
            $self->warning_message('breakpoint seq file: '.$bp_file. ' existing, Now remove it');
            unlink $bp_file;
        }
        $bp_io = Bio::SeqIO->new(
            -file   => ">>$bp_file",
            -format => 'Fasta',
        );
    }

    my $cm_aln_file = $self->cm_aln_file;
    if ($cm_aln_file) {
        #$cm_aln_fh = Genome::Sys->open_file_for_writing($cm_aln_file) or return;
        $cm_aln_fh->open(">$cm_aln_file");
    }

    my $unc_pred_file = $self->unconfirm_predict_file;
    if ($unc_pred_file) {
        #my $unc_pred_fh = Genome::Sys->open_file_for_writing($unc_pred_file) or return;
        $unc_pred_fh->open(">$unc_pred_fh");
        $self->_unc_pred_fh($unc_pred_fh);
    }

    #my $out_fh  = Genome::Sys->open_file_for_writing($out_file) or return;
    #$out_fh->print("\#CHR1\tPOS1\tCHR2\tPOS2\tORI\tSIZE\tTYPE\tHET\twASMSCORE\tTRIMMED_CONTIG_SIZE\tALIGNED\%\tNUM_SEG\tNUM_FSUB\tNUM_FINDEL\tBP_FINDEL\tMicroHomology\tMicroInsertion\tPREFIX\tASMPARM\tCopyNumber\tGene\tKnown\n");

    srand(time ^ $$);

    my $datadir = $self->intermediate_read_dir;
    my $sv_base = basename $self->sv_file;
    unless ($datadir) {
        if (defined $self->intermediate_save_dir) {
            $datadir = "/tmp/$sv_base";
            if (-d $datadir) {
                $self->error_message("Datadir $datadir already existing\n");
                return;
            }
            mkdir $datadir;
        }
        else {
            $datadir = File::Temp::tempdir("SV_Assembly_XXXXXX", DIR => '/tmp', CLEANUP => 1);
        }
    }
    $self->debug_message("Data directory: $datadir");
    $self->_data_dir($datadir);

    my $tigra_sv_cmd = Genome::Sys->swpath('tigra-sv', '0.1');
    my $tigra_sv_options = $self->_get_tigra_options;
    my $bam_files = $self->_check_bam;
    my $tigra_sv_desc_file = $self->tigra_sv_output_description_file;
    $tigra_sv_desc_file = "$datadir/output_description.tsv" unless $tigra_sv_desc_file;
    $tigra_sv_cmd .= " -D $tigra_sv_desc_file ";
    $tigra_sv_cmd .= ' '. $tigra_sv_options . $sv_file . $bam_files . " > " . $out_file;

    print "tigra-sv command: $tigra_sv_cmd\n";
    $self->debug_message("tigra-sv command: $tigra_sv_cmd");
    my $rv = Genome::Sys->shellcmd(
        cmd           => $tigra_sv_cmd,
        input_files   => [$sv_file],
        #output_files => [$snv_output_file],
        #allow_zero_size_output_files => 1,
    );
    unless ($rv) {
        $self->error_message("Running tigra_sv failed.\nCommand: $tigra_sv_cmd");
        return;
    }

    my @tigra_sv_fas = glob("$datadir/*.fa.contigs.fa");#get tigra homo ctg list
    @tigra_sv_fas = sort{(basename ($a)=~/^\S+?\.(\d+)\./)[0]<=> (basename ($b)=~/^\S+?\.(\d+)\./)[0]}@tigra_sv_fas; #sort ctg file by chr pos
    
    
    # open output file for the following resume
    my $out_fh = new IO::File;
    $out_fh->open(">>$out_file");

    my $output_map = Genome::Model::Tools::TigraSv->parse_output_description($tigra_sv_desc_file);
    $self->_tigra_sv_output_map($output_map);

    for my $tigra_sv_fa (@tigra_sv_fas) {
        my ($tigra_sv_filename) = basename $tigra_sv_fa;
        my ($tigra_sv_name) = $tigra_sv_filename =~ /^(\S+)\.fa\.contigs\.fa/;
        my ($chr1,$start,$chr2,$end,$type,$size,$ori,undef) = split /\./, $tigra_sv_name; # you get the $size from $prefix        
        if (exists $output_map->{$tigra_sv_filename}) {
            my $p = $output_map->{$tigra_sv_filename};
            $chr1 = $p->{chr1};
            $start = $p->{start};
            $chr2 = $p->{chr2};
            $end = $p->{end};
            $type = $p->{type};
            $size = $p->{size};
            $ori = $p->{ori};
        }

        if(defined $self->specify_chr){
            if($chr2 ne $self->specify_chr){
                next;
            }
        }
        my $prefix = join('.',$chr1,$start,$chr2,$end,$type,$size,$ori);

        $self->_N50size(_ComputeTigraN50($tigra_sv_fa));
        $self->_WeightAvgSize(_ComputeTigraWeightedAvgSize($tigra_sv_fa));

        
        #test homo, het contigs
        for my $ctg_type ('homo', 'het') {
            $self->_cross_match_validation($ctg_type, $tigra_sv_name, $tigra_sv_filename);
        }
        
        my $maxSV = $self->_maxSV;

        if (defined $maxSV && ($type eq 'CTX' && $maxSV->{type} eq $type ||
	        $type eq 'INV' && $maxSV->{type} eq $type ||
		    (($type eq $maxSV->{type} && $type eq 'DEL') ||
		    ($type eq 'ITX' && ($maxSV->{type} eq 'ITX' || $maxSV->{type} eq 'INS')) ||
		    ($type eq 'INS' && ($maxSV->{type} eq 'ITX' || $maxSV->{type} eq 'INS'))) &&
                # $maxSV->{size}
                   $size >= $self->min_size_of_confirm_asm_sv && (!defined $self->invalid_indel_range || abs($maxSV->{size}-$size)<=$self->invalid_indel_range))) {
            my $scarstr = $maxSV->{scarsize}>0 ? substr($maxSV->{contig},$maxSV->{bkstart}-1,$maxSV->{bkend}-$maxSV->{bkstart}+1) : '-';

            #        printf STDOUT ("%s\t%d(%d)\t%s\t%d(%d)\t%s\t%d(%d)\t%s(%s)\t%s\t%d\t%d\t%d\%\t%d\t%d\t%d\t%d\t%d\t%s\t%s\ta%d.b%d\t%s\t%s\t%s\n",$maxSV->{chr1},$maxSV->{start1},$start,$maxSV->{chr2},$maxSV->{start2},$end,$maxSV->{ori},$maxSV->{size},$size,$maxSV->{type},$type,$maxSV->{het},$maxSV->{weightedsize},$maxSV->{read_len},$maxSV->{fraction_aligned}*100,$maxSV->{n_seg},$maxSV->{n_sub},$maxSV->{n_indel},$maxSV->{nbp_indel},$maxSV->{microhomology},$scarstr,$prefix,50,100, 'NA', 'NA', 'NA');
            $out_fh->printf("%s\t%d(%d)\t%s\t%d(%d)\t%s\t%d(%d)\t%s(%s)\t%s\t%d\t%d\t%d\%\t%d\t%d\t%d\t%d\t%d\t%s\t%s\ta%d.b%d\t%s\t%s\t%s\n",$maxSV->{chr1},$maxSV->{start1},$start,$maxSV->{chr2},$maxSV->{start2},$end,$maxSV->{ori},$maxSV->{size},$size,$maxSV->{type},$type,$maxSV->{het},$maxSV->{weightedsize},$maxSV->{read_len},$maxSV->{fraction_aligned}*100,$maxSV->{n_seg},$maxSV->{n_sub},$maxSV->{n_indel},$maxSV->{nbp_indel},$maxSV->{microhomology},$scarstr,$prefix,50,100, 'NA', 'NA', 'NA');
        
            if ($bp_io) {  #save breakpoint sequence
                my $coord = join(".",$maxSV->{chr1},$maxSV->{start1},$maxSV->{chr2},$maxSV->{start2},$maxSV->{type},$maxSV->{size},$maxSV->{ori});
                my $coord_pipe = join("|",$maxSV->{chr1},$maxSV->{start1},$maxSV->{chr2},$maxSV->{start2},$maxSV->{type},$maxSV->{size},$maxSV->{ori});
                my $contigsize = $maxSV->{contiglens};
                my $seqobj = Bio::Seq->new( 
                    -display_id => "ID:$prefix,Var:$coord,Ins:$maxSV->{bkstart}\-$maxSV->{bkend},Length:$contigsize,KmerCoverage:$maxSV->{contigcovs},Strand:$maxSV->{strand},Assembly_Score:$maxSV->{weightedsize},PercNonRefKmerUtil:$maxSV->{kmerutil},Ref_start:$maxSV->{refpos1},Ref_end:$maxSV->{refpos2},Contig_start:$maxSV->{rpos1},Contig_end:$maxSV->{rpos2},CrossMatch:$coord_pipe,TIGRA",
                    -seq => $maxSV->{contig}, 
                );
                $bp_io->write_seq($seqobj);
            }

            if ($cm_aln_fh) {
                $cm_aln_fh->printf("%s\t%d(%d)\t%s\t%d(%d)\t%s\t%d(%d)\t%s(%s)\t%s\t%d\t%d\t%d\%\t%d\t%d\t%d\t%d\t%d\t%s\t%s\ta%d.b%d\n",$maxSV->{chr1},$maxSV->{start1},$start,$maxSV->{chr2},$maxSV->{start2},$end,$maxSV->{ori},$maxSV->{size},$size,$maxSV->{type},$type,$maxSV->{het},$maxSV->{weightedsize},$maxSV->{read_len},$maxSV->{fraction_aligned}*100,$maxSV->{n_seg},$maxSV->{n_sub},$maxSV->{n_indel},$maxSV->{nbp_indel},$maxSV->{microhomology},$scarstr,$prefix,'50','100');
                for my $aln (split /\,/, $maxSV->{alnstrs}) {
	                $cm_aln_fh->printf("%s\n", join("\t", split /\|/, $aln));
                }
                $cm_aln_fh->print("\n");
            }
        }
        elsif ($self->_unc_pred_fh) {
            $self->_unc_pred_fh->printf("%s\t%d\t%s\t%d\t%s\t%d\t%s\n",$chr1,$start,$chr2,$end,$type,$size,$ori);
        }
        $self->_maxSV(undef)         if $self->_maxSV; #reset for each SV
        $self->_N50size(undef)       if $self->_N50size;
        $self->_WeightAvgSize(undef) if $self->_WeightAvgSize;
    }

    if (!defined $self->intermediate_read_dir and defined $self->intermediate_save_dir) {
        my $cmd = "mv -f $datadir " . $self->intermediate_save_dir;
        system $cmd;
    }
    #elsif (!defined $self->intermediate_read_dir){
    #    File::Temp::cleanup();
    #}

    $out_fh->close;
    $cm_aln_fh->close if $cm_aln_fh;
    $self->_unc_pred_fh->close if $self->_unc_pred_fh;
    $self->debug_message('AssemblyValidation finished ok.');

    return 1;
}


sub _check_bam {
    my $self = shift;
    my $bam_files = $self->bam_files;
    my (@bam_files, $valid_bam);

    if ($bam_files =~ /\,/) {
        @bam_files = split /\,/, $bam_files; #TODO validation check each file
    }
    else {
        @bam_files = ($bam_files);
    }

    for my $bam (@bam_files) {
        unless (-s $bam) {
            $self->error_message("Bam file: $bam is not valid");
            die $self->error_message;
        }
        unless (-s $bam.'.bai') {
            $self->error_message("bam index file: $bam.bai is not valid");
            die $self->error_message;
        }
        $valid_bam .= ' ' . $bam;
    }
    return $valid_bam;
}


sub _get_tigra_options {
    my $self = shift;

    my %tigra_sv_options = (
        est_max_ins_size      => 'A',
        flank_size            => 'l',
        pad_local_ref         => 'w',
        map_qual_to_asm       => 'q',
        mismatch_limit        => 'N',
#        avg_read_depth_limit  => 'p',
        min_breakdancer_score => 'Q',
        skip_libraries        => 'L',
        skip_call            => 'z',
    );

    my $tigra_opts = '-d -r -I '. $self->_data_dir  . ' ';
    $tigra_opts .= '-b ' unless $self->custom_sv_format;
    $tigra_opts .= '-p 10000 ' if($self->asm_high_coverage);
    $tigra_opts .= '-h 300 ' if($self->asm_high_coverage);
    $tigra_opts .= '-z ' . $self->skip_call if($self->skip_call);
    $tigra_opts .= '-c ' . $self->specify_chr . ' ' if($self->specify_chr);
    $tigra_opts .= '-M ' . $self->min_size_of_confirm_asm_sv . ' ' if($self->min_size_of_confirm_asm_sv);
    $tigra_opts .= '-R ' . $self->reference_file . ' ' if($self->reference_file);
    $tigra_opts .= '-k ' . '15,25,39' . ' ' if($self->asm_high_coverage);

    for my $opt (keys %tigra_sv_options) {
        if ($self->$opt) {
            $tigra_opts .= '-'.$tigra_sv_options{$opt}.' '.$self->$opt . ' ';
        }
    }
    return $tigra_opts;
}


sub _cross_match_validation {
    my ($self, $ctg_type, $tigra_sv_name, $tigra_sv_filename) = @_;
    my ($chr1,$start,$chr2,$end,$type,$size,$ori,$regionsize,$pos1,$pos2) = split /\./, $tigra_sv_name;

    my $output_map = $self->_tigra_sv_output_map;
    if (exists $output_map->{$tigra_sv_filename}) {
        my $p = $output_map->{$tigra_sv_filename};
        $chr1 = $p->{chr1};
        $start = $p->{start};
        $chr2 = $p->{chr2};
        $end = $p->{end};
        $type = $p->{type};
        $size = $p->{size};
        $ori = $p->{ori};
        $regionsize = $p->{region_size};
        $pos1 = $p->{pos1};
        $pos2 = $p->{pos2};
    }

    my $posstr = join '_', $chr1, $pos1, $chr2, $pos2, $type, $size, $ori;
    my $head   = join '.', $chr1, $start, $chr2, $end, $type, $size, $ori;

    my $datadir = $self->_data_dir;
    my $ref_fa  = $datadir . "/$head.ref.fa"; 
#    my $ref_fa = $datadir . "/$head.1.ref.fa";
    my ($tigra_sv_fa, $cm_out);

    if ($ctg_type eq 'homo') {
        $tigra_sv_fa = $datadir . "/$tigra_sv_name.fa.contigs.fa";
        $cm_out      = $datadir . "/$tigra_sv_name.stat";
    }
    elsif ($ctg_type eq 'het') {
        $tigra_sv_fa = $datadir . "/$tigra_sv_name.fa.contigs.het.fa";
        $cm_out      = $datadir . "/$tigra_sv_name.het.stat";
    }
    else {
        $self->error_message("Wrong type: $ctg_type");
        die $self->error_message;
    }
    
    unless (-s $tigra_sv_fa) {
        $self->warning_message("tigra sv fasta: $tigra_sv_fa is not valid. Skip this $ctg_type cross_match run");
        return;
    }

    unless (-s $ref_fa){
        $self->warning_message("tigra sv ref: $ref_fa is not valid. Skip this $ref_fa cross_match run");
        return;
    }

    my $cm_cmd_opt = '-bandwidth 20 -minmatch 20 -minscore 25 -penalty '.$self->cm_sub_penalty.' -discrep_lists -tags -gap_init '.$self->cm_gap_init_penalty.' -gap_ext -1';
	my $cm_cmd = "cross_match $tigra_sv_fa $ref_fa $cm_cmd_opt > $cm_out 2>/dev/null";
	           
    my $rv = Genome::Sys->shellcmd (
        cmd           => $cm_cmd,
        input_files   => [$tigra_sv_fa, $ref_fa],
        #output_files => [$cm_out],
        #allow_zero_size_output_files => 1,
    );
    
    unless ($rv) {
        $self->error_message("Running cross_match for $tigra_sv_name homo failed.\nCommand: $cm_cmd");
        die $self->error_message;
    }
    $self->debug_message("Cross_match for $type contigs: $tigra_sv_name Done");
        
    my $makeup_size      = 0; # by default they are zero
    my $concatenated_pos = 0; # by default they are zero

    if ($size > 99999) {
        $makeup_size      = ($end - $self->flank_size) - ($start + $self->flank_size) - 1;
        $concatenated_pos = 2 * $self->flank_size + $self->pad_local_ref;
    }

    $self->debug_message("GetCrossMatchIndel for $tigra_sv_name $type");
    my $cm_indel = Genome::Model::Tools::Sv::CrossMatchForIndel->create(
        cross_match_file     => $cm_out,
        local_ref_seq_file   => $ref_fa,
        assembly_contig_file => $tigra_sv_fa,
        per_sub_rate         => $self->max_polymorphism_rate,
        ref_start_pos        => $posstr,
        ref_cat_pos          => $concatenated_pos,
        variant_size         => $makeup_size,
    );
    my $result = $cm_indel->execute;

    # unless explicitly deleted, this rather large object will stay in the UR
    # cache forever
    $cm_indel->delete;
    $cm_indel = undef;

    if ($result && $result =~ /\S+/) {
	    $self->_UpdateSVs($result,$makeup_size,$regionsize,$tigra_sv_fa,$ctg_type, $cm_out);
    }
    
    return 1;
}


sub _UpdateSVs{
    my ($self,$result,$makeup_size,$regionsize,$tigra_sv_fa,$type, $cm_out) = @_;
    my $datadir = $self->_data_dir;
    my $maxSV   = $self->_maxSV;
    my $N50size = $self->_N50size;
    my $depthWeightedAvgSize = $self->_WeightAvgSize;

    if (defined $result) {
        print STDERR "Result: $result\n";
        my ($pre_chr1,$pre_start1,$pre_chr2,$pre_start2,$ori,$pre_bkstart,$pre_bkend,$pre_size,$pre_type,$pre_contigid,$alnscore,$scar_size,$read_len,$fraction_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$strand,$microhomology,$alnstrs) = split /\s+/, $result;
        #$pre_size += $makeup_size if $n_seg >= 2;
        if (defined $pre_size && defined $pre_start1 && defined $pre_start2) {
            my ($contigseq,$contiglens,$contigcovs,$kmerutil) = _GetContig($tigra_sv_fa, $pre_contigid);
            my ($refpos1, $refpos2, $rpos1, $rpos2) = _GetRefPos($cm_out, $pre_contigid,$pre_size,$pre_type);
            $alnscore = int($alnscore*100/$regionsize); 
            $alnscore = $alnscore>100 ? 100 : $alnscore;
            if (!defined $maxSV || $maxSV->{size}<$pre_size || $maxSV->{alnscore} < $alnscore) {
	            my $N50score = int($N50size*100/$regionsize); 
                $N50score = $N50score>100 ? 100 : $N50score;
	            if ($self->assemble_mouse) {  #Mouse
	                $pre_chr1 =~ s/.*\///; 
                    $pre_chr1 =~ s/\.fasta//;
	                $pre_chr2 =~ s/.*\///; 
                    $pre_chr2 =~ s/\.fasta//;
                }
	            ($maxSV->{chr1},$maxSV->{start1},$maxSV->{chr2},$maxSV->{start2},$maxSV->{bkstart},$maxSV->{bkend},$maxSV->{size},$maxSV->{type},$maxSV->{contigid},$maxSV->{contig},$maxSV->{contiglens},$maxSV->{contigcovs},$maxSV->{kmerutil},$maxSV->{N50},$maxSV->{weightedsize},$maxSV->{alnscore},$maxSV->{scarsize},$maxSV->{a},$maxSV->{b},$maxSV->{read_len},$maxSV->{fraction_aligned},$maxSV->{n_seg},$maxSV->{n_sub},$maxSV->{n_indel},$maxSV->{nbp_indel},$maxSV->{strand},$maxSV->{microhomology},$maxSV->{refpos1},$maxSV->{refpos2},$maxSV->{rpos1},$maxSV->{rpos2}) = ($pre_chr1,$pre_start1,$pre_chr2,$pre_start2,$pre_bkstart,$pre_bkend,$pre_size,$pre_type,$pre_contigid,$contigseq,$contiglens,$contigcovs,$kmerutil,$N50score,$depthWeightedAvgSize,$alnscore,$scar_size,'50','100',$read_len,$fraction_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$strand,$microhomology,$refpos1,$refpos2,$rpos1,$rpos2);
                print STDERR "\n\nWeightedsize: $maxSV->{weightedsize}\n";
	            $maxSV->{het}     = $type;
	            $maxSV->{ori}     = $ori;
	            $maxSV->{alnstrs} = $alnstrs;

                    # add according to Ken's requirement, now cut and add it in CrossMatchIndel.pm for only INDEL
#                    if($maxSV->{strand} =~ /-/ && $maxSV->{type} =~ /DEL/i){
#                        $maxSV->{start2} -= 1;
#                        $maxSV->{bkend} -= 1;
#                    }
            }
        }
    }
    $self->_maxSV($maxSV);
    return 1;
}

# only good for small INDELs
sub _GetRefPos{
    my ($fin, $contigid, $pre_size, $pre_type) = @_;

    my ($ref_start, $ref_end, $r_start, $r_end);
    my ($ref_start1, $ref_end1, $r_start1, $r_end1);
    my ($refpos1, $refpos2, $rpos1, $rpos2);
    my ($chr, $ref_1, $ref_2);
    my $num = 0;
    my ($r1, $r2, $r3, $type, $str, $g1, $g2, $g3);
    my $check = 0;
    my $check_results = 0;
    open(CM,"<$fin") || die "unable to open $fin\n";
    while(<CM>){
        # go ahead to check DISCREPANCY contains I or D directly after find a matched ALIGNMENT
        if($check == 1 && $_ =~ /^DISCREPANCY/){
            if($pre_type =~ /DEL/){
                if($_ =~ /D-$pre_size/){
                    $check_results = 1;
                    last;
                }
            }
            elsif($pre_type =~ /INS/){
                if($_ =~ /I-$pre_size/){
                    $check_results = 1;
                    last;
                }
            }
        }
        else{
            $check = 0;
        }

        if($_ =~ /^ALIGNMENT/ && $_ =~ /$contigid/){
            my @a = split(/\s+/, $_);
            next if($a[5] ne $contigid);

            if($#a == 13){
                ($r1, $r2, $r3, $type, $str, $g1, $g2, $g3) = ($a[$#a - 7], $a[$#a - 6], $a[$#a - 5], $a[$#a - 4], $a[$#a - 3], $a[$#a - 2], $a[$#a - 1], $a[$#a]);
            }
            else{
                ($r1, $r2, $r3, $str, $g1, $g2, $g3) = ($a[$#a - 6], $a[$#a - 5], $a[$#a - 4], $a[$#a - 3], $a[$#a - 2], $a[$#a - 1], $a[$#a]);
            }
            
            #print "HHHHHHHHHHHHHH: $r1, $r2, $r3\n";
            ($chr, $ref_1, $ref_2) = ($str =~ /(\S+):(\d+)-(\d+)/);
=cut                
            #print "HHHHHHHHHHHHHHH: $g1\t$g2\t$g3\n";
            if($g1 =~ /\(/){
                ($ref_start, $ref_end) = ($g2, $g3) if($num == 0);
                ($ref_start1, $ref_end1) = ($g2, $g3) if($num == 1);
            }
            if($r1 =~ /\(/){
                ($r_start, $r_end) = ($r2, $r3) if($num == 0);
                ($r_start1, $r_end1) = ($r2, $r3) if($num == 1);
            }
            if($g3 =~ /\(/){
                ($ref_start, $ref_end) = ($g1, $g2) if($num == 0);
                ($ref_start1, $ref_end1) = ($g1, $g2) if($num == 1);
            }
            if($r3 =~ /\(/){
                ($r_start, $r_end) = ($r1, $r2) if($num == 0);
                ($r_start1, $r_end1) = ($r1, $r2) if($num == 1);
            }
=cut
            $check = 1;
#$num ++;
        }
    }

    if($check_results == 1){
        ($refpos1, $refpos2) = ($g2, $g3) if($g1 =~ /\(/);
        ($refpos1, $refpos2) = ($g1, $g2) if($g3 =~ /\(/);
        ($rpos1, $rpos2) = ($r1, $r2) if($r3 =~ /\(/);
        ($rpos1, $rpos2) = ($r2, $r3) if($r1 =~ /\(/);
    }
    else{
        print STDERR "Error, didn't find $contigid of $pre_type and $pre_size in $fin";
    }
=cut
    if($num >= 2){
        #$refpos1 = $ref_start > $ref_start1 ? $ref_start1 : $ref_start;
        #$refpos2 = $ref_end > $ref_end1 ? $ref_end : $ref_end1;
        #$rpos1 = $refpos1 == $ref_start ? $r_start : $r_start1;
        #$rpos2 = $refpos2 == $ref_end ? $r_end : $r_end1;
        $rpos1 = $r_start > $r_start1 ? $r_start1 : $r_start;
        $rpos2 = $r_end > $r_end1 ? $r_end : $r_end1;
        $refpos1 = $rpos1 == $r_start ? $ref_start : $ref_start1;
        $refpos2 = $rpos2 == $r_end ? $ref_end : $ref_end1;
    }
    else{
        $refpos1 = $ref_start;
        $refpos2 = $ref_end;
        $rpos1 = $r_start;
        $rpos2 = $r_end;
    }
=cut
    $refpos1 += $ref_1;
    $refpos2 += $ref_1;


#    print "$refpos1\t$refpos2\n";
    return ($refpos1, $refpos2, $rpos1, $rpos2);
}






sub _GetContig{
    my ($fin, $contigid) = @_;
    my $in = Bio::SeqIO->newFh(-file => $fin, '-format' => 'Fasta');
    my $sequence;
    my @info;
    while (my $seq = <$in>) {
    # do something with $seq
        next unless $seq->id eq $contigid;
        @info = split /\s+/, $seq->desc;
        $sequence = $seq->seq();
        last;
    }
    return ($sequence, $info[0], $info[1], $info[$#info]);
}


sub _ComputeTigraN50{
    my ($contigfile) = @_;
    my @sizes;
    my $totalsize = 0;

    my $fh = Genome::Sys->open_file_for_reading($contigfile) or return;
    while (my $l = $fh->getline){
        chomp $l;
        next unless $l =~ /^\>/;
        next if $l =~ /\*/;
        my ($id, $size, $depth, $ioc, @extra) = split /\s+/, $l;
        next if $size <= 50 && $depth < 3 && $ioc =~ /0/;
        push @sizes, $size;
        $totalsize += $size;
    }
    $fh->close;
    my $cum_size = 0;
    my $halfnuc_size = $totalsize/2;
    my $N50size;
    @sizes = sort {$a<=>$b} @sizes;
    while ($cum_size < $halfnuc_size) {
        $N50size = pop @sizes;
        $cum_size += $N50size;
    }
    return $N50size;
}


sub _ComputeTigraWeightedAvgSize{
    my $contigfile = shift;
    my $totalsize  = 0;
    my $totaldepth = 0;

    my $fh = Genome::Sys->open_file_for_reading($contigfile) or return;
    while (my $l = $fh->getline) {
        chomp $l;
        next unless $l =~ /^\>/;
        next if $l =~ /\*/;
        my ($id,$size,$depth,$ioc,@extra) = split /\s+/, $l;
        next if $size <= 50 && (($depth<3 && $ioc=~/0/) || $depth>500);  #skip error tips or extremely short and repetitive contigs
        $l = $fh->getline;
        chomp $l;
        #$_=<CF>; chomp;
        #next if($size<=50 && (/A{10}/ || /T{10}/ || /C{10}/ || /G{10}/));  #ignore homopolymer contig
        next if $size <= 50 && $l =~ /A{10}|T{10}|C{10}|G{10}/;
        $totalsize += $size*$depth;
        $totaldepth+= $depth;
    }
    $fh->close;
    return $totaldepth > 0 ? int($totalsize*100/$totaldepth)/100 : 0;
}

1;
