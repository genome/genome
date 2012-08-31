package Genome::Model::Tools::Annotate::1000GenomesFrequency2;

use strict;
use warnings;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Annotate::1000GenomesFrequency2 {
    is => 'Command',                       
    has => [                                # specify the command's single-value properties (parameters) <--- 
    variant_file => { },
    tabix_database => {default => '/gscmnt/gc6132/info/medseq/1000_genomes/downloads/2012-03-27/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz', },
    novel_output => { },
    known_output => { },
    allele_match => { default => '1', is_optional=>1 },
    genome_build => { default => '37' },
    liftover_chain_path => { default => '/gscmnt/sata847/info/liftOver_chain_files/' },
    ], 
};


sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "filter out sites with 1000 genomes frequency",
}

sub help_synopsis {                         # replace the text below with real examples <---
    return <<EOS
gmt annotation 1000-genomes-frequency --variant-file=<anno or bed format>  --novel_output=<bed_output> --known_output=<bed_output> --allele-match=<0 or 1> --genome-build=<36 or 37>
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
    detailed help WOULD go here. 
EOS
}

sub execute {                               # replace with real execution logic.
    my $self = shift;
    $DB::single=1;
    my ($bed_fh, $bed_path);
    my $input_fh = Genome::Sys->open_file_for_reading($self->variant_file);
    chomp(my $line = $input_fh->getline);
    my @fields = split /\t/, $line;

    ####prepare input for liftover in bed_fh temp file or use original file
    if($fields[3] =~ m/\//) {  ####bed file
        $input_fh->close;
        $bed_fh = Genome::Sys->open_file_for_reading($self->variant_file);
        $bed_path = $self->variant_file;
    }
    else {  ###anno format
        ($bed_fh, $bed_path) = Genome::Sys->create_temp_file();
        do {
            chomp($line);
            my ($chr, $start, $stop, $ref, $var, undef) = split /\t/, $line;
            $start = $stop -1;
            $bed_fh->print("$chr\t$start\t$stop\t$ref/$var\n");
        }while ($line = $input_fh->getline);
        $bed_fh->close;
        $bed_fh = Genome::Sys->open_file_for_reading($bed_path);
    }
    ######## done with prepping input

    my ($lines, $unmapped_lines);
    ###liftover to build37 if necessary###
    my $build37_bed = Genome::Sys->create_temp_file_path('build37_mapped');
    my $unmapped_bed = Genome::Sys->create_temp_file_path('build37_unmapped');
    if($self->genome_build == '36') {
        my $liftover_file = $self->liftover_chain_path . "hg18ToHg19.over.chain.gsc_chromosomes";
        Genome::Sys->shellcmd(
            cmd=> "liftOver $bed_path $liftover_file $build37_bed $unmapped_bed",
            skip_if_output_is_present=>0,
        );
        chomp($unmapped_lines = `cat $unmapped_bed | grep -v "^#" | wc -l`);
        chomp($lines = `cat $bed_path | wc -l`);
        my $percent_unmapped = sprintf( "%.2f", $unmapped_lines/$lines * 100);
        print STDERR "$unmapped_lines variants were not mapped during liftover, $percent_unmapped% of total($lines)\n";
        ($bed_fh, $bed_path) = Genome::Sys->open_file_for_reading($build37_bed);
    }
    else {
        chomp($lines = `cat $bed_path | wc -l`);
        $unmapped_lines=0;
    }

    ###########end liftover###############
    my $tabix_db = $self->tabix_database;
    my ($novel_tmp, $novel_tmp_path ) = Genome::Sys->create_temp_file('novel_variants_temp');
    my ($known_tmp, $known_tmp_path ) = Genome::Sys->create_temp_file('known_variants_temp');
    $novel_tmp->print("CHR\tSTART\tSTOP\tREF/VAR\tAF\tAFR_AF\tEUR_AF\n");
    $known_tmp->print("CHR\tSTART\tSTOP\tREF/VAR\tAF\tAFR_AF\tEUR_AF\n");
    while(my $line = $bed_fh->getline) {
        chomp($line);
        my ($chr, $start, $stop, $ref_var) = split /\t/, $line;
        my ($ref, $var) = split /\//, $ref_var;
        my @vcf_line =  `tabix $tabix_db $chr:$stop-$stop`;
        die "Multiple lines found for $line\n" if (scalar(@vcf_line) >1);
        my ($vcf_chr, $vcf_pos, $marker_id, $vcf_ref, $vcf_var, $qual, $filter, $info);

        if(scalar(@vcf_line) == 1) {
            chomp($vcf_line[0]);
            ($vcf_chr, $vcf_pos, $marker_id, $vcf_ref, $vcf_var, $qual, $filter, $info, undef) = split /\t/, $vcf_line[0];
            my @info_fields = split(/\;/, $info);
            my($AF, $EUR_AF, $AFR_AF);
            for my $info_field (@info_fields) {
                if($info_field =~m/EUR_AF/) {
                    $info_field =~ s/EUR_AF=//;
                    $EUR_AF=$info_field;
                }elsif($info_field =~m/AFR_AF/){
                    $AFR_AF=$info_field;
                    $AFR_AF=~ s/AFR_AF=//;
                } elsif($info_field =~m/^AF/) {
                    $AF = $info_field;
                    $AF =~ s/AF=//;
                }
            }
            unless(defined($AFR_AF)) {
                $AFR_AF="N/A";
            }
            unless(defined($EUR_AF)) {
                $EUR_AF = "N/A";
            }
            unless(defined($AF)) {
                die "Can't find AF Info code in line: $vcf_line[0]\n";
            }
#            if ($AFR_AF =~ EUR_AF) {
            $known_tmp->print("$vcf_chr\t$start\t$vcf_pos\t$vcf_ref/$vcf_var\t$AF\t$AFR_AF\t$EUR_AF\n");
            #           } else {
            #  $known_fh->print("$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]\t$chr[10]\t0\t$AFR_AF\t$EUR_AF\n");
#            }

        } else { #we didn't get a line back, there's no 1000 genomes record for this event
            $novel_tmp->print("$chr\t$start\t$stop\t\t\t\n");
        } 
    }
    $novel_tmp->close();
    $known_tmp->close();

    my $novel_file = $self->novel_output;
    my $known_file = $self->known_output;

    ###re liftover###
    if($self->genome_build == '36') {
        ##lift that shit back!
        my $liftover_file = $self->liftover_chain_path . "hg19ToHg18.over.chain.gsc_chromosomes";
        Genome::Sys->shellcmd(
            cmd=> "liftOver $novel_tmp_path $liftover_file $novel_file /dev/null",
            skip_if_output_is_present=>0,
        );
        Genome::Sys->shellcmd(
            cmd=> "liftOver $known_tmp_path $liftover_file $known_file /dev/null",
            skip_if_output_is_present=>0,
        );
    } else {
        ###no liftover needed just copy back to destination
        Genome::Sys->shellcmd(cmd=>"cp $known_tmp_path $known_file");
        Genome::Sys->shellcmd(cmd=>"cp $novel_tmp_path $novel_file");
    }     


    chomp(my $known_lines = `cat $known_file | wc -l`);
    $known_lines-=1; #account for header we added
    my $novel_lines = $lines - $unmapped_lines - $known_lines;
    my $total = $lines - $unmapped_lines;
    my $novel_percentage = sprintf( "%.2f", ($novel_lines)/($total) * 100);
    my $known_percentage = sprintf( "%.2f", ($known_lines)/($total) * 100);
    print STDERR "$novel_lines novel variants for $novel_percentage% of mapped(liftOvered) total($total)\n";
    print STDERR "$known_lines previously discovered variants for $known_percentage% of mappeD(liftOvered) total($total)\n";
    return 1;
}






=cut
my $AF = $info_fields[10];
            $info_fields[10] = 0 unless $AF =~ /^\S+/;
            my $col11 = $info_fields[11];
            #print "col 11 : $col11\n";           
            my $col12 = $info_fields[12];
            my $AFR_AF = $info_fields[13];
            $DB::single=1 unless ($info_fields[13]);
            $info_fields[13] = 0 unless $AFR_AF =~ /^\S+/;

            my $EUR_AF = $info_fields[14];
            $info_fields[14] = 0 unless $EUR_AF =~ /^\S+/;

            if ($info_fields[11] =~ "AFR_AF") {
                $AFR_AF = $info_fields[11];
            }
            if ($info_fields[11] =~ "EUR_AF") {
                $EUR_AF = $info_fields[11];
            }
            if ($info_fields[12] =~ "AFR_AF") {
                $AFR_AF = $info_fields[12];
            }
            if ($info_fields[12] =~ "EUR_AF") {
                $EUR_AF = $info_fields[12];
            }
            $AFR_AF = 0 unless $AFR_AF =~ /^\S+/;
            $EUR_AF = 0 unless $EUR_AF =~ /^\S+/;
            $AF = 0 unless $AF =~ /^\S+/;
=pod



1;
