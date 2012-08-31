package Genome::Model::Tools::Annotate::1000GenomesFrequency;

use strict;
use warnings;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Annotate::1000GenomesFrequency {
    is => 'Command',                       
    has => [                                # specify the command's single-value properties (parameters) <--- 
    variant_file => { },
    tabix_database => {default => '/gscmnt/sata921/info/medseq/1000genomes_vcf/ALL.2of4intersection.20100804.genotypes.vcf.gz' },
    novel_output => { },
    known_output => { },
    allele_match => { default => '1', is_optional=>1 },
    genome_build => { default => '36' },
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
    while(my $line = $bed_fh->getline) {
        chomp($line);
        my ($chr, $start, $stop, $ref_var) = split /\t/, $line;
        my ($ref, $var) = split /\//, $ref_var;
        my @vcf_line =  `tabix $tabix_db $chr:$stop-$stop`;
        die if (scalar(@vcf_line) >1);
        
        my ($vcf_chr, $vcf_pos, $marker_id, $vcf_ref, $vcf_var, $qual, $filter, $info);
        if(@vcf_line) {
            ($vcf_chr, $vcf_pos, $marker_id, $vcf_ref, $vcf_var, $qual, $filter, $info, undef) = split /\t/, $vcf_line[0];
        }
        $marker_id = '-' unless $marker_id;
        $info = '-' unless $info;
        my $new_bed_line = $line . "\t$marker_id;$info";
        if(scalar(@vcf_line) < 1) { ###not found
            $novel_tmp->print("$new_bed_line\n");
        }
        elsif($self->allele_match) {  ##check for variant match to determine novelness
            if($vcf_var eq $var ) {
                $known_tmp->print("$new_bed_line\n");
            }
            else {
                $novel_tmp->print("$new_bed_line\n");
            }
        }
        else {  ##we got a positional hit and allele_match is not on print to known
            $known_tmp->print("$new_bed_line\n");
        }
    }

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
    my $novel_lines = $lines - $unmapped_lines - $known_lines;
    my $total = $lines - $unmapped_lines;
    my $novel_percentage = sprintf( "%.2f", ($novel_lines)/($total) * 100);
    my $known_percentage = sprintf( "%.2f", ($known_lines)/($total) * 100);
    print STDERR "$novel_lines novel variants for $novel_percentage% of mapped(liftOvered) total($total)\n";
    print STDERR "$known_lines previously discovered variants for $known_percentage% of mappeD(liftOvered) total($total)\n";
    return 1;
}







1;
