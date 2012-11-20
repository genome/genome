package Genome::Model::Tools::Annovar::AnnotateVariation;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annovar::AnnotateVariation {
    is => 'Genome::Model::Tools::Annovar::Base',
    has_input => [
        buildver => {
            is => "String",
        },
        table_names => {
            is => "String",
            is_many => 1,
        },
        input_file => {
            is => 'File',
        },
        outfile => {
            is => 'Path',
        },
        annotation_type => {
            is => 'String',
            valid_values => ["geneanno", "regionanno"],
        },
    ],
    has_optional_input => [
        scorecolumn => {
            is => 'Number',
        },
    ],
};

sub execute {
    my $self = shift;
    my $input_file = $self->input_file;
    if (_is_vcf($self->input_file)) {
        $self->status_message("Converting VCF to annovar");
        my $converter_input = $self->input_file;
        if (Genome::Sys->file_is_gzipped($converter_input)) {
            $self->status_message("Uncompressing VCF file");
            $converter_input = Genome::Sys->create_temp_file_path;
            my $in = Genome::Sys->open_gzip_file_for_reading($self->input_file);
            my $out = Genome::Sys->open_file_for_writing($converter_input);
            while(my $line = <$in>) {
                $out->print($line);
            }
            $in->close;
            $out->close;
        }
        $input_file = Genome::Sys->create_temp_file_path;
        my $convert = $self->script_path."convert2annovar.pl $converter_input -format vcf4 -outfile ".$input_file;
        $self->status_message("Running convert command: $convert");
        my $rv = Genome::Sys->shellcmd(cmd => $convert);
        unless ($rv) {
            $self->error_message("Could not convert vcf to annovar input");
            return $rv;
        }
    }


    foreach my $table_name ($self->table_names) {
        my $db = Genome::Model::Tools::Annovar::Db->get_or_create(
            buildver => $self->buildver,
            table_name => $table_name,
            test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
        );

        unless ($db) {
            $self->error_message("Couldn't get ANNOVAR db result for ".$table_name." buildver ".$self->buildver);
            return;
        }

        $self->status_message("Annotating with ".$table_name." in ".$db->output_dir);

        my $cmd = $self->script_path."annotate_variation.pl --outfile ".$self->outfile." --buildver ".$self->buildver;
        if ($self->annotation_type eq "geneanno") {
            $cmd .= " --geneanno ";
        }
        elsif ($self->annotation_type eq "regionanno") {
            if (defined $self->scorecolumn) {
                $cmd .= " -scorecolumn ".$self->scorecolumn." ";
            }
            $cmd .= " -regionanno -dbtype ".$table_name;
        }
        $cmd .= " ".$input_file." ".$db->output_dir;
        $self->status_message("Executing command $cmd");
        my $rv = Genome::Sys->shellcmd(cmd => $cmd);
        unless ($rv) {
            $self->error_message("Failed to annotate with table $table_name");
            return;
        }
    }

    my %variants;

    foreach my $table_name ($self->table_names) {
        my $annotated_file = $self->outfile.".".$self->buildver."_".$table_name;
        my $in = Genome::Sys->open_file_for_reading($annotated_file);
        while (my $line = <$in>) {
            chomp $line;
            my @fields = split(/\t/, $line);
            if (defined $variants{$fields[2]}{$fields[3]}{$fields[4]}){
                $variants{$fields[2]}{$fields[3]}{$fields[4]}{$fields[0]} = $fields[1];
            }
            else {
                $variants{$fields[2]}{$fields[3]}{$fields[4]} = {$fields[0] => $fields[1]};
            }
        }
        $in->close;
    }
    

    my $in = Genome::Sys->open_file_for_reading($input_file);
    my $out = Genome::Sys->open_file_for_writing($self->outfile.".summary");
    #print header
    $out->print(join("\t", "#Chr", "Start", "Stop", "RefAllele", "VarAllele", $self->table_names)."\n");
    
    while (my $line = <$in>) {
        chomp $line;
        my @fields = split(/\t/, $line);
        $out->print(join("\t", $fields[0], $fields[1], $fields[2], $fields[3], $fields[4]));
        foreach my $table_name ($self->table_names) {
            $out->print("\t");
            if (defined ($variants{$fields[0]}{$fields[1]}{$fields[2]}{$table_name})){
                $out->print($variants{$fields[0]}{$fields[1]}{$fields[2]}{$table_name});
            }
            else {
                $out->print(".");
            }
        }
        $out->print("\n");
    }
    $out->close;
    return 1;
}

sub _is_vcf {
    my $in_file = shift;
    return ($in_file =~ /vcf$/ or $in_file =~ /vcf.gz$/);
}

1;

