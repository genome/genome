package Genome::Model::Tools::Assembly::Ace;

use strict;
use warnings;

use Genome;
use IO::File;
use Data::Dumper;

class Genome::Model::Tools::Assembly::Ace {
    is => 'Command',
    has => [ ],
};

sub help_brief {
    'Tools to export or remove contigs in ace file'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools assembly ace
EOS
}

sub xhelp_detail {
    return <<EOS
EOS
}
#validate contigs list and contig names in list
sub get_valid_contigs_from_list {
    my $self = shift;
    my $contigs_list = {}; #return a hash of contig names
    if ( $self->contigs ) {
        for my $contig ( $self->contigs ) {
            $contigs_list->{$contig} = 1;
        }
    }
    if ( $self->contigs_list) {
        my $fh = Genome::Sys->open_file_for_reading( $self->contigs_list );
        while (my $line = $fh->getline) {
            next if $line =~ /^\s+$/;
            chomp $line;
            my @ar = split(/\s+/, $line); #Contig9.99 Contig10.1
            foreach (@ar) {
                $contigs_list->{$_} = 1;
            }
        }
        $fh->close;
    }
    return $contigs_list;
}
#validate ace files provided
sub get_valid_input_acefiles {
    my $self = shift;
    #validate ace input params
    unless ($self->ace_list or $self->ace_files) {
	$self->error_message("ace or list of ace files must be supplied");
	return;
    }
    my @acefiles;
    if($self->ace_list) {
	unless (-s $self->ace_list) {
	    $self->error_message("Can not find ace list file or file is zero size: ".$self->ace_list);
	    return;
	}
	my $fh = IO::File->new("<".$self->ace_list) ||
	    die "Can not create file handle to read ace list: ".$self->ace_list."\n";
	while (my $line = $fh->getline) {
	    next if $line =~ /^\s+$/;
	    chomp $line;
	    my @ar = split(/\s+/, $line);
	    foreach (@ar) {
		my $acefile = (-s $_) ? $_ : $self->directory.'/'.$_;
		unless (-s $acefile) {
		    $self->error_message("Can not find ace file or file is zero size: ".$acefile);
		    return;
		}
		push @acefiles, $acefile;
	    }
	}
    }
    #check string of ace files
    if ($self->ace_files) {
	foreach ( $self->ace_files ) {
	    my $acefile = (-s $_) ? $_ : $self->directory.'/'.$_;
	    unless (-s $acefile) {
		$self->error_message("Can not find ace file or file is zero size: ".$acefile);
		return;
	    }
	    push @acefiles, $acefile;
	}
    }
    return \@acefiles;
}
#selectively print contig lines for exported contigs
sub filter_ace_files {
    my ($self, $acefiles, $contigs, $action) = @_;

    my @new_aces; #return ary ref of new ace names
    foreach my $acefile (@$acefiles) {
	my $ace_name = File::Basename::basename($acefile);
	my $int_file = $self->intermediate_file_name( $ace_name );
        unlink $int_file;
	my $int_ace_fh = Genome::Sys->open_file_for_writing( $int_file );
	my $export_setting = 0; #set to print ace contents if == 1
	my $contig_count = 0;
	my $total_read_count = 0;
	my $fh = Genome::Sys->open_file_for_reading( $acefile );
	while (my $line = $fh->getline) {
	    next if $line =~ /^AS\s+/;
	    if ($line =~ /^CO\s+/) {
		my $contig_name = $self->get_contig_name_from_ace_CO_line($line);
		$self->status_message("\tChecking $contig_name\n");
		if($action =~ /remove/) {
		    $export_setting = (exists $contigs->{$contig_name}) ? 0 : 1;
		}
		else {
		    $export_setting = (exists $contigs->{$contig_name}) ? 1 : 0;
		}
		#keep count of contigs and reads of those contigs that will be exported
		if ($export_setting == 1) {
		    my $read_count = $self->get_read_count_from_ace_CO_line($line);
		    $contig_count++;
		    $total_read_count += $read_count;
		}
	    }
	    last if ($line =~ /^CT{/ or $line =~ /^WA{/);
	    $int_ace_fh->print ($line) if $export_setting == 1;
	}
	$fh->close;
	$int_ace_fh->close;
	#re-write the ace file with correct contig and read counts
	my $ace_ext = ($action =~ /remove/) ? $ace_name.'.contigs_removed' : $ace_name.'.exported_contigs';
	my $final_ace;
	unless ($final_ace = $self->rewrite_ace_file($int_file, $contig_count, $total_read_count, $ace_ext)) {
	    $self->error_message("Failed to write final ace: $final_ace");
	    return;
	}
	unlink $int_file;
	push @new_aces, $final_ace;
    }
    return \@new_aces;
}
#cats multiple ace files together
sub merge_acefiles {
    my ($self, %params) = @_;

    my $int_file = $self->intermediate_file_name('merge');
    unlink $int_file;
    my $int_fh = Genome::Sys->open_file_for_writing( $int_file );
    my $contig_count = 0;     my $read_count = 0;
    #incrementing contigs numbering by 1M for each acefile
    my $increment = ( defined $params{increment} ? $params{increment} : 1000000 );
    my $inc_count = 0;
    foreach my $ace_in (@{$params{ace_files}}) {
	my $fh = Genome::Sys->open_file_for_reading( $ace_in );
	while (my $line = $fh->getline) {
	    if ($line =~ /^AS\s+/) {
		chomp $line;
		my ($c1, $c2) = $self->get_counts_from_ace_AS_line($line);
		$contig_count += $c1;
		$read_count += $c2;
		next;
	    }
	    if ($line =~ /^CO\s+/) {
		chomp $line;
		my $contig_number = $self->get_contig_number_from_ace_CO_line($line);
		#need to rename contigs so there are no duplicate names
		#first ace file will retain same contig numbering
		#following ace files will have contig numbers incremented by 1,000,000
		#TODO - need intelligent way of doing this
                my ($sc_num, $ct_num) = split(/\./, $contig_number);
		my $new_contig_number = (($increment * $inc_count) + $sc_num);
                if ( $ct_num ) {
                    $new_contig_number .= '.'.$ct_num;
                }
		$line =~ /^CO\s+(\S+)/; #just to capture $'
		$int_fh->print("CO Contig".$new_contig_number."$'"."\n");
		next;
	    }
	    last if ($line =~ /^CT{/ or $line =~ /^WA{/); #reached end of contigs .. tags not transferred
	    $int_fh->print($line);
	}
	$fh->close;
	$inc_count++;
    }
    $int_fh->close;
    #returns final ace name but not needed here
    unless ($self->rewrite_ace_file($int_file, $contig_count, $read_count, 'merged.final')){
	$self->error_message("Failed to write final ace file");
	return;
    }
    unlink $int_file;
    return 1;
}

#writes updated ace file with correct contig and read counts
sub rewrite_ace_file {
    my ($self, $ace_in, $contig_count, $read_count, $name) = @_;

    unless (-s $ace_in) {
	$self->error_message("Can't find int ace file or file is zero size: ".$ace_in);
	return;
    }
    my $ace_out = $name.'.ace';
    $ace_out = $self->directory.'/'.$ace_out if $self->directory;
    $self->status_message("Writing final ace file: $ace_out");

    my $fh_in = Genome::Sys->open_file_for_reading($ace_in) ||
	return;

    unlink $ace_out;
    my $fh_out = Genome::Sys->open_file_for_writing($ace_out) ||
	return;

    $fh_out->print("AS $contig_count $read_count\n\n");

    while (my $line = $fh_in->getline) {
	$fh_out->print($line);
    }

    # check for phdball files and create wa tags for them
    my @ball_files;
    if ( @ball_files = $self->_assembly_phdball_files ) {
        $fh_out->print("WA{\nphdBall newbler 080416:144002\n");
        $fh_out->print( map {'../phdball_dir/'.$_."\n"} @ball_files );
        $fh_out->print("}\n");
    }

    $fh_out->close;
    $fh_in->close;

    return $ace_out;
}

# check for phdball dir to write wa tags
sub _assembly_phdball_files {
    my $self = shift;
    my $phdball_dir = $self->directory.'/../phdball_dir';
    my @phdball_files;
    if ( -d $phdball_dir ) {
        my @files = glob( $phdball_dir."/*ball*" );
        for my $file ( @files ) {
            push @phdball_files, File::Basename::basename( $file );
        }
    }
    return @phdball_files;
}

#get ace contig and read counts
sub get_counts_from_ace_AS_line {
    my ($self, $line) = @_;
    $line =~ /^AS\s+(\d+)\s+(\d+)/;
    my $contig_count = $1;
    my $read_count = $2;
    unless ($contig_count =~ /^\d+$/ and $read_count =~ /^\d+$/) {
	$self->error_message("Can't get contig and read counts from ace AS line: $line");
	return;
    }
    return $contig_count, $read_count;
}
#get read count from ace CO line
sub get_read_count_from_ace_CO_line {
    my ($self, $line) = @_;
    my ($read_count) = $line =~ /CO\s+\S+\s+\d+\s+(\d+)/;
    unless ($read_count and $read_count =~ /^\d+$/) {
	$self->error_message("Can't get read count from line: $line");
	return;
    }
    return $read_count;
}
#get contig name from ace CO line
sub get_contig_name_from_ace_CO_line {
    my ($self, $line) = @_;
    my ($contig_name) = $line =~ /CO\s+(\S+)\s+\d+\s+\d+/;
    unless ($contig_name) {
	$self->error_message("Can't get contig name from line: $line");
	return
    }
    return $contig_name;
}
#get contig number from ace CO lines
sub get_contig_number_from_ace_CO_line {
    my ($self, $line) = @_;
    my $contig_name = $self->get_contig_name_from_ace_CO_line($line);
    my ($contig_number) = $contig_name =~ /Contig(\S+)/i;
    unless ($contig_number =~ /^\d+$/ or $contig_number =~ /^\d+\.\d+$/) {
	$self->error_message("Can't get contig number from contig name: $contig_name".
			     "names should look like this: Contig4 or Contig 8.9");
	return;
    }
    #TODO - see what happens with Contig0 or Contig0.0 though these shouldn't exist
    return $contig_number;
}
#name for temp intermediate files that are not functional ace files
sub intermediate_file_name {
    my ($self, $name) = @_;

    return $self->directory.'/'.$name.'.intermediate';
}

1;
