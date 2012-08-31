package Genome::Model::Tools::Assembly::AceAxe;

use strict;
use warnings;
use Genome;

our $initialized = 0;
sub init_gtk {
    return if $initialized;
    eval {
        require Gtk;
        require Gtk::GladeXML;
        require Gnome;
    };
    die $@ if $@;
    $initialized = 1;
}

my ($no_pcap, $keep_name);

class Genome::Model::Tools::Assembly::AceAxe {
    is => 'Command',
    has => 
    [ 
        #sort_con => 
        #{
        #    type => 'String',
        #    is_optional => 1,
        #    doc => "sort contig number numerically",
        #},
        no_pcap =>
        {
            type => 'Boolean',
            is_optional => 1,
            doc => "Use this flag to determine whether or not ace axe renames the contig when it moves it to a new ace file",    
        },    
        keep_name => 
        {
            type => 'Boolean',
            is_optional => 1,
            doc => "Use this flag to determine whether or not ace axe renames the contig when it moves it to a new ace file"
        },                       
    ]    
};

sub help_brief {
    ""
}

sub help_synopsis { 
    return;
}
sub help_detail {
    return <<EOS 
    ace-axe -help   print this usage\nace_axe -no_pcap  naming contig number by 1,2,3...
    ace-axe -keep-name   not change contig number
    ace-axe -no-pcap naming contig number by concatenating like .1,.2,.3...
EOS
}

my $glade;
my %acehash;
my %selecthash;
my %contigshash;
my @acelist;
sub execute
{
    my ($self) = @_;
    init_gtk() unless $initialized;
    $keep_name = $self->keep_name;
    $no_pcap = $self->no_pcap;
    init Gnome "AceAxe";
    $glade = new Gtk::GladeXML("/gsc/scripts/share/ace_axe/ace_axe.glade");
    $glade->signal_autoconnect_from_package('Genome::Model::Tools::Assembly::AceAxe');
    @acelist = `ls -1rt *.ace*`; 

    my $aceTextBox = $glade->get_widget("AvailAceFiles");
    $aceTextBox->freeze();
    foreach my $acefile (@acelist){
        $aceTextBox->append( $acefile );
    }
    my $optwidth = $aceTextBox->optimal_column_width(0);
    $aceTextBox->set_column_width(0, $optwidth);
    $aceTextBox->thaw();

    my $mainWin = $glade->get_widget("AceAxeMain");
    $mainWin->signal_connect("destroy",sub { gtk_main_quit() });
    $mainWin->show_all();

    Gtk->main();
    return 1;
}


#subs
sub gtk_main_quit {
    print "closing...\n";
    Gtk->main_quit;
}

sub on_searchAceEntry_insert_text {

    my $Ace_clist = $glade->get_widget("AvailAceFiles");
    my $search_entry = $glade->get_widget("searchAceEntry");
    my $search_value = $search_entry->get_text();
    $search_value =~ s/\./\\\./g;
    $search_value =~ s/\*/\.\*/g;
    my $refnum = 0;
    foreach my $ref (@acelist){
	if ($ref =~ /$search_value/){
	     $Ace_clist->select_row( $refnum, 0);
	     $Ace_clist->moveto( $refnum, 0, 1.0, 0.0);
	 }
	$refnum++;
    }

}

sub on_addacebutton_released {

    my $Ace_clist = $glade->get_widget("AvailAceFiles");
    my @selectedAce = $Ace_clist->selection();
    foreach my $ref (@selectedAce)
    {
	$acehash{$ref}{'ace'} = $acelist[$ref];
	#print "$acelist[$ref]";
    }
    &populate_selectedAceContentsCList;
}

sub populate_selectedAceContentsCList {

    my $contents_clist = $glade->get_widget("SelectedAceContigs");
    $contents_clist->clear();
    $contents_clist->freeze();
    foreach my $ref (keys %acehash){
	my $ace = $acehash{$ref}{'ace'};
	my %sorthash=();
	#set delimiter
	my $acenew = $ace;
	chomp $acenew;
	my $ao = Genome::Model::Tools::Pcap::Ace->new(input_file => $acenew, using_db => 1);
	my @contig_names = @{$ao->get_contig_names};
	foreach(@contig_names)
	{
		$contents_clist->append( "$_ : $ace" );	
	}
	#change delimiter back
	$/ = "\n";
	if (%sorthash) {
	    my ($x, $y);
	    foreach my $ctg(sort{($x)=$a=~/g(\S+)/; ($y)=$b=~/g(\S+)/; $x <=> $y }keys(%sorthash) ){
		$contents_clist->append( "$ctg : $sorthash{$ctg}" );
	    }
	}
    }
    my $optwidth = $contents_clist->optimal_column_width(0);
    $contents_clist->set_column_width(0, $optwidth);
    $contents_clist->thaw();
    
}

sub on_removeacebutton_released {

    my $contents_clist = $glade->get_widget("SelectedAceContigs");
    $contents_clist->clear();
    foreach my $ref (keys %acehash){
	delete($acehash{$ref});
    }
    my $Ace_clist = $glade->get_widget("AvailAceFiles");
    $Ace_clist->unselect_all();
}

sub on_addcontigbutton_released {
    my $contig_clist = $glade->get_widget("SelectedAceContigs");
    my @selectedContig = $contig_clist->selection();
    foreach my $ref (@selectedContig)
    {
	 my $text = $contig_clist->get_text( $ref, 0 );
	 $selecthash{$ref}{'ace'} = $text;
    }

    &populate_selectedContigsCList;
}

sub populate_selectedContigsCList {
    my $select_clist = $glade->get_widget("SelectedIndivContigs");
    $select_clist->clear();
    $select_clist->freeze();
    foreach my $ref (keys %selecthash){
	my $selectiondata = $selecthash{$ref}{'ace'};
	$select_clist->append( $selectiondata );
    }
    my $optwidth = $select_clist->optimal_column_width(0);
    $select_clist->set_column_width(0, $optwidth);
    $select_clist->thaw();
    
}

sub on_removecontigbutton_released {

    my $select_clist = $glade->get_widget("SelectedIndivContigs");
    my @selectedContig = $select_clist->selection();
    
    foreach my $ref (@selectedContig){
	delete($selecthash{$ref});
    }
    &populate_selectedContigsCList;
}

sub on_mergebutton_released {

    my $select_clist = $glade->get_widget("SelectedIndivContigs");
    $select_clist->select_all();
    my @selectedContig = $select_clist->selection();
    foreach my $ref (@selectedContig)
    {
	my $text = $select_clist->get_text( $ref, 0 );
	my ($contig, $ace) = split(/\s\:\s/, $text);
	chomp $contig; chomp $ace;
	
	#if ($contig=~/\./) {
	#    $contig=~s/\./\\\./g;
	#}
	push @{$contigshash{$ace}} ,$contig;
    }
    &merge_data;

    foreach my $ref (keys %selecthash){
	delete($selecthash{$ref});
    }
    $select_clist->clear();
}

# Feiyu made some changes here
sub merge_data {

    my $Ace_entry = $glade->get_widget("AceOutput");
    my $ace_output = $Ace_entry->get_text();
	unlink "ace_axe.db";    
    my $ao_writer = Genome::Model::Tools::Pcap::Ace->new(using_db => 1, db_file => "ace_axe.db");
	my $c_counter = 0;
	
    foreach my $ace_file (keys %contigshash){
		
		my $ao = Genome::Model::Tools::Pcap::Ace->new(input_file => $ace_file, using_db => 1 );
		my @contig_names = @{$contigshash{$ace_file}};
		foreach my $contig_name (@contig_names)
		{	
			my $contig = $ao->get_contig($contig_name);
			if (defined $no_pcap) 
			{
		    	$contig_name = "Contig".$c_counter;
				$c_counter++;
			}
			elsif(!(defined $keep_name))
			{ 
				$contig_name .= ".$c_counter";
				$c_counter++;		    	
			}
			$contig->name($contig_name);
			my @tags = @{$contig->tags};
			foreach my $tag (@tags)
			{
				$tag->parent($contig_name);
			}
			$contig->tags(\@tags);
			$ao_writer->add_contig($contig);	
		}				 
   }
   $ao_writer->write_file(output_file => $ace_output);    

}

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/finishing/trunk/ace_axe/ace_axe.pl $
#$Id: ace_axe.pl 9324 2006-08-23 18:01:33Z fdu $
1;
