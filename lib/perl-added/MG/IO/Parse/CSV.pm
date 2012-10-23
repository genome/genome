#####################################################################################
# perl module used to parse in CSV file format
#*** MG::IO::Parse::CSV.pm ***#
# Copyright (C) 
#    2007 Brian Dunford-Shore bshore@watson.wustl.edu
# All rights reserved!
#####################################################################################

# This parses the CSV files and produces a data structure.
# IN THE FUTURE:

package MG::IO::Parse::CSV;
#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
use warnings;
use strict;

use Data::Dumper;
use Text::CSV_XS;

#__SET THE STYLE OF OUTPUT WHEN USING THE SCRIPT IN "CHECK" MODE
$Data::Dumper::Indent = 1;

sub new {
  my ($class, %arg) = @_;
  my $self = {
							_source => $arg{source} || 'CSV',
							_keyfields => $arg{keyfields} || '',
							_tabkey => $arg{tabkey} || 0,
							_processor => $arg{processor} || {} ,#|| MG::?????->new(),
							_header_translation => $arg{header_translation} || {},
							_no_header => $arg{no_header} || 0,
							_header_skip => $arg{header_skip} || 0,
							_header_fields => $arg{header_fields} || undef,
							_ucheader_fields => $arg{ucheader_fields} || undef,
							_field_subset => $arg{field_subset} || undef,
							_record => $arg{record} || {},
                            _separator => $arg{separator} || undef,
                            _no_spaces => $arg{no_spaces} || 0,
						 };
  bless($self, $class || ref($class));

  return $self;
}

#  =========
#  parse_CSV  parse a csv file format type, i.e. a "csv" file
#  =========  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub Parse {
	my ($self,$fh,$file,%args) = @_;
	my $source = $self->{_source};
	my $processor = $self->{_processor};
	my $keyfields = $self->{_keyfields};
	my $field_subset_array = $self->{_field_subset};
	my $header_translation = $self->{_header_translation};
	my $header_skip = $self->{_header_skip};
	my $no_header = $self->{_no_header};
	my $header_fields = $self->{_header_fields};
	my $ucheader_fields = $self->{_ucheader_fields};
    my $separator = $self->{_separator};
    my $no_spaces = $self->{_no_spaces};
	my $line_number_field = 'file_line_num';
	my $line_field = 'file_line';
	my $line_num = 1;
    my $record = $self->{_record};
    my $csv = defined($separator) ? Text::CSV_XS->new({'sep_char' => $separator}) : Text::CSV_XS->new();
        

	for (my $skip = 0; $skip < $header_skip; $skip++) {
		my $skip_line = <$fh>;
		$line_num++;
	}
	if ((defined($no_header) && $no_header) &&
			!defined($header_fields)) {
		warn "The no_header argument was given but not the header_fields argument";
		return $record;
	}
	my $header;
	unless (defined($no_header) && $no_header) {
		$header = <$fh>;
		$line_num++;
	}
	if (defined($header_fields) && $header_fields) {
		$header = $header_fields;
	}
    unless(defined($separator)) {
        $header =~ s/\t/,/gx;
    }
	$csv->parse($header);
	my @header_fields = $csv->fields();
	if ($ucheader_fields) {
		@header_fields = map { uc($_) } @header_fields;
	}
    if(defined($no_spaces) && $no_spaces) {
        @header_fields = map { $_ =~ s/ /_/g } @header_fields;
    }
	unshift( @header_fields, $line_field);		# Add 'extra' fields of the input line
	unshift( @header_fields, $line_number_field);	# and the line number
	# Translate the header names, if a translation is given
	for (my $h = 0; $h <= $#header_fields; $h++) {
		$header_fields[$h] = (exists($header_translation->{$header_fields[$h]})) ?
			$header_translation->{$header_fields[$h]} : $header_fields[$h];
	}

	$keyfields ||= $line_number_field; # Default to the line number as the key field
	my (@key_fields) = split(':',$keyfields);
	my %key_fields;
	@key_fields{ @key_fields } = ( 0 .. $#key_fields );

	# Construct field name to position lookup
	my %header_fields;
	@header_fields{ @header_fields } = ( 0 .. $#header_fields );
	# Construct a subset of the fields--the default is the complete set of fields
	my %field_subset;
	if (defined($field_subset_array)) {
		@field_subset{ @{$field_subset_array} } = @header_fields{ @{$field_subset_array } };
	} else {
		@field_subset{ @header_fields } = ( 0 .. $#header_fields );
	}
	# Construct a list of fields that are not key (are values only)
	my @value_fields;
	foreach my $field (@header_fields[ (values %field_subset ) ]) {
		unless (exists($key_fields{$field})) {
			push @value_fields, ($field);
		}
	}
	#__PARSE FILE
	my $line;
	while (<$fh>) {
		chomp;
		$line = $_;
        my $temp = $_;
        unless(defined($separator)) {
            #maintain original default behavior of handling both tabs and
            #commas    
            $temp =~ s/\t/,/gx;
        }
        $_ = $temp;
        $csv->parse($_);
		my @values = $csv->fields();

        if(defined($no_spaces) && $no_spaces) {
            @values = map { $_ =~ s/ /_/g } @values;
        }
		unshift (@values, $line);		# Add 'extra' fields of the input line
		unshift (@values, $line_num++);	# and the line number

		my $sub_record;
		if ($self->{_tabkey}) {
			# Construct the 'flat' tab delimited key
			my $key = join("\t", @values[ @header_fields{@key_fields} ]);
			unless (exists($record->{ $key } )) {
				$record->{$key} = {};
			}
			$sub_record = $record->{$key};
		} else {
			# Construct the hierarchical key structure
			$sub_record = $record;
			foreach my $sub_key  (@key_fields) {
				unless (exists($sub_record->{ $values[ $header_fields{ $sub_key } ] } )) {
					$sub_record->{ $values[ $header_fields{ $sub_key } ] } = {};
				}
				$sub_record = $sub_record->{ $values[ $header_fields{ $sub_key } ] };
			}
		}
		# Get the hash array slice of the values
		@{$sub_record}{ @value_fields } = @values[ @header_fields{@value_fields} ];

    unless (defined($args{all}) && $args{all}) {
			#__DUMP PARSED RESULTS AND STOP IF JUST CHECKING
			if (defined($args{check}) && $args{check}) {
				print "DATA:\n";
				print Dumper ($record);
				$record = {};
			} else {
				unless (defined($args{no_process}) && $args{no_process}) {
					# Call processor (database loader)
					$processor->Process($record,%args);
					$record = {};
				}
			}
		}
	}

	if (defined($args{all}) && $args{all}) {
		#__DUMP PARSED RESULTS AND STOP IF JUST CHECKING
		if (defined($args{check}) && $args{check}) {
			print "DATA:\n";
			print Dumper ($record);
		} else {
			unless (defined($args{no_process}) && $args{no_process}) {
				# Call processor (database loader)
				$processor->Process($record,%args);
			}
		}
	}

	#__STOP IF JUST CHECKING
	if (defined($args{check}) && $args{check}) {
		exit;
	}

	#__RETURN DATA STRUCTS
	$self->{_record} = $record;
	return ($record);
}

=head1 NAME

 MG::IO::Parse::CSV -- parses the CSV file format.

=head1 SYNOPSIS

 use MG::IO::Parse::CSV;

 # Create a new parser object and process
 my $source = 'CSV';
 my $parser = MG::IO::Parse::CSV->new(source => $source);
 my %args = (
            check => 1,
# Define all to one to have everything processed into a single structure (and then processed)
#            all => 1, 
# Define no_process to one to have everything processed into a single structure (and not processed)
#            no_process => 1, 
            );
 # The args are also passed through to the processor (database loader)
 #__OPEN INPUT FILE FOR PARSING
 use FileHandle;
 my $fh = new FileHandle;
 unless ($fh->open (qq{$file})) {
    die "Could not open CSV file '$file' for reading";
 }
 $parser->Parse($fh,$file,%args);
 #__CLOSE INPUT FILE
 $fh->close;

 exit 1;

=head1 DESCRIPTION

This module parses the CSV file format.

=head1 Constructor and Initialization

=item (object) new (arguments)

=head2 arguments (optional)

=item source

 Sets the source.

=item keyfields

 Identifies the key fields.  The default key field, if none is given, is the implicit field 'file_line_num'.  This uses the field names given (or overridden) in the header fields.

=item tabkey

 Construct a 'flat' tab delimited key.  The default is to construct a hierarchical key.

=item processor

 Specify a parser to post process the data.  These can be chained.  See the MG::Transform::Process class for the base class definition.

=item header_translation

 Specify a hash for translating the field names.  Mostly useful if using the header within the file but just overriding the names of some fields.

=item no_header

 Specify that there is no header in the file.  The 'header_fields' argument is mandatory if the no_header argument is given.

=item header_skip

 Skip the given 'header_skip' number of lines before starting to process the file.  If you need to parse the data in these lines, then don't specify this argument--instead, process these lines prior to calling this class.

=item header_fields

 Specify, using a comma delimited list, the names of the fields in the file.  This may be used even if the file contains a header either because the header may be wrong or the header names are not optimal.

=item field_subset

 Specify the subset of fields to process.  This uses the field names in the header of the 'header_fields' argument.

=item record

 The record data structure to add data to.  This allows passing a previous data structure to a new instance.  This record data structure is added to and returned on each call to the parser.

=item separator

 The separator character of the records in the file. If none is given defaults
 to assuming fields are separated by commas or tabs. This should be specified
 explicitly if you expect that your fields may contain either of these
 characters.

=item no_spaces

 Replace any spaces within fields with an underscore. By default this is off.
 
=head2 Methods

=item (object) Parse (arguments)

=head2 arguments

=item file handle

 Open input file handle to parse.  You may parse up to the header prior to call the Parse method or, if you specify the 'header_fields' and the 'no_header' argument, up to any point in the file prior to calling the parser.

=item file

 The name of the file that is the open file handle.  This is used for error reporting.

=item args

 Optional hash array arguments.  These include: 'all' to parse all records prior to processing, 'check' to no call the processor but merely to 'check' or dump them, and 'no_process' to not call the processor.

=head2 EXPORT

None by default.

=head1 SEE ALSO

See also MG::Transform::Process.  See MG::Transform::Process::MutationCSV for an example of usage.

=head1 FILES

None

=head1 BUGS

I'm sure that you'll find some. Let me know.

=head1 AUTHOR

Brian Dunford-Shore, E<lt>bshore@watson.wustl.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2007 by Brian Dunford-Shore.  All Rights Reserved.

=cut

__END__

1;

# $Header$
