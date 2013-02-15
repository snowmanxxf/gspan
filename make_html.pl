#!/usr/bin/perl -w
use strict;
use warnings;

# Warning, this directory and files will be totally overwriten
my $html_image_dir = "images";
my $html_file = "frequent_graphs.html";

my $CLOSEGRAPH_CMD="closegraph_st";

sub process_file
{
    my $img_prefix = $_[1];

    my %dot_codes;
    my %comments;
    my %img_files;

    my $graph_id;
    my $num_graphs = 0;
    my $file_comment = "";

    open FILE, "<", $_[0] or die $!;
    while (my $line = <FILE>)
    {
	my @values = split(' ', $line);
	if (@values)
	{
	    if ($values[0] =~ "^t.*")
	    {
		$dot_codes{$graph_id} = $dot_codes{$graph_id} . "};\n\n" if ($num_graphs > 0);
		
		$graph_id = $values[2];
		$dot_codes{$graph_id} = "";
		$comments{$graph_id} = "";
		$dot_codes{$graph_id} = $dot_codes{$graph_id} . "graph" . " " . $graph_id . " {\n";
		$num_graphs += 1;
	    }
	    elsif ($values[0] =~ "^v.*")
	    {
		$dot_codes{$graph_id} = $dot_codes{$graph_id} . "\t" . $values[1] . " [label=\"" . $values[2] . "\"]\n";
	    }
	    elsif ($values[0] =~ "^e.*")
	    {
		my $s;
		if ($values[3])
		{
		    $s = "\t" . $values[1] . " -- " . $values[2] . " [label=\"" . $values[3] . "\"]\n";
		}
		else
		{
		    $s = "\t" . $values[1] . " -- " . $values[2] . "\n";
		}
		$dot_codes{$graph_id} = $dot_codes{$graph_id} . $s;
	    }
	    elsif ($values[0] =~ "^#")
	    {
		if (defined($graph_id))
		{
		    #print "|" . $line . "!";
		    $comments{$graph_id} = $comments{$graph_id} . $line;
		}
		else
		{
		    if ($num_graphs == 0)
		    {
			$file_comment = $line;
		    }
		}
	    }
	}
    }
    $dot_codes{$graph_id} = $dot_codes{$graph_id} . "};\n\n" if ($num_graphs > 0);
    close(FILE);
    
    my $dot;
    while (($graph_id, $dot) = each %dot_codes)
    {
	open DOT_FILE, ">", "/tmp/dotfile";
	print DOT_FILE $dot;
	close DOT_FILE;
	my $result_file = "$html_image_dir/" . "$img_prefix" . "$graph_id.png";
	system("echo -n make file: $result_file ... ; neato -Tpng -o $result_file /tmp/dotfile && echo ' done'; rm -f /tmp/dotfile");

	print HTML_FILE "<object data=\"$result_file\" type=\"image/png\"> </object>\n";
	print HTML_FILE "<br>" . "graph " . $graph_id . "<br />\n";

	$comments{$graph_id} =~ s/\n/<br>/g;
	$comments{$graph_id} =~ s/#mv/<br>#mv/g;
	print HTML_FILE "<br>" . $comments{$graph_id} . "<br />\n";
    }
}


my $input_file  = $ARGV[0];
my $minsupport  = $ARGV[1];
my $tr_filename;

if ($input_file =~ ".*.lg")
{
    $tr_filename = $input_file;
}
elsif ($input_file =~ ".*.tgf")
{
    $tr_filename = "/tmp/tr_graphs.lg";
    system("./tgf_to_lg.pl <$input_file > $tr_filename");
}

open HTML_FILE, ">", $html_file;
 
# --------------------------------
# process input (transactional) graphs
# --------------------------------
system("rm -rf $html_image_dir; mkdir $html_image_dir");
print HTML_FILE "<html>\n";
print HTML_FILE "<head>\n";
print HTML_FILE "<title>Frequent Graphs</title>\n";
print HTML_FILE "</head>\n";
print HTML_FILE "<body>\n";
print HTML_FILE "<h2> INPUT GRAPHS </h2>\n";

process_file $tr_filename, "graph_";

# --------------------------------
# run closegraph $minsupport < $tr_filename > pattern_file
# --------------------------------
system("cat $tr_filename | ./$CLOSEGRAPH_CMD $minsupport > /tmp/patterns.lg");

# --------------------------------
# process pattern_file
# --------------------------------
print HTML_FILE "<h2> FOUND PATTERNS, minimum support=$minsupport </h2>\n";
process_file "/tmp/patterns.lg", "pattern_";

print HTML_FILE "</body>\n";
print HTML_FILE "</html>\n";
close(HTML_FILE);
