#!/usr/bin/perl -w
use strict;
use warnings;

my $directed = 0;
my $num_graphs = 0;
my $graph_id;

my $graph_kind;
my $edgeop;

my $lab_inner_circle = "^1";
my @inner_circle_nodes;

if (scalar(@ARGV) > 0 && $ARGV[0] =~ "-dir")
{
    $directed = 1;
    $graph_kind = "digraph";
    $edgeop  = " -> ";
}
else
{
    $directed = 0;
    $graph_kind = "graph";
    $edgeop  = " -- ";
}

my $line;
while (defined($line = <STDIN>))
{
    #my $line = $_;
    my @values = split(' ', $line);
    if (@values)
    {
	if ($values[0] =~ "^t.*")
	{
	    print "};\n\n" if ($num_graphs > 0);
	    $graph_id = $values[2];
	    print $graph_kind . " " . $graph_id . " {\n";
	    $num_graphs += 1;
	}
	elsif ($values[0] =~ "^v.*")
	{
	    print "\t" . $values[1] . " [label=\"" . $values[2] . "\"]\n";

	    if ($values[2] =~ $lab_inner_circle)
	    {
		push(@inner_circle_nodes, $values[1]);
	    }
	}
	elsif ($values[0] =~ "^e.*")
	{
	    if ($values[3])
	    {
		print "\t" . $values[1] . $edgeop . $values[2] . " [label=\"" . $values[3] . "\"]\n";
	    }
	    else
	    {
		print "\t" . $values[1] . $edgeop . $values[2] . "\n";
	    }
	}
	elsif ($values[0] =~ "^#")
	{
	    print $line;
	}
    }
}

if (scalar(@inner_circle_nodes) > 0) {
    print "\t" . -1 . " [label=nil,style=invis]\n";
    print "\t graph [root=-1]\n";
    print "\t ranksep=3\n";
    foreach my $v (@inner_circle_nodes) {
	print -1 . " -- " . $v . " [style=invis] \n";
    }
}

print "}\n" if ($num_graphs > 0);

