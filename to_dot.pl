#!/usr/bin/perl -w
use strict;
use warnings;

my $directed = 0;
my $num_graphs = 0;
my $graph_id;

my $graph_kind;
my $edgeop;

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
	}
	elsif ($values[0] =~ "^e.*")
	{
	    print "\t" . $values[1] . $edgeop . $values[2] . " [label=\"" . $values[3] . "\"]\n";;
	}
	elsif ($values[0] =~ "#undirected")
	{
	    $directed = 0;
	    $graph_kind = "graph";
	    $edgeop  = " -- ";
	}
	elsif ($values[0] =~ "#directed")
	{
	    $directed = 1;
	    $graph_kind = "digraph";
	    $edgeop  = " -> ";
	}
	elsif ($values[0] =~ "#found_in")
	{
	    print "graph [comment=\"GRAPH " . $graph_id . "\"]\n"
	}
    }
}


print "};\n" if ($num_graphs > 0);

