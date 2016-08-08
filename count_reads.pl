#!/usr/bin/perl

$bam = $ARGV[0];
$bai = $ARGV[1];
$genelist = $ARGV[2];
$output = $ARGV[3];
$wsize = $ARGV[4];
$interval = $ARGV[5];
$sliding_num = $ARGV[6];
$divide = $ARGV[7];
$multiply = $ARGV[8];
$output_bedgraph = $ARGV[9];

open(GENELIST, $genelist);
open(OUTPUT, "> $output");

my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime();
my $tmp = $year,$mon,$mday,$hour,$min,$sec;

system "ln -s $bam $tmp";
system "ln -s $bai $tmp.bai";
system "mkdir bedgraph";

# Get number of total reads
system "echo $divide";
if ($divide == 1) {
  system "echo $bam";
  $stat = `samtools flagstat $tmp`;
  if ($stat =~ /([0-9]+) \+ [0-9]+ in total .*/) { $reads_total = $1; }
  print "Total Reads: $reads_total\n";
}

while ($line = <GENELIST>) {
  print $line;
  chomp($line);
  ($txstart, $chr, $strand, $name, $ids) = split(/\t+/, $line); 
  print OUTPUT "$chr\t$strand\t$txstart\t$name\t$ids";
  open(OUTPUT2, "> bedgraph/$chr:$txstart$strand.bedgraph");
  if ($strand eq '+') { $sign = 1; } else { $sign = -1; }

  # $rstart and $rend are used for executing samtools view,
  # so $start is always smaller than $end whichever $strand is.
  $rstart = $txstart - $interval * $sliding_num - ($wsize / 2);
  $rend   = $txstart + $interval * $sliding_num + ($wsize / 2);
  print "$chr\t$strand\t$txstart\t$name\t$rstart\t$rend\n";
  
  # List of start position of reads
  @reads_start = ();
  my $error_text = "";
  open(CMD, "samtools view $tmp $chr:$rstart-$rend 2>&1 |");
  while (<CMD>) {
    if ($_ =~ /chr[0-9a-zA-Z]+\t([0-9]+)/) { push(@reads_start, $1); }
    elsif ($_ =~ /main_samview/) {
      $error_text = $_;
      last;
    }
  }
  close(CMD);

  # Error
  if ($error_text ne "") {
    print "$chr\t$strand\t$txstart\t$name\t$rstart\t$rend\t$error_text";
    print OUTPUT "\t$error_text\n";
    next;
  }

  # For all windows
  for ($i = -$sliding_num; $i <= $sliding_num; $i++) {
    $wstart = $txstart + ($i * $interval * $sign) - ($wsize / 2);
    $wend = $wstart + $wsize;
    
    # Number of reads whose start position is in this window
    my $reads_count = 0;
    foreach $start (@reads_start) {
      if ($start >= $wstart && $start < $wend) { $reads_count++; }
    }

    # Per million reads
    if ($divide) {
      $value = $reads_count * $multiply / $reads_total;
    } else {
      $value = $reads_count * $multiply;
      
      # $bstart and $bend are used to create bedgraph
      $bstart = $txstart + ($i * $interval * $sign) - ($interval / 2);
      $bend = $bstart + $interval;
      print OUTPUT2 "$chr\t$bstart\t$bend\t$reads_count\n"
    }
    print "  $i\t$wstart\t$wend\t$reads_count\t$reads_million\n";
    print OUTPUT "\t$value";
  }
  print OUTPUT "\n";
  close (OUTPUT2);
}

close (GENELIST);
close (OUTPUT);

system "rm $tmp";
system "rm $tmp.bai";
system "tar czf $output_bedgraph bedgraph";
