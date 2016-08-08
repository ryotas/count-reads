#!/usr/bin/perl

$bam = $ARGV[0];
$bai = $ARGV[1];
$peaklist = $ARGV[2];
$output = $ARGV[3];

open(PEAKLIST, $peaklist);
open(OUTPUT, "> $output");

my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime();
my $tmp = $year,$mon,$mday,$hour,$min,$sec;

system "ln -s $bam $tmp";
system "ln -s $bai $tmp.bai";

# Get number of total reads
print "bam: $bam\n";
print "bai: $bai\n";
print "bed: $peaklist\n";

# Get number of total peak reads
my $total_peak_reads = 0;
while ($line = <PEAKLIST>) {
  chomp($line);
  ($chr, $start, $end, $name, $score) = split(/\t+/, $line);
  print("samtools view $tmp $chr:$start-$end | wc -l\t");
  open(CMD, "samtools view $tmp $chr:$start-$end 2>&1 | wc -l |");
  while (<CMD>) {
    # When error is caught
    if ($_ =~ /main_samview/) {
      $error_text = $_;
      print("error: $error_text\n");
      last;
    }
    # The 4th column is the start position of the read.
    # The start position will be pushed into @reads_start list.
    else {
      my $read_num = $_;
      print("read_num: $read_num");
      $total_peak_reads = $total_peak_reads + $read_num;
    }
  }
  close(CMD); 
}
print("total_reads_num: $total_peak_reads\n");
print OUTPUT $total_peak_reads;
close (PEAKLIST);
close (OUTPUT);

system "rm $tmp";
system "rm $tmp.bai";
