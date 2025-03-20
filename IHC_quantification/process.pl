#!/usr/bin/perl

my $LOCK_UN = 8;
my $LOCK_EX = 2;
my $LOCK_SH = 1;
my $LOCK_NB = 4;

my $TRUE = 1;
my $FALSE = 0;

#==================
#Error
#==================
sub quit_error {
my $error = shift;
print qq~Error: $error\n~;
exit 0;
} 

#==================
#Main
#==================
my $argv_1 = shift( @ARGV );
my $argv_2 = shift( @ARGV );
my @data = ();
my @proc = ();
my $j = 0;

my $r,$g,$b;
my @output = (0,0,0,0);

#read in RGB array data
open(FILE, "<$argv_1") || quit_error("Error with $argv_1.csv!"); 
flock(FILE, $LOCK_EX); 

while( <FILE> ) {
	#print "$j\n";
	chomp;
		$data[$j] = [ split( /\,/, $_ ) ];
		#print "$data[$j][1]\n";
		$j++;
}
flock(FILE, $LOCK_UN);
close(FILE);
#print "value of j;$j";
#;print "number=" . scalar(@data[1]) . "\n";

#process array data
for( my $i = 0; $i < $j; $i++ ) {
	#for ( my $i = 0; $i < scalar(@data[$j]; $i++ ) {
	#	$data[$j][$i] /= 255;
	#}
	$data[$i][1] /= 255; #normalize green
	$data[$i][0] *= $data[$i][1] / 255; #normalize red to green
	$data[$i][2] *= $data[$i][1] / 255; #normalize blue to green

	if  ( $data[$i][0] > $data[$i][2] ) {#r>b
		$output[1] += $data[$i][0] - $data[$i][2];
		$output[3] += $data[$i][2];
		$output[0] += $data[$i][1] - $data[$i][0];
	} else {
		$output[2] += $data[$i][2] - $data[$i][0];
		$output[3] += $data[$i][0];
		$output[0] += $data[$i][1] - $data[$i][2];
	}
}

#file out
open(FILE, ">>$argv_2") || quit_error("Could not output file: $argv_2.\n");
flock(FILE, $LOCK_EX);
print FILE "$argv_1";
for (my $i = 0; $i < scalar(@output); $i++) {
	print FILE "," . $output[$i];
}
print FILE "\n";
flock(FILE, $LOCK_UN);
close(FILE);
unlink("$argv_1");
