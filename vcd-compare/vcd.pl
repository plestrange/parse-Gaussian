#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum min max);
use POSIX qw(floor);

our $version = "1.5";

=head1 NAME

VCD Compare

=head1 DESCRIPTION

This program sifts through several Gaussian output files from VCD calculations. The 
different conformations are thermodynamically average and combined. This spectrum 
can also be statistically compared to an experimental spectrum.

=begin comment

Version Info

0.6	caps the gamma optimization at 12
0.7	optimizes the scale factor from the user input gamma value
0.8	optimizes the scale factor and then the gamma value
0.9 	????? i don't remember
1.0	adds a check for imaginary frequencies
1.1	uses a coarse and then a fine grid to find the scale factor more quickly
1.2	????? i don't remember
1.3	added more complex data structures(hashes of arrays),changed parsing of data using split function
		lorentzian plot it around the peak instead of over the whole spectrum
1.4	changed search algorithm to search for greatest absolute value of SimVCD
		default temp is now 300 K, The scale search is now done at gamma=8
1.5	checks for normal termination and changes control structure for input

@flag = [e t T r w f g s]

TO DO LIST

WRITE ERROR SUBROUTINE!!!
	Add error for improper flags

=end comment

=cut

if ( $#ARGV >= 0 ) {
    exec("perldoc", $0) if $ARGV[0] eq "-h";
}
our $start_dir = `pwd`;
#chomp($start_dir);
our $theor_dir = $start_dir;
our $exp_dir = $start_dir;
our $write_dir = $start_dir;
our($basis_func,$degfree,$avg_time);
our(@name,@eng,@bolt,@flag,@abs);
our(%freq,%rot,%weight);
our(@theor_freq,@theor_int,%theor_hash);
our(@exp_freq,@exp_int,%exp_hash,$step_size);
our($SimVCD,@Sim,$SimOld,@scale,@fine_Sim,%sim_hash);
our($output);
our $nu_low = 1000;
our $nu_high = 1500;
our $scale = 0.95;
our $g = 8;
our $temp = 300;

=head1 DEFAULTS

Command line arguments influences the types of input accepted and certain default parameters.
By default, the program reads and writes to the current directory. It expects every .out and 
.log file within the current directory to be a Gaussian VCD output file containing information 
about one of the conformations of a single compound. It expects experimental data to be 
within the current directory. The default temperature is 300 K and the default range of 
interest is 1000-1500 cm^-1. The program will optimize the scale factory and the peak width
of the theoretical data to best fit the experimental spectrum. 

=over 5

=item B<-e>

Only experimental data will be accepted. The data will be trimmed to the range of interest
and the peak intensities normalized to be between -1 and 1. 

=item B<-f>

Specify the frequency range of interest. Default is 1000-1500 cm^-1.

=item B<-g>

Switch off optimization of the peak width of the lorentzian functions to best match experiment. 

=item B<-h>

Prints help page and shows documentation for program.

=item B<-r>

Specify where the program will look for the input files.

=item B<-s> 

Switch off optimization of the frequency scale factor to best match experiment.

=item B<-t>

Only Gaussian output files will be accepted. The different spectra will be thermodynamically
weighted and combined to form one spectrum. 

=item B<-T>

Specify temperature for thermodynamic weighting. Default is 300K.  

=item B<-w> 

Specify where the program will write the output file.

=back

=head1 VERSION

1.5

=head1 AUTHOR

Patrick J Lestrange, E<lt>plestran@uw.eduE<gt>

=cut

##################
## Main Program ##
##################

#set default values for flags
$flag[0] = 1;		#experimental data
$flag[1] = 1;		#theoretical data
$flag[2] = 0;		#change temperature
$flag[3] = 0;		#change read directory
$flag[4] = 0;		#change write directory
$flag[5] = 0;		#change freq range
$flag[6] = 1;		#optimize gamma
$flag[7] = 1;		#optimize scale factor

if ( scalar @ARGV >= 0 ) {
    parse(@ARGV);
}

# check frequency flag
if ($flag[5]) {
	print "What is the range of the spectrum?\n";
	print "Low frequency (default: 1000 cm^-1) ";
	$nu_low = <STDIN>;
	chomp($nu_low);
	if ($nu_low eq "") {
		$nu_low = 1000;
	}
	print "High frequency (default: 1500 cm^-1) ";
	$nu_high = <STDIN>;
	chomp($nu_high);
	if ($nu_high eq "") {
		$nu_high = 1500;
	}
}
# check temperature flag
if ($flag[2]) {
	print "What is the temperature for this model? ";
	$temp = <STDIN>;
	chomp($temp);
}
# check read directory flag
if ($flag[3]) {
	if ($flag[0]) {
		print "Where is the experimental data located? (default: current directory) ";
		$exp_dir = <STDIN>;
		chomp($exp_dir);
		if ($exp_dir eq "") {
			$exp_dir = `pwd`;
		}
	}
	if ($flag[1]) {
		print	"Where are the Gaussian output files located? (default: current directory) ";
		$theor_dir = <STDIN>;
		chomp($theor_dir);
		if ($theor_dir eq "") {
			$theor_dir = `pwd`;
		}
	}
}
# check write directory flag	
if ($flag[4]) {
	print "Where should the output file be written? (default: current directory) ";
	$write_dir = <STDIN>;
	chomp($write_dir);
	if ($write_dir eq "") {
		$write_dir = `pwd`;
	}
}
# check exp data flag
if ($flag[0]) {
	search_exp_file();
	normalize();
}
# check theoretical data flag
if ($flag[1]) {
	search_logfiles();
	bolt_weight();
#	check whether to optimize
	if ($flag[0]) {
		which_opt();
	}
	else{
		$step_size = 2;
		print "What is the scale factor for the modeled spectra? (default: 0.98) ";
		$scale = <STDIN>;
			print "What is the gamma value for the lorentzian peak widths? (default: 8) ";
		$g = <STDIN>;
		chomp($scale, $g);
		if ($scale eq "") {
			$scale = 0.98;
		}
		if ($g eq "") {
			$g = 8;
		}
		lorentz_func();
		normalize();
	}
}
# write the output file
write_output();

#################
## Subroutines ##
#################

sub bolt_weight{

#	Constants used for calculation
	my $k 		= 1.3806488E-23;		#in J/K
	my $conv		= 4.359744344E-18;	#in J/au
	my $kT		= $temp*$k/$conv;

	my @q;
	my $q = 0;
	my $min = min @eng;
	
	for(my $i = 0; $i < scalar @eng; $i++) {
		$q[$i] = exp(-($eng[$i]-$min)/$kT);
		$q += $q[$i];
	}
	for(my $i = 0; $i < scalar @eng; $i++) {
		$bolt[$i] = $q[$i]/$q;
	}
	foreach my $key (keys %rot){
		for(my $j = 0; $j < $degfree; $j++) {
			$weight{$key}[$j] = $rot{$key}[$j]*$bolt[$key];
		}
	}
}

sub error{
	if ($_[0] == 1) {
		print "$_[1] contains imaginary frequencies.\n";
		print "It is a transition structure.\n";
		print "Other files may also be transition structures.\n";
	}
	elsif ($_[0] == 2) {
		print "File $_[1] did not terminate properly.\n";
		print "Other files may also not have terminated properly.\n";
	}
	elsif ($_[0] == 3) {
		print "There are no input files.\n";
		print "Check the directory specified and ensure that there are Gaussian .log or .out files.\n"; 
	}
	die;	
}

sub file_test{

	while (1) {
		print "What file should I write to? (default: output.txt) ";
		$output = <STDIN>;
		chomp $output;
		if ($output eq "") {
			$output = "output.txt";
		}
		if (-d $output) {
			print "No, $output is a directory.\n";
			next;
		}
		if (-e $output) {
			print "File already exists. What should I do?\n";
			print "(Enter 'r' to write to a different name, ";
			print "'o' to overwrite or";
			print "'b' to back up to $output.old)\n";
			my $choice = <STDIN>;
			chomp $choice;
			if ($choice eq "r") {
				next;
			} elsif ($choice eq "o") {
				unless (-o $output) {
					print "Can't overwrite $output, it's not yours.\n";
					next;
				}
				unless (-w $output) {
					print "Can't overwrite $output: $!\n";
            		next;
         		}
			} elsif ($choice eq "b") {
				if ( rename($output, $output.".old") ) {
        	    	print "OK, moved $output to $output.old\n";
        	 	} else {
        	    	print "Couldn't rename file: $!\n";
        	    	next;
        	 	}
			} else {
				print "I didn't understand that answer.\n";
				next;
      			}
		}
		last if open(OUT, '>', $output);
		print "I couldn't write to $output: $!\n";
	}
}

sub gamma_opt{

#	Search until Old value is better than new value
	$SimOld = 0;
	my $go = 1;
	while ($go) {				

		lorentz_func();
		normalize();
		SimVCD();
	
		if (abs($SimOld) > abs($SimVCD)) {
			$g = $g - 0.5;
			lorentz_func();
			normalize();
			SimVCD();	
			$go = 0;
		}
		else {
			$SimOld = $SimVCD;
			$g = $g + 0.5;
		}
		
		if($g > 12) {
			$g = $g - 0.5;
			lorentz_func();
			normalize();
			SimVCD();
			$go = 0;
		}
	}
}

sub lorentz_func{

	my @lorentz;			
	my $pi = 3.14159265;	

#	Forms a 4001 part array full of 0's
	for(my $i = 0; $i <= 4000; $i++) {
		$theor_int[$i] = 0;
		$lorentz[$i] = 0;
	}

	foreach my $key (keys %freq){
		for(my $j = 0; $j < $degfree; $j++) {
			my $cm = floor $freq{$key}[$j];
#			Plots peak 2000 units around rot. str.
			if ($cm > 1000) {
				for(my $i = ($cm-1000); $i <= ($cm+1000); $i++) {
					$lorentz[$i] = $weight{$key}[$j]/$pi;
					$lorentz[$i] = $lorentz[$i]*($g);
					$lorentz[$i] = $lorentz[$i]/((($i*$step_size/$scale)-$freq{$key}[$j])**2+($g)**2);
				}
			}	
			else {
				for(my $i = 0; $i <= ($cm+2000); $i++) {
					$lorentz[$i] = $weight{$key}[$j]/$pi;
					$lorentz[$i] = $lorentz[$i]*($g);
					$lorentz[$i] = $lorentz[$i]/((($i*$step_size/$scale)-$freq{$key}[$j])**2+($g)**2);
				}	
			}
			for(my $i = 0; $i <= 4000; $i++) {
				$theor_int[$i] = $theor_int[$i] + $lorentz[$i];	
			}
		}
	}
#	Scales all the frequencies by the scalar $step_size
	for(my $i = 0; $i <= 4000; $i++) {
		$theor_freq[$i] = $i*$step_size;
	}
#	Pairs the frequency and the intensity
	for(my $i = 0; $i < scalar @theor_freq; $i++) {
		$theor_hash{$theor_freq[$i]} = $theor_int[$i];
	}
}

sub normalize{

	my (@trash, @sort, $max);	
	if ($flag[1]) {	
#		Trims frequency array to range of interest
		for(my $i = 0; $i < scalar @theor_freq; $i++) {
			if ($theor_freq[$i] >= $nu_low && $theor_freq[$i] <= $nu_high) {
				@trash = (@trash, $theor_freq[$i]);
			}
		}
		@theor_freq = @trash;
		@trash = ();
	
#		Trims intensity array to frequency range of interest
		for(my $i = 0; $i < scalar @theor_freq; $i++) {
			@trash = (@trash, $theor_hash{$theor_freq[$i]});
		}
		@theor_int = @trash;
		@trash = ();

#		Finds the maximum intensity with which to normalize the peaks
		for(my $i = 0; $i < scalar @theor_int; $i++) {
			$sort[$i] = abs($theor_int[$i]);
		}
		$max = max @sort;
		for(my $i = 0; $i < scalar @theor_int; $i++) {
			$theor_int[$i] = $theor_int[$i]/$max;
		}
	}
#	The same is done if there is experimental data	
	if ($flag[0]) {	
		for(my $i = 0; $i < scalar @exp_freq; $i++) {
			if ($exp_freq[$i] >= $nu_low && $exp_freq[$i] <= $nu_high) {
				@trash = (@trash, $exp_freq[$i]);
			}
		}
		@exp_freq = @trash;
		@trash = ();
	
		for(my $i = 0; $i < scalar @exp_freq; $i++) {
			@trash = (@trash, $exp_hash{$exp_freq[$i]});
		}
		@exp_int = @trash;
		@trash = ();
	
		@sort = ();
		for(my $i = 0; $i < scalar @exp_int; $i++) {
			$sort[$i] = abs($exp_int[$i]);
		}
			
		$max = max @sort;
		for(my $i = 0; $i < scalar @exp_int; $i++) {
			$exp_int[$i] = $exp_int[$i]/$max;
		}
	}
}

sub parse{
	my @inp = @_;
	for (my $i = 0; $i < scalar @inp; $i++) {
		$_ = $inp[$i];
		if ($_ =~ /-t/){
			$flag[0] = 0;
			next;
		}
		if ($_ =~ /-e/){
			$flag[1] = 0;
			next;
		}
		if ($_ =~ /-T/){
			$flag[2] = 1;
			next;
		}
		if ($_ =~ /-r/){
			$flag[3] = 1;
			next;
		}
		if ($_ =~ /-w/){
			$flag[4] = 1;
			next;
		}
		if ($_ =~ /-f/){
			$flag[5] = 1;
			next;
		}
		if ($_ =~ /-g/){
			$flag[6] = 0;
			next;
		}
		if ($_ =~ /-s/){
			$flag[7] = 0;
			next;
		}
	}
}

sub scale_opt{
	
#	Coarse scale factor search over 0.95-1.03
	my $go = 1;
	while ($go) {

		lorentz_func();
		normalize();
		SimVCD();	
		
		push (@scale, $scale);
		push (@Sim, $SimVCD);
		$scale = $scale + 0.005;	
		$go = 0 if $scale > 1.03;
	}
	
	for(my $i = 0; $i < scalar @Sim; $i++) {
		$abs[$i] = abs($Sim[$i]);
	}
	for(my $i = 0; $i < scalar @Sim; $i++) {
		$sim_hash{$abs[$i]} = $scale[$i];
	}

	my $max = max @abs;
	my $coarse_scale = $sim_hash{$max};	
	$scale = $coarse_scale - 0.004;	
	@scale = ();
	@Sim = ();
	@abs = ();
	
#	Fine search around best match from coarse search
	$go = 1;
	while ($go) {
		
		lorentz_func();
		normalize();
		SimVCD();	
	
		push (@scale, $scale);
		push (@Sim, $SimVCD);
		$scale = $scale + 0.001;
				
		if($scale > $coarse_scale + 0.004) {
			$go = 0;
		}
	}
	
	%sim_hash = ();
	for(my $i = 0; $i < scalar @Sim; $i++) {
		$abs[$i] = abs($Sim[$i]);
	}
	for(my $i = 0; $i < scalar @Sim; $i++) {
		$sim_hash{$abs[$i]} = $scale[$i];
	}

	$max = max @abs;
	$scale = $sim_hash{$max};
}

sub search_exp_file{

	if ($exp_dir ne $start_dir) {
		chdir $exp_dir or die "Can't find exp_dir directory: $exp_dir \n"; 
	}
	print "What is the file name of the experimental data? ";
	my $exp_input = <STDIN>;
	chomp($exp_input);

#	Searches the input file for the frequency and intensity
	open(IN, '<', $exp_input) or die "Can't read input file: $!\n"; 
	while (<IN>) { 
		my @string = split /	/, $_;
		last if $string[0] > 3000;
		push (@exp_freq, $string[0]);
		push (@exp_int, $string[1]);
	}
	chomp(@exp_freq, @exp_int);
	for(my $i = 0; $i < scalar @exp_int; $i++) {
		$exp_hash{$exp_freq[$i]} = $exp_int[$i];
	}
#	Defines the step size between frequencies
	$step_size = $exp_freq[1] - $exp_freq[0];
	if ($exp_dir ne $start_dir) {
		system('pwd');
#		chdir($start_dir) or die "Can't find start_dir directory: $start_dir \n"; 
		chdir $start_dir;
		system('pwd');
		#system('cd $start_dir'); 
	}
}

sub search_logfiles{

	my $i = 0;				
	my (@days, @hours, @minutes, @seconds, @time);

	if ($theor_dir ne $start_dir) {	
		chdir $theor_dir or die "Can't find directory $theor_dir \n";
	}
	foreach (glob('*.log *.LOG *.out *.OUT')) { 
		
		my @string;
		my $normterm = 0;	
		$name[$i] = $_;
		
		open(IN, '<', $_) or die "Can't read input file: $!\n"; 
		while(<IN>){
#			Finds the number of basis functions			
			if ($_ =~ /primitive/) {
				@string = split / +/;
				$basis_func = $string[1];
			}	
#			Finds the degrees of freedom in the log file
			if ($_ =~ /Deg. of freedom/) {
				@string = split / +/;
				$degfree = $string[4];
			}
#			Finds the rotational strengths in the log file
			if ($_ =~ /Rot. str.   --/) {
				@string = split / +/;
				chomp @string;
				push (@{$rot{$i}}, $string[4], $string[5], $string[6]);
			}
#			Finds the frequencies in the log file
			if ($_ =~ /Frequencies --/) {
				@string = split / +/;
				chomp @string;
#				Checks for imaginary frequencies		
				if ($string[3] < 0 or $string[4] < 0 or $string[5] < 0) {
					error(1, $name[$i]);
				}
				push (@{$freq{$i}}, $string[3], $string[4], $string[5]);
			}
#			Finds the energy of the conformer
			if ($_ =~ /Sum of electronic and thermal Free Energies=/) {
				@string = split / +/;
				push (@eng, $string[8]); 
			}
#			Finds the CPU time in log file
			if ($_ =~ /Job cpu time/) {
				@string = split / +/;
				push (@days, $string[4]);
				push (@hours, $string[6]);
				push (@minutes, $string[8]);
				push (@seconds, $string[10]);
			}
#		Checks for normal termination of job
		$normterm = 1 if $_ =~ /Normal termination/;
		}
		if ($normterm == 0) { error(2, $name[$i]); }
		$normterm = 0;
		$i++;
	}
#	Computes the CPU time in hours
	if ($#days > 0) {
		for (my $i = 0; $i < scalar @days; $i++) {
			$days[$i] = $days[$i]*24;
			$minutes[$i] = $minutes[$i]/60;
			$seconds[$i] = $seconds[$i]/3600;
			$time[$i] = $days[$i]+$hours[$i]+$minutes[$i]+$seconds[$i];
		}
	} else { error(3); }
	$avg_time = sum(@time)/@time;  
	chomp($basis_func,$degfree,@eng,$avg_time);	
	if ($theor_dir ne $start_dir) {	
		chdir $start_dir or die "Can't find directory: $start_dir \n"; 
	}
}

sub SimVCD{

	my $Icc = 0;
	my $Ioo = 0;
	my $Ico = 0;

	for(my $i = 0; $i < scalar @exp_freq-1; $i++) {
		$Icc = $Icc + ($theor_int[$i]*$theor_int[$i]);
		$Ioo = $Ioo + ($exp_int[$i]*$exp_int[$i]);
		$Ico = $Ico + ($exp_int[$i]*$theor_int[$i]);
	}

	$SimVCD = $Ico/($Icc + $Ioo - abs($Ico));
	
	print "SimVCD (scale = $scale, gamma = $g): 	$SimVCD\n";
}

sub which_opt{

	$SimOld = 0;
	if ($flag[6] && $flag[7]){
		$g 			= 4;
		$scale 		= 0.95;
		scale_opt();
		gamma_opt();	
	}
	elsif ($flag[6] && $flag[7] eq 0){
		print "What is the scale factor for the modeled spectra? (default: 0.98) ";
		$scale = <STDIN>;
		chomp($scale);
		if ($scale eq "") {
			$scale = 0.98;
		}
		$g = 4;
		gamma_opt();
	}
	elsif ($flag[6] eq 0 && $flag[7]){
		print "What is the gamma value for the lorentzian peak widths? (default: 8) ";
		$g = <STDIN>;
		chomp($g);
		if ($g eq "") {
			$g = 8;
		}
		$scale = 0.95;
		scale_opt();
	}
	else {
		print "What is the scale factor for the modeled spectra? (default: 0.98) ";
		$scale = <STDIN>;
		print "What is the gamma value for the lorentzian peak widths? (default: 8) ";
		$g = <STDIN>;
		chomp($scale, $g);
		if ($scale eq "") {
			$scale = 0.98;
		}
		if ($g eq "") {
			$g = 8;
		}
		lorentz_func();
		normalize();
		SimVCD();
	}
}

sub write_output{

	if ($write_dir ne $start_dir) {
		system('pwd');	
		chdir $write_dir; 
		#chdir($write_dir) or die "Can't find directory: $write_dir \n";
		system('pwd');	
	}
	file_test();

	open(OUT, '>', $output);
	if ($flag[1]) {	
		for(my $i = 0; $i < scalar @name; $i++) {
			print OUT "Name of file:\t$name[$i]\t\t\t"
		}
		print OUT "\n";
		for(my $i = 0; $i < scalar @name; $i++) {
			print OUT "Energy in Hartrees:\t$eng[$i]\t\t\t";
		}
		print OUT "\n";
		for(my $i = 0; $i < scalar @name; $i++) {
			print OUT "Boltzmann Weight:\t$bolt[$i]\t\t\t";
		}
		print OUT "\n";
		for(my $i = 0; $i < scalar @name; $i++) {
			print OUT "Frequency\tRot.Str.\tWeighted Rot.Str.\t\t";
		}
		print OUT "\n";
		for (my $j = 0; $j < $degfree; $j++) {
			for(my $i = 0; $i < scalar @name; $i++) {
				print OUT "$freq{$i}[$j]\t$rot{$i}[$j]\t$weight{$i}[$j]\t\t";
			}
			print OUT "\n";
		}
		print OUT "\n\nScalie Factor\t$scale\n";
		print OUT "Gamma\t$g\n";
		print OUT "Avg CPU Time (hrs)\t$avg_time\n";
		print OUT "Basis Functions\t$basis_func\n";
		if ($flag[0]) {
			print OUT "SimVCD\t$SimVCD\n";
			print OUT "Frequency\tExperiment\tTheory\n";
			for (my $i = 0; $i < scalar @theor_freq; $i++) {
				print OUT "$exp_freq[$i]\t$exp_int[$i]\t$theor_int[$i]\n";
			}
		} else {
			print OUT "Frequency\tTheory\n";
			for (my $i = 0; $i < scalar @theor_freq; $i++) {
				print OUT "$theor_freq[$i]\t$theor_int[$i]\n";
			}
		}	
	} 
	if ($flag[0]) {
		print OUT "Frequency\tExp Int\n";
		for(my $i = 0; $i < scalar @exp_freq; $i++) {
			print OUT "$exp_freq[$i]\t$exp_int[$i]\n";
		}
	} 
	close OUT; 
	print "Output file $output made in directory: $write_dir \n";
	if ($write_dir ne $start_dir) {	
		system('pwd');	
		#chdir($start_dir) or die "Can't find directory: $write_dir \n";
		chdir $start_dir;
		system('pwd');	
	}
}

