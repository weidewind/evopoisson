		#my $str = <GULP>; # skipping first ">iteration" line 
		while (<GULP>){
			if ($_=~ /^[^>]/) {next;}
			print ">iteration number $iteration_number\n";
			my $max_depth;
			my @str_array;
			$str = $_; 

			
			my %sums;
			my %hash;
			
			my $test_obs_summ;
			my $test_exp_summ;
			print "str is now $str\n";
			
	#		while ($str =~ /^[^>]/){ 

			
				if ($str =~ /^site/){
					
					unless ($test_obs_summ - $test_exp_summ < 0.0001 && -$test_obs_summ + $test_exp_summ < 0.0001){
						print "summtest failed! $simsite $simnode obssum $test_obs_summ, expsum $test_exp_summ\n";
					}
		
					$test_obs_summ = 0;
					$test_exp_summ = 0;	
						
					my @str_array = split(/\s+/, $str);
					$simsite = $str_array[1]; # careful
					$simnode = $str_array[3];
					$max_depth = $str_array[5];
					$simmutnum = $str_array[7];

						
					print "CONC ".$simsite." site,".$simnode." node,".$max_depth." maxdepth\n";
					$str = <GULP>; #5.02
				}
				
				my @bin_data;
				while (!($str =~ /^site/)){
					@str_array = split(/,/, $str);
					$test_obs_summ += $str_array[1];
					$test_exp_summ += $str_array[2];
					push @bin_data, \@str_array;
					$str = <GULP>;			
				}
				# now $str is site..
				
				#print " Maxdepth $max_depth itnum $iteration_number bin ".$str_array[0]." exp ".$str_array[2]." obs ".$str_array[1]." \n";
				foreach my $md(@maxdepths){
					if ($max_depth > $md){
						foreach my $group_number(0..scalar @groups-1){
							my $label;
							if (exists $label_hashes{$md}[$group_number]) {$label = $label_hashes{$md}[$group_number];}
							else {$label = 0;}
							if ($counter_hashes{$md}[$group_number]{$simnode}){ # change group_hashes to counter_hashes
								my $mutnums = $counter_hashes{$md}[$group_number]{$simnode};
								for (my $i = 0; $i < scalar @{$mutnums}; $i++){
									my $diff = abs($simmutnum - $mutnums->[$i])/$mutnums->[$i];
									if ($diff <= $mutnum_control){
										## ! we need bin loop here, not outside!
										## hash hash hash
										print "CONC "."group number $group_number md $md node name $simnode simmutnum $simmutnum realmutnum".$mutnums->[$i]."\n";
										foreach my $bindat (@bin_data){
											$sums{$md}[$group_number][$label] += $bindat->[1];
											$hash{$md}[$group_number][$label]{$bindat->[0]}[1] += $bindat->[2];
											$hash{$md}[$group_number][$label]{$bindat->[0]}[0] += $bindat->[1];
										}
										splice (@{$mutnums}, $i, 1);
										last; #we cant use one subtree more than once (in the same group)
									}
								}
								if (scalar @{$mutnums} == 0){
									"CONC no more mutnums for group number $group_number md $md node name $simnode, deleting this node\n";
									delete $counter_hashes{$md}[$group_number]{$simnode};
								}
								if (scalar keys %{$counter_hashes{$md}[$group_number]} == 0){
									"CONC no more nodes for group number $group_number md $md, starting new collection\n";
									$label_hashes{$md}[$group_number] += 1; # label corresponds to number of full collections (undef, if there is 0, 1, if we got one. But collection with this label is not full!)
									$counter_hashes{$md}[$group_number] = \%{ dclone($group_hashes{$md}[$group_number]) };
								}

							}
						}
					}
				}

				
	# because we iterate through site earlier $str = <GULP>;
				
	# from while it		}
			

			# maxbins are different for every iteration. Find maximum and use it.

			$maxbin = max($maxbin, $str_array[0]);
print "CONC "."maxbin $maxbin\n";	
#print "CONC "."sum50 $sum50 sum100 $sum100 sum150 $sum150 norm 50 $norm50 norm 100 $norm100 norm 150 $norm150\n";		
			
			foreach my $md(@maxdepths){ 
				foreach my $group_number(0..scalar @groups-1){
					
					print "For ".$group_names[$group_number]." ".$md." we needed ";
					foreach my $snode (keys %{$counter_hashes{$md}[$group_number]}){
						print "$snode :";
						foreach my $mnum (@{$counter_hashes{$md}[$group_number]{$snode}}){
							print "$mnum".",";
						}
						print "; ";
					}
					print "\n";
					
					
					
					
					unless (! exists $label_hashes{$md}[$group_number]){
						print "found ".$label_hashes{$md}[$group_number]." labels for group ".$group_names[$group_number]."\n";
						foreach my $label(0..($label_hashes{$md}[$group_number]-1)){ # last hash (with current label) is uncomplete
							#print "maxdepth $md group number $group_number \n";
							if ($sums{$md}[$group_number][$label] == 0){
								foreach my $bin(1..$maxbin){
									$hash{$md}[$group_number][$label]{$bin}[0] = "NA";
									$hash{$md}[$group_number][$label]{$bin}[1] = "NA";
								}
							}
							else {
								print "CONC "."norm ".$norms{$md}[$group_number]."\n";
								print "CONC "."sum ".$sums{$md}[$group_number][$label]."\n";
								print "CONC "."in hash, bin 12: ".$hash{$md}[$group_number][$label]{12}[0]."\n";
								
								foreach my $bin(1..$maxbin){
									$hash{$md}[$group_number][$label]{$bin}[0] = $hash{$md}[$group_number][$label]{$bin}[0]*$norms{$md}[$group_number]/$sums{$md}[$group_number][$label];
									$hash{$md}[$group_number][$label]{$bin}[1] = $hash{$md}[$group_number][$label]{$bin}[1]*$norms{$md}[$group_number]/$sums{$md}[$group_number][$label];
								}
							}
							
							my $filehandle = $filehandles{$md}{$group_number};
							print "CONC "."going to print something\n";
							foreach my $bin(1..$maxbin){
								print $filehandle $hash{$md}[$group_number][$label]{$bin}[0].",".$hash{$md}[$group_number][$label]{$bin}[1].",";
							}
							print $filehandle "\n";
						}
					}
				}
			}
			

			
			$iteration_number++;
			if ($iteration_number%50 == 0){
				print "Iteration number ".$iteration_number."\n";
			}
		}