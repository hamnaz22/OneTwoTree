equal = "rates=equal"
propinv = "rates=propinv"
gamma = "rates=gamma"
invgamma = "rates=invgamma"
Prset = "Prset applyto=(1)"
Lset = "Lset applyto=(1)"
nst1 = "nst=1"
nst2 = "nst=2"
nst6 = "nst=6"
shapepr = "shapepr=Uniform(0.1,50.0)"
pinvarpr = "pinvarpr=Uniform(0.0,1.0)"
eqfreqpr = "statefreqpr=Fixed(Equal)"
uneqfreqpr = "statefreqpr=Dirichlet(1.0,1.0,1.0,1.0)"
revmatpr = "revmatpr=Dirichlet(1.0,1.0,1.0,1.0,1.0,1.0)"
tratiopr = "tratiopr=Beta(1.0,1.0)"

def mb_block_dict(x):
	return {"JC":   	" " + Lset + " " + nst1 + " " + equal + ";\n"  	+ " " + Prset + " " + eqfreqpr + " " + ";",
			"JC+I": 	" " + Lset + " " + nst1 + " " + propinv+ ";\n" 	+ " " + Prset + " " + eqfreqpr + " " + pinvarpr + ";",
			"JC+G": 	" " + Lset + " " + nst1 + " " + gamma+ ";\n"   	+ " " + Prset + " " + eqfreqpr + " " + shapepr + ";",
			"JC+I+G": 	" " + Lset + " " + nst1 + " " + invgamma+ ";\n" 	+ " " + Prset + " " + eqfreqpr + " " + shapepr + " " + pinvarpr + ";",
			"F81": 		" " + Lset + " " + nst1 + " " + equal+ ";\n" 		+ " " + Prset + " " + uneqfreqpr + ";",
			"F81+I": 	" " + Lset + " " + nst1 + " " + propinv+ ";\n" 	+ " " + Prset + " " + uneqfreqpr + " " + pinvarpr + ";",
			"F81+G": 	" " + Lset + " " + nst1 + " " + gamma+ ";\n" 		+ " " + Prset + " " + uneqfreqpr + " " + shapepr + ";",
			"F81+I+G": 	" " + Lset + " " + nst1 + " " + invgamma+ ";\n" 	+ " " + Prset + " " + uneqfreqpr + " " + shapepr + " " + pinvarpr + ";",
			"K80": 		" " + Lset + " " + nst2 + " " + equal+ ";\n" 		+ " " + Prset + " " + eqfreqpr + " " + tratiopr + ";",
			"K80+I": 	" " + Lset + " " + nst2 + " " + propinv+ ";\n" 		+ " " + Prset + " " + eqfreqpr + " " + tratiopr + " " + pinvarpr + ";",
			"K80+G": 	" " + Lset + " " + nst2 + " " + gamma+ ";\n" 		+ " " + Prset + " " + eqfreqpr + " " + tratiopr + " " + shapepr + ";",
			"K80+I+G": 	" " + Lset + " " + nst2 + " " + invgamma+ ";\n" 	+ " " + Prset + " " + eqfreqpr + " " + tratiopr + " " + shapepr + " " + pinvarpr + ";",
			"HKY": 		" " + Lset + " " + nst2 + " " + equal + ";\n" 	+ " " + Prset + " " + uneqfreqpr + " " + tratiopr + " " + ";",
			"HKY+I": 	" " + Lset + " " + nst2 + " " + propinv + ";\n" 	+ " " + Prset + " " + uneqfreqpr + " " + tratiopr + " " + pinvarpr + ";",
			"HKY+G": 	" " + Lset + " " + nst2 + " " + gamma + ";\n" 	+ " " + Prset + " " + uneqfreqpr + " " + tratiopr + " " + shapepr + ";",
			"HKY+I+G": 	" " + Lset + " " + nst2 + " " + invgamma + ";\n" 	+ " " + Prset + " " + uneqfreqpr + " " + tratiopr + " " + shapepr + " " + pinvarpr + ";",
			"SYM": 		" " + Lset + " " + nst6 + " " + equal + ";\n" 	+ " " + Prset + " " + revmatpr + " " + eqfreqpr + ";",
			"SYM+I": 	" " + Lset + " " + nst6 + " " + propinv + ";\n" 	+ " " + Prset + " " + revmatpr + " " + eqfreqpr + " " + pinvarpr + ";",
			"SYM+G": 	" " + Lset + " " + nst6 + " " + gamma + ";\n" 	+ " " + Prset + " " + revmatpr + " " + eqfreqpr + ";",
			"SYM+I+G": 	" " + Lset + " " + nst6 + " " + invgamma + ";\n" 	+ " " + Prset + " " + revmatpr + " " + eqfreqpr + " " + shapepr + " " + pinvarpr + ";",
			"GTR":   	" " + Lset + " " + nst6 + " " + equal + ";\n" 	+ " " + Prset + " " + revmatpr + " " + uneqfreqpr + ";",
			"GTR+I": 	" " + Lset + " " + nst6 + " " + propinv + ";\n" 	+ " " + Prset + " " + revmatpr + " " + uneqfreqpr + " " + pinvarpr + ";",
			"GTR+G":	" " + Lset + " " + nst6 + " " + gamma + ";\n" 	+ " " + Prset + " " + revmatpr + " " + uneqfreqpr + " " + shapepr + ";",
			"GTR+I+G":	" " + Lset + " " + nst6 + " " + invgamma + ";\n" 	+ " " + Prset + " " + revmatpr + " " + uneqfreqpr + " " + shapepr + " " + pinvarpr + ";"
	}[x]
