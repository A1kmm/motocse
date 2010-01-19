all: motocse which_different
motocse: motocse.hs
	ghc --make motocse.hs
which_different: which_different.hs
	ghc --make which_different.hs