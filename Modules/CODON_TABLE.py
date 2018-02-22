#!/usr/bin/env python

CODON_TABLE = {"TTT":{'3letter':'Phe','1letter':'F','full':'Phenylalanine'},
"TTC":{'3letter':'Phe','1letter':'F','full':'Phenylalanine'},

"TCT":{'3letter':'Ser','1letter':'S','full':'Serine'},
"TCC":{'3letter':'Ser','1letter':'S','full':'Serine'},
"TCA":{'3letter':'Ser','1letter':'S','full':'Serine'},
"TCG":{'3letter':'Ser','1letter':'S','full':'Serine'},

"TTA":{'3letter':'Leu','1letter':'L','full':'Leucine'},
"TTG":{'3letter':'Leu','1letter':'L','full':'Leucine'},
"CTT":{'3letter':'Leu','1letter':'L','full':'Leucine'},
"CTC":{'3letter':'Leu','1letter':'L','full':'Leucine'},
"CTA":{'3letter':'Leu','1letter':'L','full':'Leucine'},
"CTG":{'3letter':'Leu','1letter':'L','full':'Leucine'},

"ATT":{'3letter':'Ile','1letter':'I','full':'Isoleucine'},
"ATC":{'3letter':'Ile','1letter':'I','full':'Isoleucine'},
"ATA":{'3letter':'Ile','1letter':'I','full':'Isoleucine'},

"ATG":{'3letter':'Met','1letter':'M','full':'Methionine'},

"GTT":{'3letter':'Val','1letter':'V','full':'Valine'}, 
"GTC":{'3letter':'Val','1letter':'V','full':'Valine'},
"GTA":{'3letter':'Val','1letter':'V','full':'Valine'},
"GTG":{'3letter':'Val','1letter':'V','full':'Valine'},

"CCT":{'3letter':'Pro','1letter':'P','full':'Proline'},
"CCC":{'3letter':'Pro','1letter':'P','full':'Proline'},
"CCA":{'3letter':'Pro','1letter':'P','full':'Proline'},
"CCG":{'3letter':'Pro','1letter':'P','full':'Proline'},

"ACT":{'3letter':'Thr','1letter':'T','full':'Threonine'},
"ACC":{'3letter':'Thr','1letter':'T','full':'Threonine'},
"ACA":{'3letter':'Thr','1letter':'T','full':'Threonine'},
"ACG":{'3letter':'Thr','1letter':'T','full':'Threonine'},

"GCT":{'3letter':'Ala','1letter':'A','full':'Alanine'},
"GCC":{'3letter':'Ala','1letter':'A','full':'Alanine'},
"GCA":{'3letter':'Ala','1letter':'A','full':'Alanine'},
"GCG":{'3letter':'Ala','1letter':'A','full':'Alanine'},

"TAT":{'3letter':'Tyr','1letter':'Y','full':'Tyrosine'}, 	
"TAC":{'3letter':'Tyr','1letter':'Y','full':'Tyrosine'},

"TAA":{'3letter':'***','1letter':'*','full':'Stop (Ochre)'},

"TAG":{'3letter':'***','1letter':'*','full':'Stop (Amber)'},

"CAT":{'3letter':'His','1letter':'H','full':'Histidine'},
"CAC":{'3letter':'His','1letter':'H','full':'Histidine'},

"CAA":{'3letter':'Gln','1letter':'Q','full':'Glutamine'},
"CAG":{'3letter':'Gln','1letter':'Q','full':'Glutamine'},

"AAT":{'3letter':'Asn','1letter':'N','full':'Asparagine'},
"AAC":{'3letter':'Asn','1letter':'N','full':'Asparagine'},

"AAA":{'3letter':'Lys','1letter':'K','full':'Lysine'},
"AAG":{'3letter':'Lys','1letter':'K','full':'Lysine'},

"GAT":{'3letter':'Asp','1letter':'D','full':'Aspartic acid'},
"GAC":{'3letter':'Asp','1letter':'D','full':'Aspartic acid'},

"GAA":{'3letter':'Glu','1letter':'E','full':'Glutamic acid'},
"GAG":{'3letter':'Glu','1letter':'E','full':'Glutamic acid'},

"TGT":{'3letter':'Cys','1letter':'C','full':'Cysteine'},
"TGC":{'3letter':'Cys','1letter':'C','full':'Cysteine'},

"TGA":{'3letter':'***','1letter':'*','full':'Stop (Opal)'}, 	

"TGG":{'3letter':'Trp','1letter':'W','full':'Tryptophan'},

"CGT":{'3letter':'Arg','1letter':'R','full':'Arginine'}, 	
"CGC":{'3letter':'Arg','1letter':'R','full':'Arginine'}, 	
"CGA":{'3letter':'Arg','1letter':'R','full':'Arginine'}, 	
"CGG":{'3letter':'Arg','1letter':'R','full':'Arginine'},
"AGA":{'3letter':'Arg','1letter':'R','full':'Arginine'}, 	
"AGG":{'3letter':'Arg','1letter':'R','full':'Arginine'}, 	
	
"AGT":{'3letter':'Ser','1letter':'S','full':'Serine'},
"AGC":{'3letter':'Ser','1letter':'S','full':'Serine'}, 	

"GGT":{'3letter':'Gly','1letter':'G','full':'Glycine'},	
"GGC":{'3letter':'Gly','1letter':'G','full':'Glycine'}, 	
"GGA":{'3letter':'Gly','1letter':'G','full':'Glycine'}, 	
"GGG":{'3letter':'Gly','1letter':'G','full':'Glycine'}}	

