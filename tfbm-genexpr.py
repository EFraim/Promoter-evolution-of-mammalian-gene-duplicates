import pickle
import numpy as np
from row_matchers import one_to_one_matches
from pandas import concat
export=pickle.load(open('/data/db/import/save/mouse-export.pkl', 'rb'))

branch_order = ["Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Gnathostomata", "Euteleostomi", "Sarcopterygii", "Tetrapoda", "Amniota", "Mammalia", "Theria", "Eutheria", "Boreoeutheria", "Euarchontoglires", "Glires", "Rodentia", "Myomorpha", "Muroidea", "Muridae", "Murinae", "Mus", "Mus musculus strain reference (CL57BL6)"]
exprcols=[c for c in export.columns if c.endswith('-tpm')]
exprnorm = export.assign(avgexpr = np.log(export[exprcols].mean(axis=1)))


concat(concat([less, full]) for less, full in (one_to_one_matches(exprnorm[(exprnorm['Branch point'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 0)], exprnorm[(exprnorm['Branch point'] == b)&(exprnorm['avgexpr'].notna()) & (exprnorm['CpG-ness'] == 2)], 'avgexpr', 0.1) for b in branch_order))
			   
