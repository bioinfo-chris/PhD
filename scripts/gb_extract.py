from Bio import SeqIO
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description='parser script for retrieving records from a GenBank file')
parser.add_argument('-i', type=str, help='path to input GenBank file')
args = parser.parse_args()
gb_file = Path(args.i)

# header line
print('accession\tvirus\thost\tcollection_date\tsubmission_date\tcountry\tlatlon\treference_title')

# loop through each record
for record in SeqIO.parse(open(gb_file, 'r'), 'genbank'):

	# set all fields initially to 'unknown'; will overwrite if record is present
	accession, reference_title, organism, host, col_date, sub_date, country, latlon = 'unknown','unknown','unknown','unknown','unknown','unknown','unknown','unknown'
	
	# accession
	accession = record.id

	# submission date
	sub_date = record.annotations['date']
	
	# get paper title, if present
	refs = record.annotations['references']
	for ref in refs:
			if ref.title != "Direct Submission":
				reference_title = ref.title
	
	# get source features
	# virus name
	feat = record.features['type'=='source']
	if 'organism' in feat.qualifiers:
		organism = feat.qualifiers['organism']
	
	# host species
	if 'host' in feat.qualifiers:
		host = feat.qualifiers['host']
	
	# isolation date
	if 'collection_date' in feat.qualifiers:
		col_date = feat.qualifiers['collection_date']
	
	# country
	if 'country' in feat.qualifiers:
		country = feat.qualifiers['country']
	
	# latlon
	if 'latlon' in feat.qualifiers:
		latlon = feat.qualifiers['latlon']
	
	# print fields, one line per record iteration
	print(accession,'\t',''.join(organism),'\t',''.join(host),'\t',''.join(col_date),'\t',sub_date,'\t',''.join(country),'\t',''.join(latlon),'\t',reference_title, sep='')

# END
