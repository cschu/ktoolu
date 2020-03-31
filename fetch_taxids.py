#!/usr/bin/env python
# coding=utf-8
import sys
import csv

from Bio import Entrez

from ktoolu_taxonomy import ktTaxonomyTree 

whoami = None
assert whoami is not None
Entrez.email = whoami

def fetchTaxonomyData(ids, db='nuccore'):
    handle = Entrez.efetch(db=db, id=','.join(ids), retmode='xml')#, rettype='uilist')
    records = Entrez.read(handle)
    return records

def fetchTaxGILink(ids):
    handle = Entrez.elink(dbfrom='protein', db="Taxonomy", id=','.join(ids), retmode='xml')
    records = Entrez.read(handle)
    return records




transcripts = dict()
for row in csv.reader(sys.stdin, delimiter='\t'):
    acc = row[2].split('|')[3]
    transcripts[acc] = transcripts.get(acc, list())
    transcripts[acc].append(row[1])

"""
[{'GBQualifier_value': 'Alistipes sp. CAG:831', 'GBQualifier_name': 'organism'}, {'GBQualifier_value': 'derived from human gut metagenome', 'GBQualifier_name': 'isolation_source'}, {'GBQualifier_value': 'taxon:1262698', 'GBQualifier_name': 'db_xref'}, {'GBQualifier_name': 'environmental_sample'}, {'GBQualifier_value': 'contig 18', 'GBQualifier_name': 'note'}]


    tax = list(set([GBQual['GBQualifier_value'].split(':')[1] for GBQual in feature['GBFeature_quals'] if GBQual['GBQualifier_value'].startswith('taxon:')]))[0]
KeyError: 'GBQualifier_value'
"""

db = sys.argv[1]

acc_tax_table = dict()
taxData = fetchTaxonomyData(list(transcripts), db=db)
for rec in taxData:
    acc = rec['GBSeq_accession-version']
    
    for feature in rec['GBSeq_feature-table']:
        if feature['GBFeature_key'] == 'source':
            #try:
            tax = list(set([GBQual['GBQualifier_value'].split(':')[1] for GBQual in feature['GBFeature_quals'] if GBQual.get('GBQualifier_value', '').startswith('taxon:')]))[0]
            #except:
            #    # print(rec, file=sys.stderr)
            #    print([GBQual for GBQual in feature['GBFeature_quals']], file=sys.stderr)
            #    sys.exit(1)
            org = list(set([GBQual['GBQualifier_value'] for GBQual in feature['GBFeature_quals'] if GBQual['GBQualifier_name'] == 'organism']))[0]
            break

    acc_tax_table[acc] = acc_tax_table.get(acc, set())
    acc_tax_table[acc].add((tax, org))
    # print(acc, tax, sep='\t')


ktree = ktTaxonomyTree('/Users/schudomc')
branch = ktree.getDescendents(33208)
plants = ktree.getDescendents(33090)
viruses = ktree.getDescendents(10239)
bacteria = ktree.getDescendents(2)
insects = ktree.getDescendents(6960)


for acc in acc_tax_table:
    for tax, org in acc_tax_table[acc]:
        tax = int(tax)
        flag = 'OK' if tax in branch else 'CHECK'
        kingdom = 'other_metazoa' if flag == 'OK' else 'other'
        if kingdom == 'other_metazoa' and tax in insects:
            kingdom = 'insects'
        elif tax in plants:
            kingdom = 'plants'
        elif tax in viruses:
            kingdom = 'viruses'
        elif tax in bacteria:
            kingdom = 'bacteria'
        for transcript in transcripts[acc]:
            print(transcript, acc, tax, org, kingdom, flag, sep='\t')

"""
todo:
- find reasonable ancestor node for insects and other organisms
- flag all transcripts not assigned to the branch rooted at that node as "TO_BE_CHECKED"
    



"""
"""
rec = fetchTaxonomyData(list(transcripts)[:5])[0]
print(rec)
for key in rec.keys(): print(key, rec[key])
acc = rec['GBSeq_accession-version']
# tax = [GBQual['GBQualifier_value'].split(':')[1] for GBQual in rec['GBSeq_feature-table']['GBFeature_quals'] if GBQual['GBQualifier_value'].startswith('taxon:')][0]

for feature in rec['GBSeq_feature-table']:
    if feature['GBFeature_key'] == 'source':
        tax = list(set([GBQual['GBQualifier_value'].split(':')[1] for GBQual in feature['GBFeature_quals'] if GBQual['GBQualifier_value'].startswith('taxon:')]))[0]
        break

print(acc, tax)
"""
"""
GBSeq_feature-table [
 {'GBFeature_intervals': [{'GBInterval_to': '76', 'GBInterval_from': '1', 'GBInterval_accession': 'AFK48864.1'}], 
  'GBFeature_key': 'source', 'GBFeature_location': '1..76', 
  'GBFeature_quals': [{'GBQualifier_name': 'organism', 'GBQualifier_value': 'Medicago truncatula'}, {'GBQualifier_name': 'db_xref', 'GBQualifier_value': 'taxon:3880'}, {'GBQualifier_name': 'clone', 'GBQualifier_value': 'JCVI-FLMt-21I20'}]
 }, 
 {'GBFeature_intervals': [{'GBInterval_to': '76', 'GBInterval_from': '1', 'GBInterval_accession': 'AFK48864.1'}], 
  'GBFeature_key': 'Protein', 
  'GBFeature_location': '1..76', 
  'GBFeature_quals': [{'GBQualifier_name': 'product', 'GBQualifier_value': 'unknown'}, {'GBQualifier_name': 'calculated_mol_wt', 'GBQualifier_value': '8655'}]
 }, 
 {'GBFeature_intervals': [{'GBInterval_to': '76', 'GBInterval_from': '1', 'GBInterval_accession': 'AFK48864.1'}], 
  'GBFeature_key': 'CDS', 
  'GBFeature_location': '1..76', 
  'GBFeature_quals': [{'GBQualifier_name': 'coded_by', 'GBQualifier_value': 'BT149070.1:213..443'}]}
]
"""






