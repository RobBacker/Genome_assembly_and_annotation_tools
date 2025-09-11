#!/usr/bin/env python3

import csv
import argparse
from pathlib import Path
from Bio import SeqIO
import re

"""
Integrates EnTAP and DeTEnGA results into a structural annotation (GFF3) and optional FASTAs.
Adds TE classification and functional annotations (PFAM, IPR, EggNOG, UniProt) into GFF3 `Note` 
and `Dbxref` fields, and into FASTA headers.

Usage:
    python annotate_gff_entap_detenga.py \
        --entap annotated.tsv \
        --detenga TE_summary.csv \
        --gff annotation.gff3 \
        cds.fa proteins.fa ...

Inputs:
    --entap     EnTAP annotation table (annotated.tsv)
    --detenga   DeTEnGA TE summary file (*.csv, semicolon-delimited)
    --gff       Structural annotation (GFF3)
    fastas      (optional) one or more FASTA files to annotate

Outputs:
    - Annotated GFF3 (suffix: .annotated.gff3)
    - Annotated FASTA(s) (suffix: .annotated.fa)
    Author: Robert Backer
"""

# === Annotation header block ===
ANNOTATION_HEADER_BLOCK = """
# Annotation keys used in Note= field:
#   TE_status: DeTEnGA TE classification of the transcript.
#   function: Top BLAST or DIAMOND match description (SeqSearch_Description)
#   PFAM: Pipe-separated Pfam domain identifiers (Database_UniProt_Protein_Domains)
#   IPR: InterPro protein descriptions (Database_InterProScan_Protein_Description)
#   EggNOG_desc: Functional description from EggNOG (Database_EggNOG_Description)
#   Dbxref: Cross-references (InterPro, UniProt, KEGG, EggNOG)
"""

def clean(val):
    v = val.strip()
    return '' if not v or v.upper() in ('NA','NAN') else v


def is_valid_eggnog_id(og):
    for prefix in ('COG','KOG','PTHR','ENOG','NOG','OG'):
        if og.upper().startswith(prefix) and '@' not in og:
            return True
    return False


def parse_entap_tsv(entap_file):
    entap_data = {}
    with open(entap_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            tid = row['Query_Sequence']
            func = clean(row.get('SeqSearch_Description',''))
            egg_desc = clean(row.get('Database_EggNOG_Description',''))

            # UniProt accession
            uni_acc = ''
            subj = row.get('SeqSearch_Subject_Sequence','')
            if subj and '|' in subj:
                parts = subj.split('|')
                if len(parts) > 1:
                    uni_acc = parts[1]

            # PFAM domains from UniProt_Protein_Domains
            pfam_raw = clean(row.get('Database_UniProt_Protein_Domains',''))
            pfam = ''
            if pfam_raw:
                pfam = pfam_raw.replace(';','|').replace(',', '|').strip('|')

            # InterPro IDs (with optional names) from InterProScan_InterPro_ID
            ipr_raw = clean(row.get('Database_InterProScan_InterPro_ID',''))
            ipr_with_names = []
            if ipr_raw:
                for entry in re.split(r'[;,|]', ipr_raw):
                    entry = entry.strip()
                    if entry:
                        ipr_with_names.append(entry)

            # InterPro description for Note
            ipr_desc = clean(row.get('Database_InterProScan_Protein_Description',''))

            # KEGG KOs from EggNOG_KO
            kegg_raw = clean(row.get('Database_EggNOG_KEGG_KO',''))
            ko_ids = [k.strip().replace('ko:','') 
                      for k in kegg_raw.replace('|',',').split(',') 
                      if k.strip().upper().startswith('K')]

            # EggNOG OGs
            og_raw = clean(row.get('Database_EggNOG_Member_OGs',''))
            ogs = [og.strip() for og in og_raw.replace('|',',').split(',') 
                   if is_valid_eggnog_id(og.strip())]

            entap_data[tid] = {
                'function': func,
                'PFAM': pfam,
                'IPR_desc': ipr_desc,
                'InterPro_with_names': ipr_with_names,
                'EggNOG_desc': egg_desc,
                'UniProt_acc': uni_acc,
                'KEGG_KOs': ko_ids,
                'EggNOG_OGs': ogs,
            }
    return entap_data


def parse_detenga_tsv(detenga_file):
    te_data = {}
    with open(detenga_file, newline='') as f:
        reader = csv.DictReader(f, delimiter=';')
        for row in reader:
            tid = row['Transcript_ID']
            status = clean(row.get('Interpro_status',''))
            if status:
                te_data[tid] = status
    return te_data


def parse_attributes(attr_str):
    d = {}
    for field in attr_str.split(';'):
        if '=' in field:
            k, v = field.split('=',1)
            d[k] = v
    return d


def format_attributes(d):
    return ';'.join(f"{k}={v}" for k,v in d.items())


def annotate_gff(in_gff, entap, te):
    p = Path(in_gff)
    out = p.with_name(p.stem + ".annotated" + p.suffix)
    with open(in_gff) as fin, open(out,'w') as fout:
        lines = fin.readlines()
        hdr = max(i for i, ln in enumerate(lines) if ln.startswith('#'))
        fout.writelines(lines[:hdr+1])
        fout.write(ANNOTATION_HEADER_BLOCK + "\n")
        for ln in lines[hdr+1:]:
            if ln.startswith('#') or not ln.strip():
                fout.write(ln)
                continue
            cols = ln.rstrip('\n').split('\t')
            if len(cols) != 9:
                fout.write(ln)
                continue
            attrs = parse_attributes(cols[8])
            tid = attrs.get('ID') or attrs.get('Parent')
            if not tid:
                fout.write(ln)
                continue

            # Build Note=
            notes = []
            if tid in te:
                notes.append(f"TE_status={te[tid]}")
            e = entap.get(tid, {})
            for tag, val in (('function', e.get('function')),
                             ('PFAM', e.get('PFAM')),
                             ('IPR', e.get('IPR_desc')),
                             ('EggNOG_desc', e.get('EggNOG_desc'))):
                if val:
                    notes.append(f"{tag}={val}")
            if notes:
                existing = attrs.get('Note', '')
                attrs['Note'] = existing + '|' + '|'.join(notes) if existing else '|'.join(notes)

            # Build Dbxref=
            dbx = []
            for iprn in e.get('InterPro_with_names', []):
                ipr_id = iprn.split('(')[0].strip()
                if ipr_id.startswith("IPR"):
                    dbx.append(f"InterPro:{ipr_id}")
            if e.get('UniProt_acc'):
                dbx.append(f"UniProt:{e['UniProt_acc']}")
            dbx += [f"KEGG:{k}" for k in e.get('KEGG_KOs', [])]
            dbx += [f"EggNOG:{og}" for og in e.get('EggNOG_OGs', [])]
            if dbx:
                existing = attrs.get('Dbxref', '')
                attrs['Dbxref'] = existing + ',' + ','.join(dbx) if existing else ','.join(dbx)

            cols[8] = format_attributes(attrs)
            fout.write('\t'.join(cols) + "\n")
    print(f"Annotated GFF written to: {out}")


def annotate_fasta(fasta, entap):
    p = Path(fasta)
    out = p.with_name(p.stem + ".annotated" + p.suffix)
    with open(fasta) as fin, open(out,'w') as fout:
        for rec in SeqIO.parse(fin, 'fasta'):
            e = entap.get(rec.id, {})
            parts = []
            for tag, val in (('function', e.get('function')),
                             ('PFAM', e.get('PFAM')),
                             ('IPR', e.get('IPR_desc')),
                             ('EggNOG_desc', e.get('EggNOG_desc'))):
                if val:
                    parts.append(f"{tag}={val}")
            rec.description = f"{rec.id} " + " ".join(parts) if parts else rec.id
            SeqIO.write(rec, fout, 'fasta')
    print(f"Annotated FASTA written to: {out}")


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--entap', required=True)
    p.add_argument('--detenga', required=True)
    p.add_argument('--gff', required=True)
    p.add_argument('fastas', nargs='*')
    args = p.parse_args()

    entap = parse_entap_tsv(args.entap)
    te = parse_detenga_tsv(args.detenga)
    annotate_gff(args.gff, entap, te)
    for f in args.fastas:
        annotate_fasta(f, entap)

if __name__ == '__main__':
    main()


