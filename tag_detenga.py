#!/usr/bin/env python3
import csv, argparse, re
from pathlib import Path

"""
Attach ONLY DeTEnGA TE status to an existing GFF3 (GFF3-safe, no URL-encoding).

- Adds to mRNA/transcript features ONLY:
    * Note : "TE_status: <value>" (appended; deduped)
- Leaves gene/exon/CDS etc. unchanged.
- Preserves ID/Parent exactly as in the input (no encoding/decoding).
"""

# ---------- helpers ----------
def clean(val):
    if val is None: return ''
    v = str(val).strip()
    return '' if v == '' or v.upper() in ('NA','NAN') else v

# sanitize a single attribute VALUE so it won't break GFF3 parsing
def sanitize_value(v: str) -> str:
    if v is None: return ''
    v = str(v)
    v = v.replace('\t',' ').replace('\n',' ').replace('\r',' ')
    v = v.replace(';',' ').replace('=',' ').replace(',',' ')
    v = re.sub(r'\s+', ' ', v).strip()
    return v

def parse_detenga_csv(detenga_file):
    te = {}
    with open(detenga_file, newline='') as f:
        # DeTEnGA export commonly uses ';' delimiter; change to ',' if yours is CSV-comma
        reader = csv.DictReader(f, delimiter=';')
        # Accept a few common column names
        id_keys = ('Transcript_ID','transcript_id','ID','Query_Sequence')
        status_keys = ('Interpro_status','TE_status','te_status','Status')
        for row in reader:
            # find transcript ID
            tid = ''
            for k in id_keys:
                if k in row:
                    tid = clean(row[k]); 
                    if tid: break
            if not tid: 
                continue
            # find status
            status = ''
            for k in status_keys:
                if k in row:
                    status = clean(row[k]); 
                    if status: break
            if tid and status:
                te[tid] = status
    return te

# Attribute parsing that preserves key order and multiplicity, with NO decoding
def parse_attributes(attr_str):
    order = []
    amap = {}  # key -> list of values
    for field in attr_str.strip().split(';'):
        field = field.strip()
        if not field: continue
        k, sep, v = field.partition('=')
        if not sep:
            k, v = field, ''
        k = k.strip()
        v = v.strip()
        if k not in amap:
            order.append(k)
            amap[k] = []
        vals = v.split(',') if v != '' else ['']
        amap[k].extend(vals)
    return order, amap

def format_attributes(order, amap):
    seen = set()
    final_keys = []
    for special in ('ID', 'Parent'):
        if special in amap:
            final_keys.append(special); seen.add(special)
    for k in order:
        if k not in seen:
            final_keys.append(k); seen.add(k)
    for k in amap.keys():
        if k not in seen:
            final_keys.append(k); seen.add(k)

    parts = []
    for k in final_keys:
        vs = [v for v in amap[k] if v is not None]
        if len(vs)==0:
            parts.append(f"{k}=")
        else:
            parts.append(f"{k}=" + ",".join(vs))
    return ';'.join(parts)

def add_note_items(amap, items):
    """Add a list of free-text items to Note= as separate comma entries (dedup)."""
    if not items: return
    if 'Note' not in amap: amap['Note'] = []
    seen = set(amap['Note'])
    for it in items:
        it = sanitize_value(it)
        if it and it not in seen:
            amap['Note'].append(it)
            seen.add(it)

# ---------- core ----------
ANNOTATION_HEADER_BLOCK = (
    "# Annotation keys appended:\n"
    "#   Note : TE_status: <class>\n"
)

def annotate_gff(in_gff, te):
    p = Path(in_gff)
    out = p.with_name(p.stem + ".detenga" + p.suffix)

    with open(in_gff) as fin, open(out, 'w') as fout:
        lines = fin.readlines()

        # keep the initial header block (top-run of comment lines)
        hdr_end = -1
        for i, ln in enumerate(lines):
            if ln.startswith('#'):
                hdr_end = i
            else:
                break

        # original header + our info block
        if hdr_end >= 0:
            fout.writelines(lines[:hdr_end + 1])
        fout.write(ANNOTATION_HEADER_BLOCK)

        # process feature lines
        for ln in lines[hdr_end + 1:]:
            if ln.startswith('#') or not ln.strip():
                fout.write(ln)
                continue

            cols = ln.rstrip('\n').split('\t')
            if len(cols) != 9:
                fout.write(ln)
                continue

            ftype = cols[2].lower()
            # Only annotate transcript-level rows; leave genes/exons/CDS untouched
            if ftype not in ('mrna', 'transcript'):
                fout.write(ln)
                continue

            order, amap = parse_attributes(cols[8])

            # pick transcript identifier: prefer ID, else Parent
            tid = None
            for key in ('ID', 'Parent'):
                if key in amap and len(amap[key]) > 0 and amap[key][0]:
                    tid = amap[key][0]
                    break
            if not tid:
                fout.write(ln)
                continue

            # If TE_status-like custom keys already exist, fold them into Note then drop
            existing_items = []
            for k in ('TE_status','tE_status','te_status','Te_status'):
                if k in amap:
                    for v in amap[k]:
                        v = sanitize_value(v)
                        if v:
                            existing_items.append(f"TE_status: {v}")
                    del amap[k]
            add_note_items(amap, existing_items)

            # DeTEnGA TE class from table
            if tid in te and te[tid]:
                add_note_items(amap, [f"TE_status: {sanitize_value(te[tid])}"])

            # reserialize attributes (no encoding; preserve ID/Parent and order)
            cols[8] = format_attributes(order, amap)
            fout.write('\t'.join(cols) + "\n")

    print(f"Annotated GFF written to: {out}")

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description="Attach ONLY DeTEnGA TE status to GFF3 (Note=TE_status: ...).")
    ap.add_argument('--detenga', required=True, help="DeTEnGA TE_summary.csv (semicolon-delimited by default)")
    ap.add_argument('--gff', required=True, help="Input structural annotation (GFF3)")
    args = ap.parse_args()

    te = parse_detenga_csv(args.detenga)
    annotate_gff(args.gff, te)

if __name__ == '__main__':
    main()
