#!/usr/bin/env bash
# Pangenome Utilization (HAL) + Variant extraction (WAVE-style VCF)
# ------------------------------------------------------------------
# Generic reference script — not optimized for production use.
# Will not run unless the expected naming conventions and directory
# layout are followed exactly (see Config section).
# Requires: bcftools, bedtools, bgzip/tabix, awk, conda, apptainer (cactus), odgi (optional).
# ------------------------------------------------------------------

set -euo pipefail

#################################
# Config
#################################

# Reference genome name as used inside the HAL file
refPrefix="REF_GENOME_ID"

# Assemblies (basenames) that were added to the HAL (used for PAV)
# Example: assemblies=(ASM001 ASM002 ASM003 ASM004 ASM005 ASM006)
assemblies=(ASM001 ASM002 ASM003)

# Input paths
pansn_dir="/path/to/pansn"            # <base>/<base>.(primary|alternate).longest.pansn.gff3 (directory containing assembly gff3 files (single isoform models only))
HAL_FILE="/path/to/pangenome.full.hal"
WAVE_VCF="/path/to/pangenome.wave.vcf.gz"

# Optional graph (for anchoring with ODGI)
ODGI_GRAPH="/path/to/pangenome.og"    # leave as-is if not using anchor_loci_to_graph
threads=16

# Apptainer/Singularity image for cactus/hal
APPTAINER_CACTUS_SIF="/path/to/cactus.sif"

# Conda environments
CONDA_ENV_AGAT="AGAT"
CONDA_ENV_VARIANT="VariantCalling"

# Output roots
HAL_OUT="/path/to/out/hal_pangenome"  # will contain lift/, tmp/, PAV files, etc.
VAR_OUT="${HAL_OUT}/sv"               # variant products

# Locus name prefix (generic)
pangenePrefix="pangene_"

#################################
# Pangenome Utilization (HAL)
#################################

# HAL genome name mapping
# ref = $refPrefix ; primary = "<asm>.1" ; alternate = "<asm>.2"
hal_genome_name() {
  local asm="$1" ctx="$2"
  if [[ "$ctx" == "ref" ]]; then echo "$refPrefix"
  elif [[ "$ctx" == "primary" ]]; then echo "${asm}.1"
  else echo "${asm}.2"
  fi
}

# Normalize BED chrom: keep last token after '#' and strip trailing '#<digits>'
# "ASM001#2#ctg0310#0" -> "ctg0310"
normalize_bed_seqnames_last_token() {
  local inb="$1" outb="$2"
  awk 'BEGIN{FS=OFS="\t"}
       function strip_rank(x){ sub(/#[0-9]+$/, "", x); return x }
       {
         n=split($1,a,/#/); c=a[n]; strip_rank(c);
         $1=c; print
       }' "$inb" > "$outb"
}

# GFF3 -> BED6; put ASM/CTX inside "name" as ID|||ASM|||CTX (HAL-safe)
gff3_to_bed6_with_meta_in_name() {
  local gff="$1" out="$2" asm="$3" ctx="$4"
  awk -v OFS="\t" -v ASM="$asm" -v CTX="$ctx" '
       BEGIN{FS=OFS="\t"}
       /^#/ {next}
       ($3=="mRNA"||$3=="transcript"){
         id="."; if (match($9,/ID=([^;]+)/,m)) id=m[1];
         name=id "|||" ASM "|||" CTX;
         print $1,$4-1,$5,name,0,$7
       }' "$gff" > "$out"
}

# Expand BED6 (name = ID|||ASM|||CTX) -> BED8: chrom start end id score strand asm ctx
expand_bed6_name_to_bed8() {
  local inb="$1" outb="$2"
  awk -v OFS="\t" '
    BEGIN{FS=OFS="\t"}
    {
      n = split($4, t, /\|\|\|/);     # split literal "|||"
      id  = (n>=1 ? t[1] : $4);
      asm = (n>=2 ? t[2] : ".");
      ctx = (n>=3 ? t[3] : ".");
      print $1,$2,$3,id,$5,$6,asm,ctx
    }' "$inb" > "$outb"
}

# halLiftover: source -> reference; input must be BED6
hal_liftover_bed6_to_ref() {
  local hal="$1" src="$2" ref="$3" inb="$4" outb="$5"
  apptainer exec --nv "$APPTAINER_CACTUS_SIF" \
    halLiftover "$hal" "$src" "$inb" "$ref" "$outb"
  awk 'BEGIN{FS=OFS="\t"} ($2<$3)&&($6=="+"||$6=="-")' "$outb" > "${outb}.tmp" && mv "${outb}.tmp" "$outb"
}

# Cluster lifted intervals into loci on the reference (strand-aware),
# then collapse each gene to a single span per (chr,strand,gene_id,asm,ctx).
# Input must be BED8: chrom start end id score strand asm ctx
cluster_lifted_to_loci() {
  local all_ref_bed="$1" loci_bed="$2" members="$3"

  conda deactivate || true
  conda activate "${CONDA_ENV_AGAT}"

  # 0) Sanity sort
  bedtools sort -i "$all_ref_bed" > "${loci_bed}.sorted.bed8"

  # 1) Collapse per gene to a single union span
  awk 'BEGIN{FS=OFS="\t"}
       {
         chr=$1; s=$2; e=$3; gid=$4; str=$6; asm=$7; ctx=$8;
         k=chr SUBSEP str SUBSEP gid SUBSEP asm SUBSEP ctx;
         if(!(k in seen)){minS[k]=s; maxE[k]=e; CHR[k]=chr; STR[k]=str; GID[k]=gid; ASM[k]=asm; CTX[k]=ctx}
         if(s<minS[k]) minS[k]=s;
         if(e>maxE[k]) maxE[k]=e;
         seen[k]=1
       }
       END{
         for(k in seen){
           print CHR[k], minS[k], maxE[k], GID[k], 0, STR[k], ASM[k], CTX[k]
         }
       }' "${loci_bed}.sorted.bed8" \
  | bedtools sort -i - > "${loci_bed}.collapsed_gene_spans.bed8"

  # 2) Cluster those gene spans (strand-aware)
  bedtools cluster -s -i "${loci_bed}.collapsed_gene_spans.bed8" \
  > "${loci_bed}.clust.bed8"   # cols: chrom start end id score strand asm ctx cluster_id

  # 3) Locus BED = union per cluster_id (generic prefix)
  awk -v OFS="\t" -v PFX="${pangenePrefix}" '
       {
         cl=$NF
         if(!(cl in seen)){chr[cl]=$1; s[cl]=$2; e[cl]=$3; st[cl]=$6}
         if($2<s[cl]) s[cl]=$2;
         if($3>e[cl]) e[cl]=$3;
         seen[cl]=1
       }
       END{
         for(cl in seen) print chr[cl], s[cl], e[cl], PFX cl, 0, st[cl]
       }' "${loci_bed}.clust.bed8" \
  | sort -k1,1 -k2,2n > "$loci_bed"

  # 4) Members: one row per (cluster_id × gene_id × asm × ctx) using the collapsed span
  awk -v OFS="\t" -v PFX="${pangenePrefix}" '
       {
         chr=$1; s=$2; e=$3; gid=$4; asm=$7; ctx=$8; cl=$9;
         print PFX cl, chr ":" s "-" e, gid, asm, ctx
       }' "${loci_bed}.clust.bed8" \
  | sort -u -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 > "$members"

  # Cleanup
  rm -f "${loci_bed}.sorted.bed8" "${loci_bed}.collapsed_gene_spans.bed8" "${loci_bed}.clust.bed8"
}

# Build PAV (1/0) and gene-ID map from members.tsv (locus, region, gene, ASM, CTX)
build_pav_and_geneid_maps() {
  local members="$1"
  local outdir; outdir="$(dirname "$members")"
  local pav="$outdir/PAV_matrix.tsv"
  local gid="$outdir/pangene_gene_ids.tsv"

  # columns list file
  {
    printf "%s\n" "${refPrefix}.1"
    for base in "${assemblies[@]}"; do
      printf "%s\n" "${base}.1"
      printf "%s\n" "${base}.2"
    done
  } > "${outdir}/.cols.txt"

  # ---------- PAV (build to tmp, sort rows, keep header first) ----------
  local tmp_pav="$pav.tmp"

  awk -v OFS="\t" -v REF="${refPrefix}" '
    function ctx2hap(c){ return (c=="primary") ? "1" : ((c=="alternate") ? "2" : c) }
    {
      # members.tsv: locus  region  gene  ASM  CTX
      loc=$1; asm=$4; ctx=$5;
      key=(asm==REF ? asm ".1" : asm "." ctx2hap(ctx));
      seen[loc]=1; hit[loc SUBSEP key]=1
    }
    END{
      while((getline c < colsfile)>0){cols[++n]=c}
      printf "pangene_id";
      for(i=1;i<=n;i++) printf OFS "%s", cols[i];
      printf "\n";
      for(loc in seen){
        printf "%s", loc;
        for(i=1;i<=n;i++){
          k=loc SUBSEP cols[i];
          printf OFS ((k in hit)?1:0)
        }
        printf "\n"
      }
    }' colsfile="${outdir}/.cols.txt" "$members" > "$tmp_pav"

  { head -n1 "$tmp_pav"; tail -n +2 "$tmp_pav" | sort -k1,1; } > "$pav"
  rm -f "$tmp_pav"

  # ---------- Gene IDs (same header handling) ----------
  local tmp_gid="$gid.tmp"

  awk -v OFS="\t" -v REF="${refPrefix}" '
    function ctx2hap(c){ return (c=="primary") ? "1" : ((c=="alternate") ? "2" : c) }
    {
      loc=$1; gid=$3; asm=$4; ctx=$5;
      key=(asm==REF ? asm ".1" : asm "." ctx2hap(ctx));
      seen[loc]=1;
      if((loc SUBSEP key) in G){
        if(G[loc SUBSEP key] !~ ("(^|;)" gid "($|;)")) G[loc SUBSEP key]=G[loc SUBSEP key] ";" gid;
      } else {
        G[loc SUBSEP key]=gid;
      }
    }
    END{
      while((getline c < colsfile)>0){cols[++n]=c}
      printf "pangene_id";
      for(i=1;i<=n;i++) printf OFS "%s", cols[i];
      printf "\n";
      for(loc in seen){
        printf "%s", loc;
        for(i=1;i<=n;i++){
          k=loc SUBSEP cols[i];
          printf OFS ((k in G)?G[k]:"-")
        }
        printf "\n"
      }
    }' colsfile="${outdir}/.cols.txt" "$members" > "$tmp_gid"

  { head -n1 "$tmp_gid"; tail -n +2 "$tmp_gid" | sort -k1,1; } > "$gid"
  rm -f "$tmp_gid"

  echo "PAV: $pav"
  echo "Gene IDs: $gid"
}

# Optional: anchor loci to graph via reference path (convert ref chrom -> graph path)
anchor_loci_to_graph() {
  local loci_bed="$1" out_gbed="$2"
  local tmp="${out_gbed}.path.bed6"
  awk -v OFS="\t" -v R="$refPrefix" '{ $1=R"#0#" $1; print }' "$loci_bed" > "$tmp"
  apptainer exec --nv "$APPTAINER_CACTUS_SIF" \
    odgi position -t $threads -i "$ODGI_GRAPH" -b "$tmp" > "$out_gbed"
  rm -f "$tmp"
  echo "Graph-anchored loci: $out_gbed"
}

# Driver
hal_route_build_pangenes() {
  [[ -s "$HAL_FILE" ]] || { echo "ERROR: HAL missing: $HAL_FILE" >&2; return 1; }
  local ref="$refPrefix"

  # --- Reference (optional: no liftover needed if already in ref coords)
  local ref_gff="${pansn_dir}/${ref}.longest.pansn.gff3"
  local ref_bed6="${HAL_OUT}/tmp/${ref}.ref.longest.mrna.bed6"
  local ref_bed6_norm="${HAL_OUT}/tmp/${ref}.ref.longest.mrna.norm.bed6"
  local ref_bed8="${HAL_OUT}/lift/${ref}.ref.longest.mrna.bed8"

  gff3_to_bed6_with_meta_in_name "$ref_gff" "$ref_bed6" "$ref" "ref"
  normalize_bed_seqnames_last_token "$ref_bed6" "$ref_bed6_norm"
  expand_bed6_name_to_bed8 "$ref_bed6_norm" "$ref_bed8"

  # --- Each assembly: BED6 -> normalize -> halLiftover -> expand to BED8
  for base in "${assemblies[@]}"; do
    for ctx in primary alternate; do
      local gff="${pansn_dir}/${base}/${base}.${ctx}.longest.pansn.gff3"
      [[ -s "$gff" ]] || { echo "skip: $gff"; continue; }

      local inbed6="${HAL_OUT}/tmp/${base}.${ctx}.mrna.bed6"
      local inbed6_norm="${HAL_OUT}/tmp/${base}.${ctx}.mrna.norm.bed6"
      local lifted_bed6="${HAL_OUT}/tmp/${base}.${ctx}.to_${ref}.bed6"
      local lifted_bed8="${HAL_OUT}/lift/${base}.${ctx}.to_${ref}.bed8"

      gff3_to_bed6_with_meta_in_name "$gff" "$inbed6" "$base" "$ctx"
      normalize_bed_seqnames_last_token "$inbed6" "$inbed6_norm"
      hal_liftover_bed6_to_ref "$HAL_FILE" "$(hal_genome_name "$base" "$ctx")" "$ref" "$inbed6_norm" "$lifted_bed6"
      expand_bed6_name_to_bed8 "$lifted_bed6" "$lifted_bed8"
    done
  done

  # --- Cluster on reference -> loci + members (use BED8 inputs)
  mkdir -p "${HAL_OUT}/"{lift,tmp}
  cat "${HAL_OUT}"/lift/*.bed8 > "${HAL_OUT}/all_on_ref.bed8"
  local loci_bed="${HAL_OUT}/pangene_loci.ref.bed6"
  local members="${HAL_OUT}/pangene_members.tsv"
  cluster_lifted_to_loci "${HAL_OUT}/all_on_ref.bed8" "$loci_bed" "$members"

  # --- PAV + gene-ID map
  build_pav_and_geneid_maps "$members"

  # --- Optional: anchor loci to the graph
  #anchor_loci_to_graph "$loci_bed" "${HAL_OUT}/pangene_loci.graph.bed6"
}

make_gene_to_locus() {
  local members="${HAL_OUT}/pangene_members.tsv"
  local out="${HAL_OUT}/gene_to_locus.tsv"
  # members cols: locus  region  gene_id  ASM  CTX
  awk 'BEGIN{FS=OFS="\t"} {print $3, $1, $4, $5, $2}' "$members" \
  | awk 'BEGIN{OFS="\t"} NR==1{print "gene_id","pangene_id","asm","ctx","region"; next} {print}' > "$out"
  echo "gene_to_locus: $out"
}

# Catalog PAV columns and their accession/haplotype parsing
make_pav_columns_index() {
  local pav="${HAL_OUT}/PAV_matrix.tsv"
  local out="${HAL_OUT}/pav_columns.tsv"

  awk -F'\t' -v OFS='\t' '
    BEGIN { print "label","accession","hap" }
    NR==1 {
      for (i=2; i<=NF; i++) {
        lab=$i; acc=lab; hap="0"
        if (lab ~ /[.#][0-9]+$/) {              # Acc.1 or Acc#2
          acc = lab; sub(/[.#][0-9]+$/, "", acc)
          hap = gensub(/^.*[.#]([0-9]+)$/, "\\1", "g", lab)
        } else if (lab ~ /[.#](ref)$/i) {       # Acc.ref or Acc#ref
          acc = lab; sub(/[.#](ref)$/i, "", acc)
          hap = "0"
        } else if (tolower(lab) == "ref") {     # lone "ref"
          acc = "ref"; hap = "0"
        } else {
          acc = lab; hap = "0"
        }
        print lab, acc, hap
      }
      exit
    }
  ' "$pav" > "$out"

  # Sanity check
  local n_header_cols n_rows
  n_header_cols=$(head -n1 "$pav" | awk -F'\t' '{print NF-1}')
  n_rows=$(($(wc -l < "$out") - 1))
  if [[ "$n_header_cols" -ne "$n_rows" ]]; then
    echo "WARN: pav_columns.tsv rows ($n_rows) != PAV header columns ($n_header_cols)" >&2
  fi
  echo "pav_columns: $out"
}

#################################
# Variant extraction (WAVE VCF) + pangene intersections
#################################

mkdir -p "${VAR_OUT}"/{vcf,bed,intersect,summary}

# Helpers
build_chrom_sizes_from_vcf_header() {
  local vcf="$1" out="$2"
  [[ -s "$vcf" ]] || { echo "ERROR: VCF missing: $vcf" >&2; return 1; }
  zcat -f "$vcf" | awk '
    BEGIN{OFS="\t"}
    /^##contig=<ID=/{
      id=""; len="";
      if (match($0, /ID=([^,>]+)/, m)) id=m[1];
      if (match($0, /length=([0-9]+)/, n)) len=n[1];
      if (id!="" && len!="") print id, len
    }' > "$out"
  echo "chrom.sizes → $out"
}

# Prepare & split VCF
wave_prepare_vcfs() {
  local in_vcf="$1"
  [[ -s "$in_vcf" ]] || { echo "ERROR: VCF missing: $in_vcf" >&2; return 1; }

  bcftools index -t -f "$in_vcf" 2>/dev/null || true
  bcftools sort -Oz -o "${VAR_OUT}/vcf/wave.all.sorted.vcf.gz" "$in_vcf"
  bcftools index -t -f "${VAR_OUT}/vcf/wave.all.sorted.vcf.gz"

  # SVs >= 50 bp by |LEN| OR any record with INV flag
  bcftools view -i 'INFO/INV=1 || INFO/LEN<=-50 || INFO/LEN>=50' \
    -Oz -o "${VAR_OUT}/vcf/wave.sv50.vcf.gz" "${VAR_OUT}/vcf/wave.all.sorted.vcf.gz"
  bcftools index -t -f "${VAR_OUT}/vcf/wave.sv50.vcf.gz"

  # SNP+MNP only
  bcftools view -i 'INFO/TYPE="snp" || INFO/TYPE="mnp"' \
    -Oz -o "${VAR_OUT}/vcf/wave.snp_mnp.vcf.gz" "${VAR_OUT}/vcf/wave.all.sorted.vcf.gz"
  bcftools index -t -f "${VAR_OUT}/vcf/wave.snp_mnp.vcf.gz"

  echo "VCFs:"
  echo "  ALL : ${VAR_OUT}/vcf/wave.all.sorted.vcf.gz"
  echo "  SV50: ${VAR_OUT}/vcf/wave.sv50.vcf.gz"
  echo "  SNV : ${VAR_OUT}/vcf/wave.snp_mnp.vcf.gz"
}

# VCF → BED8 (one line per ALT)  |  BED8: chrom start end name score strand type len
vcf_to_bed8() {
  local invcf="$1" outbed="$2"
  [[ -s "$invcf" ]] || { echo "WARN: empty VCF: $invcf"; : > "$outbed"; return 0; }

  bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/TYPE\t%INFO/LEN\t%INFO/INV\n' "$invcf" \
  | awk -v OFS="\t" '
      function ABS(x){return x<0?-x:x}
      {
        chrom=$1; pos=$2; id=$3; ref=$4;
        n_alt=split($5,ALT,",");
        n_typ=split($6,TY,",");
        n_len=split($7,LN,",");
        inv = ($8=="1" || tolower($8)=="true") ? 1 : 0;

        if(id=="." || id=="") id=chrom ":" pos ":" ref ">" ALT[1];

        maxn = n_alt; if(n_typ>maxn) maxn=n_typ; if(n_len>maxn) maxn=n_len;
        for(i=1;i<=maxn;i++){
          t = (i<=n_typ && TY[i]!="") ? toupper(TY[i]) : "OTHER";
          L = (i<=n_len && LN[i]!="") ? ABS(LN[i]+0) : 0;

          if(inv==1) t="INV";                       # INV flag overrides

          start=pos-1; end=start+1; score=L; strand=".";
          if(t=="SNP" || t=="INS" || t=="OTHER"){ end=start+1; }
          else if(t=="MNP"){ end=start+length(ref); }
          else if(t=="DEL" || t=="COMPLEX"){ if(L<=0)L=length(ref); end=start+(L>0?L:1); score=L; }
          else if(t=="INV"){ end=start+(L>0?L:1); }

          if(end<start) end=start;                  # safety
          name = id (n_alt>1 ? (":" i) : "");
          print chrom, start, end, name, score, strand, t, L
        }
      }' \
  | bedtools sort -i - > "$outbed"
}

# Intersect pangene loci (BED6) with BED8 variants and summarize
intersect_and_summarize() {
  local loci_bed="$1" var_bed="$2" tag="$3"
  local raw_out="${VAR_OUT}/intersect/pangene_${tag}.intersect.tsv"
  local sum_out="${VAR_OUT}/summary/pangene_${tag}.counts.tsv"

  [[ -s "$loci_bed" ]] || { echo "ERROR: loci BED missing: $loci_bed" >&2; return 1; }
  [[ -s "$var_bed"  ]] || { echo "WARN: variant BED empty: $var_bed" >&2; : > "$raw_out"; : > "$sum_out"; return 0; }

  bedtools intersect -wa -wb -a "$loci_bed" -b "$var_bed" \
  | awk -v OFS="\t" '
      { # A(1-6) + B(7-14)
        pangene=$4;
        print pangene, $1, $2, $3, $7, $8, $9, $10, $13, $14
      }' \
  | awk 'BEGIN{OFS="\t"; print "pangene_id","locus_chr","locus_start","locus_end","var_chr","var_start","var_end","var_id","var_type","var_len"}1' \
  > "$raw_out"

  # Per pangene counts
  awk -v OFS="\t" 'NR>1{
        p=$1; t=tolower($9); if(t=="") t="other";
        seen[p]=1; C[p,t]++; T[p]++;
      }
      END{
        print "pangene_id","snp","mnp","ins","del","inv","complex","other","total";
        for(p in seen){
          s=C[p,"snp"]+0; m=C[p,"mnp"]+0; i=C[p,"ins"]+0; d=C[p,"del"]+0;
          v=C[p,"inv"]+0; c=C[p,"complex"]+0;
          o=T[p]-s-m-i-d-v-c;
          printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", p,s,m,i,d,v,c,o,T[p];
        }
      }' "$raw_out" | sort -k1,1 > "$sum_out"

  echo "Intersect: $raw_out"
  echo "Summary  : $sum_out"
}

# Route for variant prep + intersections
wave_route_variants() {
  local loci_bed="${HAL_OUT}/pangene_loci.ref.bed6"

  conda deactivate || true
  conda activate "${CONDA_ENV_VARIANT}"

  build_chrom_sizes_from_vcf_header "$WAVE_VCF" "${VAR_OUT}/chrom.sizes"
  wave_prepare_vcfs "$WAVE_VCF"

  # BED8s
  vcf_to_bed8 "${VAR_OUT}/vcf/wave.all.sorted.vcf.gz" "${VAR_OUT}/bed/wave.all.bed8"
  vcf_to_bed8 "${VAR_OUT}/vcf/wave.sv50.vcf.gz"      "${VAR_OUT}/bed/wave.sv50.bed8"
  vcf_to_bed8 "${VAR_OUT}/vcf/wave.snp_mnp.vcf.gz"   "${VAR_OUT}/bed/wave.snp_mnp.bed8"

  # Intersects + summaries
  intersect_and_summarize "$loci_bed" "${VAR_OUT}/bed/wave.sv50.bed8"     "sv50"
  intersect_and_summarize "$loci_bed" "${VAR_OUT}/bed/wave.all.bed8"      "all"
  intersect_and_summarize "$loci_bed" "${VAR_OUT}/bed/wave.snp_mnp.bed8"  "snp_mnp"
}

# Prepare output dirs
mkdir -p "${HAL_OUT}"/{lift,tmp}
mkdir -p "${VAR_OUT}"/{vcf,bed,intersect,summary}


hal_route_build_pangenes
make_gene_to_locus
make_pav_columns_index
wave_route_variants
