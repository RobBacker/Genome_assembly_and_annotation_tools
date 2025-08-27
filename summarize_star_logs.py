import os
import glob
import pandas as pd
import sys

# Get directory from command-line argument
if len(sys.argv) < 2:
    print("Usage: python summarize_star_logs.py /path/to/logs/")
    sys.exit(1)

log_dir = sys.argv[1].rstrip("/")  # strip trailing slash if present
log_files = glob.glob(os.path.join(log_dir, "*.STAR.Log.final.out"))

if not log_files:
    print(f"No STAR logs found in: {log_dir}")
    sys.exit(1)

def parse_star_log(filepath):
    metrics = {}
    with open(filepath) as f:
        for line in f:
            if "|" not in line:
                continue
            key, val = line.strip().split("|")
            metrics[key.strip()] = val.strip()

    return {
        "sample": os.path.basename(filepath).replace(".STAR.Log.final.out", ""),
        "input_reads": int(metrics.get("Number of input reads", 0)),
        "unique_map_pct": float(metrics.get("Uniquely mapped reads %", "0%").strip("%")),
        "multi_map_pct": float(metrics.get("% of reads mapped to multiple loci", "0%").strip("%")),
        "too_short_pct": float(metrics.get("% of reads unmapped: too short", "0%").strip("%")),
        "total_splices": int(metrics.get("Number of splices: Total", 0)),
    }

results = [parse_star_log(log) for log in log_files]
df = pd.DataFrame(results)

# Flag low-quality libraries
df["flag_low_quality"] = (
    (df["unique_map_pct"] < 50) |
    (df["input_reads"] < 5_000_000) |
    (df["too_short_pct"] > 40) |
    (df["total_splices"] < 1_000_000)
)

# Output
# Save CSV to the same directory as the logs
output_csv = os.path.join(log_dir, "star_log_summary.csv")
df.to_csv(output_csv, index=False)
print(f"Summary saved to {output_csv}")
print(df[df["flag_low_quality"]])