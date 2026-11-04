import pandas as pd
import requests
import time

# --------------------------------------------------
# FILE SETTINGS
# --------------------------------------------------

input_file = "/storage/atedder/DESeq2_24h_rhythmic_genes_with_phase_amplitude.csv"
output_file = "/storage/atedder/DESeq2_24h_rhythmic_genes_with_phase_amplitude_filled.csv"

print("\nLoading CSV file...")

df = pd.read_csv(input_file)

print(f"Successfully loaded {len(df)} rows.")
print(f"Columns found: {list(df.columns)}")


# --------------------------------------------------
# IDENTIFY MISSING GENE NAMES
# --------------------------------------------------

mask = (
    df["gene_name"].isna() |
    (df["gene_name"].astype(str).str.strip().str.upper() == "NA") |
    (df["gene_name"].astype(str).str.strip() == "")
)

missing_rows = df[mask].index.tolist()

print("\nIdentifying missing gene names...")
print(f"Total rows: {len(df)}")
print(f"Rows requiring annotation: {len(missing_rows)}")
print(f"Rows already annotated: {len(df) - len(missing_rows)}")


# --------------------------------------------------
# ENSEMBL API SETTINGS
# --------------------------------------------------

url = "https://rest.ensembl.org/lookup/id/"

headers = {
    "Content-Type": "application/json"
}

# Counters
successful = 0
not_found = 0
errors = 0

# Avoid querying the same ID multiple times
gene_cache = {}


# --------------------------------------------------
# RETRIEVE MISSING GENE NAMES
# --------------------------------------------------

print("\nStarting Ensembl queries...\n")

for count, idx in enumerate(missing_rows, start=1):

    gene_id = df.loc[idx, "gene_id"]

    print(f"[{count}/{len(missing_rows)}] Querying: {gene_id}")

    # Check whether this ID has already been queried
    if gene_id in gene_cache:

        gene_name = gene_cache[gene_id]

        print("  Using previously retrieved annotation.")

    else:

        try:

            response = requests.get(
                url + gene_id,
                headers=headers,
                timeout=15
            )

            if response.status_code == 200:

                data = response.json()

                gene_name = data.get("display_name")

                gene_cache[gene_id] = gene_name

            else:

                gene_name = None
                gene_cache[gene_id] = None

        except requests.exceptions.RequestException as e:

            print(f"  CONNECTION ERROR: {e}")
            errors += 1
            time.sleep(1)
            continue

    # Update dataframe if gene name was retrieved
    if gene_name:

        df.loc[idx, "gene_name"] = gene_name

        successful += 1

        print(f"  SUCCESS: {gene_id} -> {gene_name}")

    else:

        not_found += 1

        print("  NOT FOUND: No gene name returned.")

    # Small delay between API requests
    time.sleep(0.1)


# --------------------------------------------------
# SAVE UPDATED FILE
# --------------------------------------------------

print("\nSaving updated CSV file...")

df.to_csv(output_file, index=False)


# --------------------------------------------------
# FINAL SUMMARY
# --------------------------------------------------

print("\n" + "=" * 50)
print("ANNOTATION COMPLETE")
print("=" * 50)

print(f"Total rows in original file: {len(df)}")
print(f"Rows requiring annotation: {len(missing_rows)}")
print(f"Successfully annotated: {successful}")
print(f"IDs not found: {not_found}")
print(f"Connection errors: {errors}")

print(f"\nUpdated file saved as: {output_file}")

print("\nFinished!")
