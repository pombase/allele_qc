import pandas
import glob

files = glob.glob("allele_*.tsv") + [
    "description_in_db_doesnt_match_03082023.tsv",
    "manual_cannot_fix_new_03082023.tsv",
    "unknowns_with_correct_descriptions_10082023.tsv",
]
data = pandas.DataFrame()

for file in files:
    print(file)
    df = pandas.read_csv(file, sep="\t", na_filter=False)
    unique_columns = df[["allele_name", "allele_description"]]
    if "change_description_to" in df.columns:
        unique_columns.loc[:, "description_changed"] = df["change_description_to"] != ""
    else:
        unique_columns.loc[:, "description_changed"] = False
    if "change_name_to" in df.columns:
        unique_columns.loc[:, "name_changed"] = df["change_name_to"] != ""
    data = pandas.concat([data, unique_columns])

data = data.drop_duplicates()

print("> alleles fixed (at most)", len(data))
print("> alleles with description changed (at most)", data["description_changed"].sum())
print("> alleles with name changed (at most)", data["name_changed"].sum())
print()
files = glob.glob("protein_*.tsv")
data = pandas.DataFrame()

for file in files:
    df = pandas.read_csv(file, sep="\t", na_filter=False)
    unique_columns = df[["systematic_id", "sequence_position"]]
    data = pandas.concat([data, unique_columns])

data = data.drop_duplicates()

print("> modifications fixed (at most)", len(data))
