import deweydatapy as ddp
import pandas
from functools import cache


API_KEY_PATH = "/home/hillandr/creds/2026_03_03_dewey_api_key.txt"


def get_api_key():
    with open(API_KEY_PATH, "r") as f:
        return f.read().rstrip("\n")

# API Key
API_KEY = get_api_key()

# Product Path URL
PRODUCT_PATH="https://api.deweydata.io/api/v1/external/data/prj_qjr4grcp__cdst_zfqg76mugq9dtycr"

if __name__ == "__main__":
    meta = ddp.get_meta(
        apikey=API_KEY,
        product_path=PRODUCT_PATH,
        print_meta=True
    )

    files_df = ddp.get_file_list(
        apikey=API_KEY,
        product_path=PRODUCT_PATH,
        start_date="2019-01-01",
        end_date="2026-04-15"
    )

    ddp.download_files(
        files_df=files_df,
        dest_folder="/biostats_share/hillandr/data/WW_Mobility_2026_04_20/neigh_patterns_plus/",
        skip_exists=True
    )

    print("Done!")


