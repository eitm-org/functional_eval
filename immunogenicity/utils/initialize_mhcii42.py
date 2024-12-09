import os
from setuptools import setup, find_packages
import requests
import tarfile

def install_mhc2_req(): 
    mhcii_42_url = "https://www.dropbox.com/scl/fi/ukm72ba6uz6bapiv6ckuc/mhcii_42.tar.gz?rlkey=l3rkcz0hh29tnevhdh76zpw91&st=lblu04w3&dl=1"
    immunogenicity_utils_dir = os.path.join(os.path.dirname(__file__), 'immunogenicity', 'utils')
    os.makedirs(immunogenicity_utils_dir, exist_ok=True)
    tar_fpath = os.path.join(immunogenicity_utils_dir, 'temp_mhcii42.tar.gz')

    print(f"Downloading {mhcii_42_url}...")
    response = requests.get(mhcii_42_url)

    if response.status_code == 200:
        with open(tar_fpath, 'wb') as f:
            f.write(response.content)
        print(f"File downloaded to {tar_fpath}")

        # Unpack the tar file
        print(f"Unpacking {tar_fpath}...")
        with tarfile.open(tar_fpath, 'r') as tar:
            tar.extractall(path=immunogenicity_utils_dir)
        print(f"File unpacked to {immunogenicity_utils_dir}")
        # Optionally, delete the tar file after unpacking (if not needed)
        os.remove(tar_fpath)
        print(f"Tar file removed: {tar_fpath}")
    else:
        print(f"Failed to download file. HTTP status code: {response.status_code}")


    return 