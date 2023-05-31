'''
Created on May 28, 2023

@author: wku
'''
import requests

import re



def get_mirbase_version():

    try:

        url = "https://www.mirbase.org/ftp/CURRENT/README"

        response = requests.get(url)

        response.raise_for_status()  # Raise an exception for HTTP errors

        readme_text = response.text

        # Search for the version information in the README file

        version_match = re.search(r"Release (\d+\.\d+)", readme_text)

        if version_match:

            miRBase_version = version_match.group(1)

            return miRBase_version

    except requests.RequestException as e:

        print("An error occurred while making the request:", e)

    except re.error:

        print("An error occurred while searching for the version information.")

    

    return None
# Call the function to retrieve the miRBase version

mirbase_version = get_mirbase_version()



if mirbase_version:

    print("Current version of miRBase:", mirbase_version)

else:

    print("Unable to retrieve the current version of miRBase.")
