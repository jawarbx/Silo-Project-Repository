# How to use this

(Optional) Create a virtual environment OR use (Ana)conda

1. Install required packages
    `pip install numpy pandas biopython`

2. In silo_utils.py, be sure to change the `email` to your email such that Entrez can email you if the job may be too much for their servers.

4. In silo_rewrite.py, ret_num refers to the number of articles you wish to fetch. Adjust as accordingly.
   I've currently set the default to 40, but the max for this given query is 622.

5. Run `python silo_rewrite.py`

6. The desired information will be in the file labeled "goal.csv"
