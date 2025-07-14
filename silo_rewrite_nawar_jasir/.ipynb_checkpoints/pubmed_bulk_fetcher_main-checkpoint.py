from Bio import Entrez
import pandas as pd
from silo_utils_TESTING import *

#TESTING COPY
# Rewrite of silo_transplant by Timur Seckin to leverage the Entrez api for better reproducibilty
"""
 Original URL: "https://pubmed.ncbi.nlm.nih.gov/?term=%28Journal+of+the+American+Society+of+Nephrology%5BJournal%5D+OR+American+Journal+of+Kidney+Diseases%5BJournal%5D+OR+Clinical+Journal+of+the+American+Society+of+Nephrology%5BJournal%5D+OR+American+Journal+of+Transplantation%5BJournal%5D+OR+Transplantation%5BJournal%5D+OR+Nephrology+Dialysis+Transplantation%5BJournal%5D%29+AND+%28clinical+research+OR+original+research+OR+basic+research+OR+clinical+investigation+OR+clinical+epidemiology+OR+original+investigation+OR+original+article+OR+original+basic+science+OR+original+clinical+science%29+AND+2023%5BPublication+Date%5D&size=200"
"""

silo_query = "(Journal of the American Society of Nephrology[Journal] OR American Journal of Kidney Diseases[Journal] OR Clinical Journal of the American Society of Nephrology[Journal] OR American Journal of Transplantation[Journal] OR Transplantation[Journal] OR Nephrology Dialysis Transplantation[Journal]) AND (clinical research OR original research OR basic research OR clinical investigation OR clinical epidemiology OR original investigation OR original article OR original basic science OR original clinical science) AND 2023[Publication Date]"
# Kidney International?
ret_num = 622
probe = search(silo_query,ret_num)

scraped_article_list = []

for article_id in probe['IdList']:
    print("Scraping article: ", article_id)
    scraped_article_list.append(article_scraper(retrieve_data(article_id)))

df = pd.DataFrame(scraped_article_list)
df.to_csv("full.csv", index=False)
print("Data extraction complete. Check full.csv for results.")