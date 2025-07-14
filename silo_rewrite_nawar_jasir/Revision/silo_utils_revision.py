from Bio import Entrez
import traceback

#TESTING COPY
email = 'jn2691@nyu.edu'
"""
    dict of journal impact factors
"""
Entrez.email = email
JOURNAL_IMPACT_FACTORS = {
    "Journal of the American Society of Nephrology : JASN": 13.6,
    "American journal of kidney diseases : the official journal of the National Kidney Foundation": 13.2,
    "Clinical journal of the American Society of Nephrology : CJASN": 9.8,
    "American journal of transplantation : official journal of the American Society of Transplantation and the American Society of Transplant Surgeons": 8.8,
    "Transplantation": 6.2,
    "Nephrology, dialysis, transplantation : official publication of the European Dialysis and Transplant Association - European Renal Association": 6.1,
}

def search(query, ret_num):
  """
    Search for a string in pubmed 
    Input:
        str query: term to be searched on pubmed. Can be "Pubmed style" conditional search
        int ret_num: Amount of results to be retreived
    Output:
        dict results: results of search in xml format
"""
  Entrez.email = email
  handle = Entrez.esearch(
      db=DB,
      retmax=ret_num,
      term=query,
      sort='relevance',
  )
  results = Entrez.read(handle)
  handle.close()
  return results

def retrieve_data(input_id, data_type):
  """
  Retrieve article data based on article uid \\
  Input:
    str article_id: article_id to be scraped
    retmode: specify the return mode
    rettype: specify return type
  Output:
    XML format of desired data
  """

    match data_type:
        case "article":
            handle = Entrez.efetch(
                db="pubmed",
                id=input_id,
            )
        case "summary":
            handle = Entrez.esummary(
                db="pubmed",
                id=input_id,
            )
        case "full_text":
            handle = Entrez.esummary(
                db="pmc",
                id=input_id,
            )
        case _:
            except:
                print(f"Data Type {data_type} not available")
                return 
    results = handle.read()
    handle.close()
    return results

   

def article_scraper(bulk_article_data):
    title = abstract = journal_name = grants_and_funding = conflict_of_interest = impact_factor = mesh_terms = pmid = pmcid = pub_type_list = doi = ""
    authors = affiliations = grant_list = mesh_list_descriptors = []

    try:
       title = bulk_article_data['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle']
    except:
        print("Could not retrieve article title:", traceback.print_exc())
        pass

    try:
       for author in bulk_article_data['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']:
          authors.append(author['ForeName'] + " " + author['LastName'])
    except:
       print("Could not retrieve authors:", traceback.print_exc())
       pass

    try:
       abstract = bulk_article_data['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
    except:
       print("Could not retrieve abstract:", traceback.print_exc())
       pass

    try:
       unique_affiliations = set()
       for author in bulk_article_data['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']:
           for aff_info in author['AffiliationInfo']:
               aff = aff['Affiliation']
               if aff not in unique_affiliations:
                   unique_affiliations.add(aff)
    except:
        print("Could not retrieve affiliations:", traceback.print_exc())
        pass

    try:
       grant_list = bulk_article_data['PubmedArticle'][0]['MedlineCitation']['Article']['GrantList']
    except:
       print("Could not retrieve grant list: ", traceback.print_exc())
       pass
    else:
       for grant in grant_list:
        id = "{NO ID}"
        acrnym = "{NO ACRNYM}"
        agency = "{NO AGENCY}"
        country = "{NO COUNTRY}"
        
        try:
           id = grant['GrantID']
        except:
           pass
        
        try:
           acrnym = grant['Acronym']
        except:
           pass
        
        try:
           agency = grant['Agency']
        except:
           pass
        
        try:
           country = grant['Country']
        except:
           pass        
        
        
        grants_and_funding += '||'.join((id, acrnym, agency, country)) + "[]"
        pass
       
    try: 
       conflict_of_interest = bulk_article_data['PubmedArticle'][0]['MedlineCitation']['CoiStatement']
    except:
       print("Could not retrieve coi statement: ", traceback.print_exc())
       pass

    try:
       journal_name = bulk_article_data['PubmedArticle'][0]['MedlineCitation']['Article']['Journal']['Title']
    except:
       print("Could not retrieve journal name: ", traceback.print_exc())
       pass
    else:
       if journal_name in JOURNAL_IMPACT_FACTORS:
             impact_factor = str(JOURNAL_IMPACT_FACTORS[journal_name])
       pass
    
    try:
       mesh_terms = bulk_article_data['PubmedArticle'][0]['MedlineCitation']['MeshHeadingList']
    except:
       print("Could not retrieve mesh_terms", traceback.print_exc())
       pass
    else:
       for term in mesh_terms:
          mesh_list_descriptor_ids.append(term['DescriptorName'].attributes['UI'])
          mesh_list_descriptors.append(term['DescriptorName'])
       pass
    
    try:
       pub_type_list = bulk_article_data['PubmedArticle'][0]['MedlineCitation']['Article']['PublicationTypeList']
    except:
       print("Could not get pub type", traceback.print_exc())
       pass
    try:
       pmid = bulk_article_data['PubmedArticle'][0]['PubmedData']['ArticleIdList'][0]
    except:
       print("Could not get pmid", traceback.print_exc())
       pass
    try:
        for str_elem in bulk_article_data['PubmedArticle'][0]['PubmedData']['ArticleIdList']:
            if str_elem.attributes['IdType'] == 'pmc':
                pmcid = int(str_elem[3:])
            if str_elem.attributes['IdType'] == 'doi':
                doi = str(str_elem)
    except:
        print("No other ids")
        pass
        
    return {

        "Title": title,

        "PMID": pmid,

        "PMCID" : pmcid,

        "DOI" : doi,

        "Authors": ','.join(authors),

        "Affiliations": '////'.join(affiliations),

        "Abstract": abstract,

        "Journal": journal_name,

        "Grants_and_Funding": grants_and_funding,

        "Conflict_of_Interest": conflict_of_interest,
        
        "Impact_Factor": impact_factor,

        "Mesh Terms" : ','.join(mesh_list_descriptors),

        "MeSH Term Descriptor IDs": ','.join(mesh_list_descriptor_ids),
        
        "Total Affiliations": len(affiliations),
        
        "Total Authors": len(authors),      

        "Pub Type" : ",".join(pub_type_list),

        
    }
def retrieve_full_text(pmc_id,retmode="text"):
  """
  Retrieve article data based on article uid 
  Input:
    str article_id: article_id to be scraped
  Output:
    dict results: results of search in dict format
  """
  Entrez.email = email
  handle = Entrez.efetch(
      db="pmc",
      id=pmc_id,
      retmode=retmode
  )
    
  results = handle.read()
  handle.close()
  return results

def scrape_full_text(bulk_pmc_article):
    """
    Scrape the article for extra pubmed variables
    Input: pmc article
    Output: dict hopefully
    """
    return