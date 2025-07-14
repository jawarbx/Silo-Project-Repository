from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service as ChromeSerivce
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import pandas as pd
import time

# Path to your Firefox profile directory
profile_path ='/Users/jasirnawar/Library/Application Support/Google/Chrome/Default'

# Configure Firefox options
options = Options()
options.add_argument(
            r"user-data-dir=/Users/jasirnawar/Library/Application Support/Google/Chrome/Default"
        )
options.add_argument("--profile-directory=Default")
options.add_argument("--start-maximized")
# Initialize the Firefox driver with the profile

df = pd.read_csv('Article_scraped.csv').fillna("")
blanks = df.loc[df['ARTICLES'] == ""]

for doi in blanks['DOIS']: 
    doi_url = doi  # Replace with your DOI URL
    driver = webdriver.Firefox(service=ChromeSerivce(ChromeDriverManager().install()), options=options)
    try:
        # Navigate to the URL
        driver.get(doi_url)
        
        # Wait for the final URL to load
        WebDriverWait(driver, 10).until(lambda d: d.execute_script('return document.readyState') == 'complete')
        final_url = driver.current_url
        
        # Navigate to the final URL (if needed)
        driver.get(final_url)
        
        # Ensure the page is fully loaded
        WebDriverWait(driver, 10).until(lambda d: d.execute_script('return document.readyState') == 'complete')
        time.sleep(3)

        if "sciencedirect" in final_url:
            element = driver.find_element(By.ID, "body")
            script = """
            arguments[0].removeChild(arguments[0].querySelector('[id^="ack"]'))
            list = arguments[0].querySelectorAll('[id^="p0"]')
            text = []
            for (let i=0; i<list.length; i++) {
                text.push(list[i].textContent)
            }
            return text.toString()
            """
            df['ARTICLES'].loc[df['DOIS'] == doi] = driver.execute_script(script, element)
        elif "journals.lww.com" in final_url:
            element = driver.find_element(By.ID, "ArticleBody")
            script = """
            list = arguments[0].querySelectorAll('[id^="JCL-P"]')
            text = []
            for (let i=0; i<list.length; i++) {
                text.push(list[i].textContent)
            }
            return text.toString()
            """
            df['ARTICLES'].loc[df['DOIS'] == doi] = driver.execute_script(script, element)
    except:
       print("Error for: ", doi_url)
       pass
    finally:
        # Close the browser
        driver.quit()
df.to_csv("Article_scraped_Elsevier.csv")
