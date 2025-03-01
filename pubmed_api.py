from Bio import Entrez
from xml.etree import ElementTree
import logging

class PubMedClient:
    def __init__(self, email, api_key=None):
        self.email = email
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.logger = logging.getLogger(__name__)

    def search(self, query, retmax=10, retstart=0):
        """Search PubMed articles by query with pagination support"""
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax, retstart=retstart)
            record = Entrez.read(handle)
            handle.close()
            self.logger.info(f"Search completed for query: {query}")
            return record.get("IdList", [])
        except Exception as e:
            self.logger.error(f"Search failed: {e}")
            return []

    def fetch_articles(self, pmid_list):
        """Fetch article details by PubMed IDs"""
        articles = []
        try:
            handle = Entrez.efetch(db="pubmed", id=",".join(pmid_list), retmode="xml")
            tree = ElementTree.parse(handle)
            root = tree.getroot()
            self.logger.info(f"Fetching {len(pmid_list)} articles")

            for article in root.findall(".//PubmedArticle"):
                title = article.find(".//ArticleTitle").text if article.find(".//ArticleTitle") is not None else "No Title"
                abstract_text = " ".join([ab.text for ab in article.findall(".//AbstractText") if ab.text]) or "No Abstract"

                articles.append({
                    "title": title,
                    "abstract": abstract_text.strip()
                })
            self.logger.info("Article fetch completed")
        except Exception as e:
            self.logger.error(f"Fetch failed: {e}")
        return articles

    def fetch_and_save_multiple_queries(self, queries, filename, retmax=10):
        """Fetch articles for multiple queries and save results to a text file"""
        try:
            with open(filename, 'w', encoding='utf-8') as file:
                for query in queries:
                    self.logger.info(f"Processing query: {query}")
                    pmid_list = self.search(query, retmax=retmax)
                    articles = self.fetch_articles(pmid_list)
                    
                    file.write(f"Query: {query}\n")
                    for article in articles:
                        file.write(f"Title: {article['title']}\n")
                        file.write(f"Abstract: {article['abstract']}\n")
                        file.write("\n")
                    file.write("="*50 + "\n")
            self.logger.info(f"Results saved to {filename}")
        except Exception as e:
            self.logger.error(f"Failed to save articles: {e}")


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


if __name__ == "__main__":
    # Initialize the PubMedClient
    client = PubMedClient(email="your_email@example.com")

    # List of queries
    queries = [
        "hypertension",
        "regular medication for diabetes",
        "allergy to penicillin",
        "chest pain"
    ]

    # Output filename
    filename = "abstract.txt"

    # Fetch and Save Articles
    client.fetch_and_save_multiple_queries(queries, filename, retmax=5)

    print(f"Articles saved to {filename}")
