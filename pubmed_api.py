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
                full_text_link = self._get_doi(article)

                articles.append({
                    "title": title,
                    "abstract": abstract_text.strip(),
                    "full_text_link": full_text_link
                })
            self.logger.info("Article fetch completed")
        except Exception as e:
            self.logger.error(f"Fetch failed: {e}")
        return articles

    def _get_doi(self, article):
        """Extract DOI link if available"""
        for article_id in article.findall(".//ArticleId"):
            if article_id.attrib.get("IdType") == "doi":
                return f"https://doi.org/{article_id.text}"
        return "Not Available"


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


