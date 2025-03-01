import logging
import chromadb
from sentence_transformers import SentenceTransformer
from pubmed_api import PubMedClient
email = "xyzcoc1234@gmail.com"

class VectorDatabase:
    def __init__(self, collection_name="loaded_articles"):
        self.client = chromadb.PersistentClient(path="./chroma_db")
        self.collection = self.client.get_or_create_collection(name=collection_name)
        self.model = SentenceTransformer('sentence-transformers/all-MiniLM-L6-v2')
        self.logger = logging.getLogger(__name__)

    def chunk_text(self, text, chunk_size=512):
        """Splits text into smaller chunks for better vectorization"""
        chunks = [text[i:i + chunk_size] for i in range(0, len(text), chunk_size)]
        self.logger.info(f"Text split into {len(chunks)} chunks")
        return chunks

    def embed_and_store(self, articles):
        """Converts articles into vector embeddings and stores them in the vector database"""
        for article in articles:
            content = f"{article['title']}. {article['abstract']}"
            chunks = self.chunk_text(content)
            for i, chunk in enumerate(chunks):
                vector = self.model.encode(chunk).tolist()
                metadata = {
                    "title": article['title'],
                    "chunk_index": i,
                    "full_text_link": article['full_text_link']
                }
                self.collection.add(documents=[chunk], metadatas=[metadata], embeddings=[vector], ids=[f"{article['title']}-{i}"])
                self.logger.info(f"Stored chunk {i} for article: {article['title']}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    vector_db = VectorDatabase()

    sample_articles = [
        {
            "title": "Psoriasis Study 1",
            "abstract": "This study investigates the effects of psoriasis on skin inflammation...",
            "full_text_link": "https://doi.org/10.1234/example1"
        },
        {
            "title": "Psoriasis Study 2",
            "abstract": "A detailed look into genetic factors influencing psoriasis...",
            "full_text_link": "https://doi.org/10.1234/example2"
        }
    ]

    vector_db.embed_and_store(sample_articles)
