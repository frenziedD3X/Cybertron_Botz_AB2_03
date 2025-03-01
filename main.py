import os
import numpy as np
from transformers import AutoTokenizer, AutoModel
from groq import Groq
import chromadb
from chromadb.utils import embedding_functions
import torch
import asyncio
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Load BioBERT model for embeddings
MODEL_NAME = "dmis-lab/biobert-base-cased-v1.1"
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
model = AutoModel.from_pretrained(MODEL_NAME)

# Function to create embeddings using Hugging Face model with mean pooling
def get_embedding(text):
    inputs = tokenizer(text, return_tensors="pt", padding=True, truncation=True, max_length=512)
    with torch.no_grad():
        outputs = model(**inputs).last_hidden_state
    embeddings = outputs.mean(dim=1).squeeze().numpy()
    return embeddings.tolist()

# Load abstracts from text file
def load_abstracts(file_path):
    with open(file_path, 'r') as f:
        abstracts = f.readlines()
    return [abstract.strip() for abstract in abstracts if abstract.strip()]

# Initialize ChromaDB
chroma_client = chromadb.PersistentClient(path="./chroma_db")
collection = chroma_client.get_or_create_collection(name="medical_abstracts")

# Store embeddings in ChromaDB
def store_embeddings(abstracts):
    ids = []
    embeddings = []
    for i, abstract in enumerate(abstracts):
        ids.append(str(i))
        embeddings.append(get_embedding(abstract))
    collection.add(ids=ids, embeddings=embeddings, documents=abstracts)

# Retrieve top K abstracts
def retrieve(text, k=5):
    query_embedding = get_embedding(text)
    results = collection.query(query_embeddings=[query_embedding], n_results=k)
    return results["documents"][0]

# Initialize Groq API
def analyze_with_groq(query, model="deepseek-r1-distill-llama-70b"):
    client = Groq(api_key=os.environ.get("GROQ_API_KEY"))
    messages = [
        {"role": "system", "content": "You are an expert doctor."},
        {"role": "user", "content": query}
    ]
    chat_completion = client.chat.completions.create(
        messages=messages,
        model=model
    )
    return chat_completion.choices[0].message.content

async def generate_streaming_response(patient_info, history, context):
    prompt = f"""
    Patient Info: {patient_info}
    
    Patient History: {history}
    
    Relevant Research Context:
    {context}
    
    Based on the above information, provide a diagnosis and treatment suggestion:
    """
    response = analyze_with_groq(prompt)
    for chunk in response.split():
        yield chunk + " "
        await asyncio.sleep(0.1)

if __name__ == "__main__":
    # Load and embed abstracts
    abstracts = load_abstracts("abstract.txt")
    store_embeddings(abstracts)

    # Example patient data
    patient_data = "Male, 45, diabetes, fever, chest pain, allergy to penicillin"
    patient_history = "Past surgeries, hypertension, regular medication for diabetes"
    top_abstracts = retrieve(patient_data)

    # Generate diagnosis and treatment plan with streaming
    async def main():
        async for chunk in generate_streaming_response(patient_data, patient_history, '\n'.join(top_abstracts)):
            print(chunk, end="")

    asyncio.run(main())
