from flask import Flask, request, jsonify, render_template
import sqlite3
from brain import load_abstracts, store_embeddings, generate_streaming_response, retrieve
from pubmed_api import PubMedClient
app = Flask(__name__)
 
client = PubMedClient(email="your_email@example.com")

# Function to retrieve data from SQLite DB
def get_user_data(name, user_id):
    try:
        conn = sqlite3.connect("patient_data.db")
        cursor = conn.cursor()
        query = "SELECT patient_id FROM patients WHERE name = ? AND id = ?"
        cursor.execute(query, (name, user_id))
        data = cursor.fetchone()
        if data:
            print(data[0])
            query = "SELECT age, gender FROM patients WHERE name = ? AND id = ?"
            cursor.execute(query, (name, user_id))
            age = cursor.fetchone()
            query = "SELECT allergy FROM allergies WHERE patient_id = ?"
            cursor.execute(query, (str(data[0]).strip(),))
            allergys = cursor.fetchall()
            print(allergys)
            query = "SELECT condition FROM conditions WHERE patient_id = ?"
            cursor.execute(query, (str(data[0]).strip(),))
            conditons = cursor.fetchall()
            print(conditons[0])
            query = "SELECT medications FROM medications WHERE patient_id = ?"
            cursor.execute(query, (str(data[0]).strip(),))
            medications = cursor.fetchall()
            print(medications[0][0])
            return age, allergys, conditons, medications
        else:
            print("No data found")
            return None
    except sqlite3.Error as e:
        print(f"Database error: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None
    finally:
        if conn:
            conn.close()



# Function to process the data and prompt
def process_data(data, allergy, conditions, medications, name, prompt):
    if data:
        patient_data = f"{name} age: {data[0]} gender: {data[1]}"
        allergies = ""
        if allergy:
            for allerg in allergy:
                allergies += ", " + allerg[0]
        condition_s = ""
        if conditions:
            for condition in conditions:
                condition_s += ", " + condition[0]
        
        medication_s = ""
        if medications:
            for medication in medications:
                medication_s += ", " + medication[0]
        patient_data += "\n" + f"Allergies: {allergies} \nCondtions: {condition_s} \nMedications: {medication_s} \n Current Query{prompt}"
        print (patient_data)
        return patient_data
    return "No user found with provided credentials."

# New function to generate text output
def generate_text_output(text):
    abstracts = load_abstracts("abstract.txt")
    store_embeddings(abstracts)
    top_abstracts = retrieve(text)
    return generate_streaming_response(text, top_abstracts)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/get_data", methods=["POST"])
def get_data():
    try:
        user_input = request.json
        name = user_input.get("name")
        user_id = user_input.get("id")
        prompt = user_input.get("prompt")
        
        if not name or not user_id or not prompt:
            return jsonify({"error": "Name or ID, and prompt are required"}), 400
        
        user_data, allergys, conditons, medications = get_user_data(name, user_id)
        result = process_data(user_data, allergys, conditons, medications, name, prompt)
        query = result.split('\n')
        client.fetch_and_save_multiple_queries(queries = query, filename='abstract.txt')
        generated_text = generate_text_output(prompt)
        
        return jsonify({"result": result, "generated_text": generated_text})
    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True)