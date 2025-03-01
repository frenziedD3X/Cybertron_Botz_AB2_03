import sqlite3
import pandas as pd
import random
from datetime import datetime

# Connect to SQLite database (or create if not exists)
conn = sqlite3.connect('patient_data.db')
c = conn.cursor()

# Create tables
c.execute('''
CREATE TABLE IF NOT EXISTS patients (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    patient_id TEXT UNIQUE,
    name TEXT NOT NULL,
    birthdate TEXT NOT NULL,
    age INTEGER,
    gender TEXT
)
''')

c.execute('''
CREATE TABLE IF NOT EXISTS allergies (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    patient_id TEXT,
    allergy TEXT,
    FOREIGN KEY (patient_id) REFERENCES patients(patient_id)
)
''')

c.execute('''
CREATE TABLE IF NOT EXISTS conditions (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    patient_id TEXT,
    condition TEXT,
    FOREIGN KEY (patient_id) REFERENCES patients(patient_id)
)
''')

c.execute('''
CREATE TABLE IF NOT EXISTS medications (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    patient_id TEXT,
    medications TEXT,
    FOREIGN KEY (patient_id) REFERENCES patients(patient_id)
)
''')

# Function to calculate age
def calculate_age(birthdate):
    try:
        birth_date = datetime.strptime(birthdate, '%Y-%m-%d')
    except ValueError:
        birth_date = datetime.strptime(birthdate, '%d-%m-%Y')
    today = datetime.today()
    age = today.year - birth_date.year - ((today.month, today.day) < (birth_date.month, birth_date.day))
    return age

# Function to generate random names
def generate_random_name():
    first_names = ['John', 'Jane', 'Alex', 'Chris', 'Pat', 'Taylor', 'Jordan', 'Morgan', 'Casey', 'Dana',
               'Sam', 'Jamie', 'Charlie', 'Cameron', 'Riley', 'Logan', 'Avery', 'Blake', 'Quinn', 'Skylar',
               'Finley', 'Peyton', 'Dakota', 'Emerson', 'Hayden', 'Hunter', 'Reese', 'Sawyer', 'Rowan', 'Sydney']

    last_names = ['Smith', 'Johnson', 'Brown', 'Taylor', 'Anderson', 'Thomas', 'Jackson', 'White', 'Harris', 'Martin',
              'Clark', 'Lewis', 'Walker', 'Allen', 'Young', 'King', 'Wright', 'Scott', 'Green', 'Adams',
              'Baker', 'Gonzalez', 'Nelson', 'Carter', 'Mitchell', 'Perez', 'Roberts', 'Campbell', 'Edwards', 'Collins']

    return f"{random.choice(first_names)} {random.choice(last_names)}"

# Function to insert patient data
def insert_patient(patient_id, birthdate, gender):
    if patient_id.strip():
        name = generate_random_name()
        age = calculate_age(birthdate)
        c.execute('''
        INSERT OR IGNORE INTO patients (patient_id, name, birthdate, age, gender)
        VALUES (?, ?, ?, ?, ?)
        ''', (patient_id.strip(), name, birthdate, age, gender))
        conn.commit()
    else:
        print(f"Invalid patient_id: {patient_id}")

# Function to insert allergy data
def insert_allergy(patient_id, allergy):
    c.execute('''
    INSERT INTO allergies (patient_id, allergy)
    VALUES (?, ?)
    ''', (patient_id.strip(), allergy))
    conn.commit()

# Function to insert condition data
def insert_condition(patient_id, condition):
    c.execute('''
    INSERT INTO conditions (patient_id, condition)
    VALUES (?, ?)
    ''', (patient_id.strip(), condition))
    conn.commit()

# Function to insert medication data
def insert_medication(patient_id, medications):
    c.execute('''
    INSERT INTO medications (patient_id, medications)
    VALUES (?, ?)
    ''', (patient_id.strip(), medications))
    conn.commit()

# Function to load data from CSVs and insert data
def load_csvs_and_insert(patients_file, allergies_file, conditions_file, medications_file):
    patients_df = pd.read_csv(patients_file, usecols=['Id', 'BIRTHDATE', 'GENDER'])
    for index, row in patients_df.iterrows():
        insert_patient(str(row['Id']).strip(), row['BIRTHDATE'], row['GENDER'])

    allergies_df = pd.read_csv(allergies_file, usecols=['Id', 'DESCRIPTION'])
    for index, row in allergies_df.iterrows():
        insert_allergy(str(row['Id']).strip(), row['DESCRIPTION'])

    conditions_df = pd.read_csv(conditions_file, usecols=['Id', 'DESCRIPTION'])
    for index, row in conditions_df.iterrows():
        insert_condition(str(row['Id']).strip(), row['DESCRIPTION'])

    medications_df = pd.read_csv(medications_file, usecols=['Id', 'DESCRIPTION'])
    for index, row in medications_df.iterrows():
        insert_medication(str(row['Id']).strip(), row['DESCRIPTION'])

# Example CSV Load
# load_csvs_and_insert('patients.csv', 'allergies.csv', 'conditions.csv', 'medications.csv')

# Close connection
# conn.close()


# Example CSV Load
load_csvs_and_insert('Final_CSV/patients.csv', 'Final_CSV/allergies.csv', 'Final_CSV/conditions.csv', 'Final_CSV/medications.csv')

# Close connection
conn.close()
