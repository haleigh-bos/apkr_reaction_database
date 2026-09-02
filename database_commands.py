# Initializes and updates the SQLite Database containing Reaction data
import sqlite3
from pathlib import Path
import pandas as pd
import re
import os
from datetime import datetime

def build_database(database_name='apkr_reactions.db'):
    """Either builds database from scratch, or opens the database if it exists.

    Args: 
        database_name: default is 'apkr_reactions.db'

    Fails if database could not be opened.
    """
    try:
        with sqlite3.connect(database_name) as conn:
            print(f"Opened SQLite database with version {sqlite3.sqlite_version} successfully.")

    except sqlite3.OperationalError as e:
        print("Failed to open database:", e)

################################################################################################################
def clean_slate(database_name='apkr_reactions.db'):
    """Gets rid of tables that are within all_data. Simultaneously removes these table names from 'added_sheets' table.
    Does not remove the individual sheets that have already been generated.

    Args:
        database_name: default is 'apkr_reactions.db'
    """
    try:
        with sqlite3.connect(database_name) as conn:
            cur = conn.cursor()
            cur.execute("DROP TABLE IF EXISTS all_data")
            cur.execute("DROP TABLE IF EXISTS added_sheets")
            conn.commit()
            print("Dropped tables from 'all_data' and removed entries from 'added_sheets'.")
    except sqlite3.OperationalError as e:
        print("Failed to wipe 'all_data':", e)

################################################################################################################
def update_database(sheet_list, database_name='apkr_reactions.db', data_path='../0.database_excels/'):
    """Takes the list of excel sheets given and appends the 'all_data' table in the database.
    Does not update 'all_data' with tables that have already been added. Names of sheets that have already been
    added have been recorded in 'added_sheets' table in database.

    Args:
        sheet_list: list of sheets from a chosen folder (currently 0.database_excels)
        database_name: default is 'apkr_reactions.db'
        data_path: default is '../0.database_excels/'

    Returns:
        List of excel sheets added or not added to the 'all_data' table in the database
        Adds the names of sheets that have been added to the 'added_sheets' table
    """
    # Loads excel sheets into the sqlite3 database 'apkr_reactions.db'
    # Makes separtate tables for each sheet, the code will combine them later
    try:
        conn = sqlite3.connect(database_name)
        data_files_added = 0
        for sheet in sheet_list:
            data_path = f'../0.database_excels/{sheet}'
            new_data_df = pd.read_excel(data_path)
            try:
                new_data_df.to_sql(f'{sheet}', conn, if_exists='fail', index=False)

            except ValueError as e:     # raises error if the table already exists in the database
                if "already exists" in str(e):
                    print(f"Skipping '{sheet}': a table with that name already exists in {database_name}.")
                    continue
                raise
            conn.commit()
            data_files_added += 1

    finally:
        print(f'Successfully added {data_files_added} to {database_name}')
        conn.close()

    # combines all newly generated tables into on table in the database
    # records all sheet names that have been added to the complete database
    conn = sqlite3.connect(database_name)
    cur = conn.cursor() # defines cursor
    try:
        # generates table that keeps track of data sheets that have been added
        cur.execute("""
            CREATE TABLE IF NOT EXISTS added_sheets (
                sheet_name TEXT PRIMARY KEY,
                rows_added INTEGER,
                added_at TEXT DEFAULT CURRENT_TIMESTAMP
            )
        """)
        # check which sheets have already been logged and skips them 
        already_added = {row[0] for row in
                        cur.execute("SELECT sheet_name FROM added_sheets").fetchall()}

        source_tables = [row[0] for row in cur.execute("""
            SELECT name FROM sqlite_master
            WHERE type='table' 
            AND name LIKE '%.xlsx' 
            AND name NOT IN ('all_data', 'added_sheets')
            """).fetchall()]

        for t in source_tables:
            if t in already_added:
                print(f"Skipping '{t}' -- already in 'all_data'")
                continue

            # how many rows need to be added?
            n_rows = cur.execute(f'SELECT COUNT(*) FROM "{t}"').fetchone()[0]

            # append rows into 'all_data' (create it from the first table if needed)
            table_exists = cur.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='all_data'"
            ).fetchone()
            if table_exists:
                cur.execute(f'INSERT INTO all_data SELECT * FROM "{t}"')
            else:
                cur.execute(f'CREATE TABLE all_data AS SELECT * FROM "{t}"')

            # log process 
            cur.execute(
                "INSERT INTO added_sheets (sheet_name, rows_added) VALUES (?, ?)",
                (t, n_rows)
            )
            try:
                # ... all the create / insert / combine work...
                conn.commit()
                print(f"Added '{t}' ({n_rows} rows)")
            except Exception as e: 
                print(f"Something went wrong: {e}")
                raise

    finally:
        conn.close()

#############################################################################################################
def check_if_updated(database_name='apkr_reactions.db'):
    """Checks if the database is completely updated. Runs update_database function
    if the database needs to be updated.
    """
    data_path = Path("../0.database_excels/")       # directory where excel sheets are located
    sheet_list = [f.name for f in data_path.iterdir() if f.is_file()]
    sheet_list.remove('template.xlsx')

    conn = sqlite3.connect(database_name)
    cur = conn.cursor()

    cur.execute("SELECT sheet_name FROM added_sheets")
    already_added = {row[0] for row in cur.fetchall()}
    conn.close()

    yet_to_upload = []
    for entry in sheet_list:
        if entry not in already_added:
            yet_to_upload.append(entry)

    if yet_to_upload:
        print("Not yet in database:", sorted(yet_to_upload))
        update_database(sheet_list=yet_to_upload)
    else:
        print("Good to go! All sheets are in database.")

    conn = sqlite3.connect(database_name)
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM all_data")
    reaction_num = cur.fetchone()[0]
    conn.close()
    print(f'{reaction_num} in database.')

################################################################################################################
# This makes more sense as a class
# Maybe fix this if you have time
# Function contains all troubleshooting scripts that have been written
def troubleshooting(type, sheet_list, database_name='apkr_reactions.db'):
    """Identifies the error more specificially.

    Args: 
        type: -- 'mismatching row count': OperationalError: table all_data has 27 columns but 26 values were supplied
              -- 'reload data': from an updated Excelsheet, removes the previous data and updates the database
        database_name: default is 'apkr_reactions.db'
        sheet_list: list of sheets that need to be troubleshooted

    Returns: 
        Useful troubleshooting information. You'll thank me. 
    """
    if type == 'mismatching row count':
        # OperationalError: table all_data has 27 columns but 26 values were supplied
        conn = sqlite3.connect(database_name)
        cur = conn.cursor()

        target_cols = [r[1] for r in cur.execute("PRAGMA table_info('all_data')")]
        print(f"all_data: {len(target_cols)} columns")
        print(target_cols)

        source_tables = [r[0] for r in cur.execute(
            "SELECT name FROM sqlite_master WHERE type='table' "
            "AND name NOT IN ('all_data', 'added_sheets')"
        )]

        for t in source_tables:
            cols = [r[1] for r in cur.execute(f"PRAGMA table_info('{t}')")]
            if cols != target_cols:
                print(f"\n{t}: {len(cols)} columns — MISMATCH")
                print(f"  missing from this sheet: {set(target_cols) - set(cols)}")
                print(f"  extra in this sheet:     {set(cols) - set(target_cols)}")

        conn.close()

    if type == 'reload data':
        # reloads data for a particular sheet so that all_data is updated
        conn = sqlite3.connect(database_name)
        for sheet in sheet_list:
            data_path = f'../0.database_excels/{sheet}'
            new_data_df = pd.read_excel(data_path)
            new_data_df.to_sql(f'{sheet}', conn, if_exists='replace', index=False)

        update_database(sheet_list)

# This should get added to the troubleshooting class 
# # pinpoints which file has the error
# conn = sqlite3.connect(database_name)
# conn.row_factory = lambda cursor, row: row[0]
# cur = conn.cursor()

# cur.execute('SELECT * FROM all_data WHERE "Substrate SMILES" = ?',
#             ('CC(=O)OC=C=C(C)CCCC#Cc1ccccc1  [(−)-9a, 98:2 er]',))
# for row in cur.fetchall():
#     print(row)

######################################################################################
def safe_ident(name):
    """SQLite can't paramaterize identifiers, so validate before interpolating.

    Returns: cleaned name without characters that SQLite cannot read (ex. "-")
    """
    cleaned = re.sub('r\W', '_', name)
    print(cleaned)
    if not re.fullmatch(r'[A-Za-z_]\w*', cleaned):
        raise ValueError(f'unusable table name: {name!r}')
    return cleaned

###########################################################################################
def make_all_data_copy(database_name='apkr_reactions.db'):
    """Makes a copy of the 'all_data' table as a back up for future use if required.
    """
    now = datetime.now()
    timestamp_string = now.strftime("%m.%d.%Y")
    print(timestamp_string)

    # Make a csv from the current database
    conn = sqlite3.connect(database_name)

    # makes 'all_data' table from database into a csv
    current_all_data_database = conn.execute(''' 
                                SELECT * FROM all_data
                                ''').fetchall()
    conn.close()

    data_df = pd.DataFrame(current_all_data_database)

    data_df.to_csv(f"../5.backup_all_data_databases/database_{timestamp_string}")






