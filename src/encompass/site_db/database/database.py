# database.py
import keyring
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, declarative_base
from dotenv import load_dotenv, find_dotenv


load_dotenv(find_dotenv())

# Get non-sensitive config from .env
db_user = os.getenv('SQL_DB_USER')
db_host = os.getenv('SQL_DB_HOST')
db_port = os.getenv('SQL_DB_PORT')
db_name = os.getenv('SQL_DB_NAME')

# Get password from keyring
db_password = keyring.get_password("encompass_sql_db", "db_password")

DATABASE_URL = f"postgresql://{db_user}:{db_password}@{db_host}:{db_port}/{db_name}"

Base = declarative_base()

engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(bind=engine)
