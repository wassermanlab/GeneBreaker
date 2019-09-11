from flask import Flask
from flask_cors import CORS
UPLOAD_FOLDER = '/tmp/simulator/'

app = Flask(__name__)
app.config.from_mapping(
    SECRET_KEY='dev',
    UPLOAD_FOLDER = UPLOAD_FOLDER
)
CORS(app)
import MenDelSIM.app.routes

