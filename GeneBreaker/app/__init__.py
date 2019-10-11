from flask import Flask
UPLOAD_FOLDER = '/tmp/simulator/'

app = Flask(__name__)
app.config.from_mapping(
    SECRET_KEY='dev',
    UPLOAD_FOLDER = UPLOAD_FOLDER
)
import GeneBreaker.app.routes

