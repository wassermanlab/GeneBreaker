from flask import Flask
from flask_cors import CORS
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

app = Flask(__name__)
app.config.from_pyfile('config.py')
limiter = Limiter(
    app,
    key_func=get_remote_address,
    default_limits=["500 per day", "5 per second"]
)
cors = CORS(app, resources={r"/*": 
{"origins": ["http://localhost:3000", "http://genebreaker.cmmt.ubc.ca"]}})


import GeneBreaker.app.routes

