import os

from flask import Flask
from flask_bootstrap import Bootstrap

app = Flask(__name__) 
boostrap = Bootstrap(app)

from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
 
def establish_GUD_session(
        db_user="ontarget_r",
        db_pass=None,
        db_host="ontarget.cmmt.ubc.ca",
        db_port=5506,
        db_name="tamar_test"
    ):

        if not db_pass: db_pass = ""
    
        gud_db = "mysql://{}:{}@{}:{}/{}".format(
            db_user,
            db_pass,
            db_host,
            db_port,
            db_name
        )

        # Establish a MySQL session
        try:
            engine = create_engine(
                gud_db,
                echo=False,
                pool_pre_ping=True
            )
            session = Session(engine)
        except:
            raise ValueError(
                "Could not connect to GUD db: %s" \
                % gud_db
            )

        return session

from . import routes

# def create_app(test_config=None):
#     # create and configure the app
#     app = Flask(__name__, instance_relative_config=True)
#     app.config.from_mapping(
#         SECRET_KEY='dev',
#     )

#     if test_config is None:
#         # load the instance config, if it exists, when not testing
#         app.config.from_pyfile('config.py', silent=True)
#     else:
#         # load the test config if passed in
#         app.config.from_mapping(test_config)

#     # ensure the instance folder exists
#     try:
#         os.makedirs(app.instance_path)
#     except OSError:
#         pass

#     # a simple page that says hello
#     @app.route('/hello')
#     def hello():
#         return 'Hello, World!'

#     return app
