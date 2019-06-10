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