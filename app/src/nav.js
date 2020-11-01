import React from 'react';
import { Link } from "react-router-dom";
function Nav(props) {

  return (<React.Fragment>
    <nav className="navbar fixed-top navbar-expand-lg navbar-light bg-light">
      <Link className="navbar-brand" to="/">
      <img src={'/genebreaker_logo_full.svg'} width="200"  alt="GeneBreaker" />
      </Link>
      <button className="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarNav" aria-controls="navbarNav" aria-expanded="false" aria-label="Toggle navigation">
        <span className="navbar-toggler-icon"></span>
      </button>
      <div className="collapse navbar-collapse" id="navbarNav">
        <ul className="navbar-nav">
          <li className="nav-item active">
            <Link className="nav-link" to="/variants">
              Variant Designer
            </Link>
          </li>
          <li className="nav-item active">
            <Link className="nav-link" to="/inheritance_testing_cases">
              Inheritance Testing Cases
            </Link>
          </li>
          <li className="nav-item active">
            <Link className="nav-link" to="/patient_scenarios">
              Patient Scenarios
            </Link>
          </li>
          <li className="nav-item active">
            <Link className="nav-link" to="/more_info">
              More Info
            </Link>
          </li>
        </ul>
      </div>
    <span className="navbar-text mr-auto">
      beta
    </span>

    </nav>



  </React.Fragment>)
}

export default Nav;