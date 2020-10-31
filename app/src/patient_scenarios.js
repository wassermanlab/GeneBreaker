import React from 'react';
import './home.css'
import Nav from './nav';

function PatientScenarios(props) {
  let data = require('./SimulationScenarioInformation.json');

  return (
    <React.Fragment>
      <Nav />
      <div className="container navbar-fix">
        <div className="row">
          <div className="col-sm">
            <h1>Patient Scenarios</h1>
            <p>Under development</p>

          </div>
        </div>
      </div>
    </React.Fragment>)
}

export default PatientScenarios;