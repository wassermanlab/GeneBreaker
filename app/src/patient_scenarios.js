import React from 'react';
import './home.css'
import Nav from './nav';

function PatientScenarios(props) {
  let data = require('./PatientScenarios.json');

  return (
    <React.Fragment>
      <Nav />
      <div className="container navbar-fix">
        <div className="row">
          <div className="col-sm">
            <h1>Patient Scenarios</h1>
            {
              data.map(({ title, clinical_family_history, hpo, ped, phenopacket_json, vcf }, i) => (

                <div class="accordion" id="accordionExample">
                  <div class="card">
                    <div class="card-header" id="headingOne">
                      <h2 class="mb-0">
                        <button class="btn btn-link btn-block text-left" type="button" data-toggle="collapse" data-target={"#collapse" + i} aria-expanded="true" aria-controls={"#collapse" + i}>
                          {title}
                        </button>
                      </h2>
                    </div>

                    <div id={"collapse" + i} class="collapse" aria-labelledby="headingOne" data-parent="#accordionExample">
                      <div class="card-body">
                        <h4>Clinical and Family History</h4>
                        <p>{clinical_family_history}</p>
                        <br></br>
                        <h4>HPO Terms</h4>
                        {hpo.map(({ term, id }, j) =>
                          <div><a href={"https://hpo.jax.org/app/browse/term/" + id}>{term}</a> <br /></div>)}
                        <br></br>
                        <h4>Links</h4>
                        <div><a href={ped}>PED file</a> <br /></div>
                        <div><a href={phenopacket_json}>Phenopacket JSON</a> <br /></div>
                        <div><a href={vcf}>VCF</a> <br /></div>
                      </div>
                    </div>
                  </div>
                </div>
              ))}
          </div>
        </div>
      </div>
    </React.Fragment>)
}

export default PatientScenarios;