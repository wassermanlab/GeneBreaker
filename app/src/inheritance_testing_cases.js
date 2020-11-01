import React from 'react';
import './home.css'
import Nav from './nav';

function InheritanceTestingCases(props) {
  let data = require('./inheritance_testing.json');

  return (
    <React.Fragment>
      <Nav />
      <div className="container navbar-fix">
        <div className="row">
          <div className="col-sm">
            <h1>Inheritance Testing Cases</h1>
            <p></p>

            <table className="table">
              <thead>
                <tr>
                  <th>Proband Name</th>
                  <th>Proband Sex</th>
                  <th>Inheritance</th>
                  <th>Gene</th>
                  <th>HPO terms</th>
                  <th>Reference Genome</th>
                  <th>Simulations (.vcf .vcfi .ped)</th>
                </tr>
              </thead>
              <tbody>
                {

                  data.map(({ patient, inheritance, gene, hpo, sex, genome, vcf, vcfi, ped, phenopacket_json }, i) => (
                    <tr>
                      <td>{patient}</td>
                      <td>{sex}</td>
                      <td>{inheritance}</td>
                      <td>{gene}</td>
                      <td>
                        {hpo.map((term, j) =>
                          <div><a href={"https://hpo.jax.org/app/browse/term/" + term}>{term}</a> <br /></div>)}
                      </td>
                      <td>{genome}</td>
                      <td>
                        <a href={vcf}>VCF</a><br/>
                        <a href={vcfi}>VCFI </a><br/>
                        <a href={ped}>PED </a><br/>
                        <a href={phenopacket_json}>Phenopacket JSON</a>
                      </td>
                    </tr>

                  ))}


              </tbody>
            </table>
          </div>
        </div>
      </div>
    </React.Fragment>)
}

export default InheritanceTestingCases;