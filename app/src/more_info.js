import React from 'react';
import './home.css'
import Nav from './nav';

function MoreInfo(props) {

  return (
    <React.Fragment>
      <Nav />
      <div className="container navbar-fix">
        <div className="row">
          <div className="col-sm">
            <h1>About</h1>
            <p>GeneBreaker was developed by Phillip Richmond and 
              Tamar Av-Shalom in Wyeth Wasserman’s lab at the BC Children’s 
              Hospital Research Institute, University of British Columbia. </p>
          </div>
        </div>
        <div className="row">
          <div className="col-sm">
            <h1>Manuscript</h1>
            <p>The manuscript describing GeneBreaker can be found at 
              [INSERT BIORXIV LINK WHEN AVAILABLE]</p>
          </div>
        </div>
        <div className="row">
          <div className="col-sm">
            <h1>Code and Data</h1>
            <p>GeneBreaker is an open source project, and the code underlying 
              the tool can be found at [INSERT NEW GENEBREAKER GITHUB] 
              Background datasets are hosted by the Zenodo repository and 
              can be found at: https://zenodo.org/record/3829960#.XsWb6S8ZNTY </p>
          </div>
        </div>
        <div className="row">
          <div className="col-sm">
            <h1>Downstream Benchmarkign</h1>
            <p>For more information about using GeneBreaker for downstream 
              benchmarking, visit the Github repository listed above. </p>
          </div>
        </div>
      </div>
      </React.Fragment >)
}

export default MoreInfo;