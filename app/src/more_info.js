import React from 'react';
import './home.css'
import Nav from './nav';

function MoreInfo(props) {

  return (
    <React.Fragment>
      <Nav />
      <div className="container navbar-fix">
        <div className="row row-padding">
          <div className="col-sm">
            <h1>Video Tutorial</h1>
            <div className="embed-responsive embed-responsive-16by9">
              <iframe title="youtubeTutorial" className="embed-responsive-item" src="https://www.youtube.com/embed/_jzI4WY8rUc" frameBorder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture;" allowFullScreen ></iframe>
            </div>
          </div>
        </div>
        <div className="row row-padding">
          <div className="col-sm">
            <h1>About</h1>
            <p>GeneBreaker was developed by Phillip Richmond and
            Tamar Av-Shalom in Wyeth Wasserman’s lab at the BC Children’s
              Hospital Research Institute, University of British Columbia. </p>
          </div>
        </div>
        <div className="row row-padding">
          <div className="col-sm">
            <h1>Manuscript</h1>
            <p>The manuscript describing GeneBreaker can be found on 
              <a href="https://www.biorxiv.org/content/10.1101/2020.05.29.124495v1"> bioRxiv.</a></p>
          </div>
        </div>
        <div className="row row-padding">
          <div className="col-sm">
            <h1>Code and Data</h1>
            <p>GeneBreaker is an open source project, and the code underlying
            the tool can be found on 
            <a href="https://github.com/wassermanlab/GeneBreaker"> GitHub.</a> 
             Please submit any code issues directly on our GitHub Repo. 
             Background datasets are hosted by
              <a href="https://zenodo.org/record/3829960#.XsWb6S8ZNTY"> Zenodo.</a> </p>
          </div>
        </div>
        <div className="row row-padding">
          <div className="col-sm">
            <h1>Downstream Benchmarking</h1>
            <p>For more information about using GeneBreaker for downstream
              benchmarking, visit the Github repository listed above. </p>
          </div>
        </div>
      </div>
    </React.Fragment >)
}

export default MoreInfo;