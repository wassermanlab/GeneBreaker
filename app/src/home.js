import React from 'react';
import './home.css'
import { Link } from "react-router-dom";
import Nav from './nav';

function Home(props) {

  return (
    <div>
      <Nav/>
      <div className="center" width="100%">
        <h1 className="header-title">GeneBreaker</h1>
        <p className="header-subtitle">design variants 	&middot; build inheritance &middot; simulate cases</p>
        <Link to="/variants">
          <button type="button" className="btn btn-outline-dark">Start Designing</button>
        </Link>
        <img src='/genebreaker_infographic.png'
         style={{
          paddingTop: '30px',
          paddingbottom: '30px',
        }}
         width="90%" alt="GeneBreaker" />
      </div>
    </div>)
}

export default Home;