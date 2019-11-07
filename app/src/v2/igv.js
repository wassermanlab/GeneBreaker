import React, { Component } from 'react';
import igv from 'igv';

class IGV extends Component {
  constructor(props) {
    super(props);
    this.state = {browser: null};
  }

  componentDidMount() {
    let currentComponent = this;
    let igvOptions = {genome: this.props.genome};
    // let igvOptions = { genome: 'hg19' };
    let igvContainer = document.getElementById('igv-div');
    igv.createBrowser(igvContainer, igvOptions)
      .then(function (browser) {
        currentComponent.setState({browser: browser})
      })
  }

  componentDidUpdate(prevProps) {
    // Typical usage (don't forget to compare props):
    if (this.props.genome !== prevProps.genome) {
      if (this.browser !== null){
        this.state.browser.loadGenome({genome: this.props.genome});
      }
    }
    if (this.props.start !== prevProps.start) {
      console.log(this.props)
      if (this.browser !== null){
        const locus = "chr"+this.props.chrom + ":" + this.props.start + "-" + this.props.end;
        console.log(locus);
        this.state.browser.search(locus);
      }
    }
  }

  render() {
    return (
      <div id="igv-div"></div>
    );
  }
}
export default IGV;
