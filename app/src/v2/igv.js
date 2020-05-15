import React, { Component } from 'react';
import igv from 'igv';

class IGV extends Component {
  constructor(props) {
    super(props);
    this.state = { browser: null };
    this.reset = this.reset.bind(this)

  }

  getGenome(g) {
    if (g === "grch38") {
      return "hg38";
    }
    else {
      return "hg19";
    }
  }

  componentDidMount() {
    let currentComponent = this;
    const g = this.getGenome(this.props.genome);
    let igvOptions = { genome: g };
    let igvContainer = document.getElementById('igv-div');
    igv.createBrowser(igvContainer, igvOptions)
      .then(function (browser) {
        currentComponent.setState({ browser: browser })
      })
  }

  componentDidUpdate(prevProps) {
    // Typical usage (don't forget to compare props):
    const g = this.getGenome(this.props.genome)
    if (this.props.genome !== prevProps.genome) {
      if (this.state.browser !== null) {
        this.state.browser.loadGenome({ genome: g });
      }
    }
    if (this.props.start !== prevProps.start) {
      if (this.state.browser !== null) {
        const locus = "chr" + this.props.chrom + ":" + this.props.start + "-" + this.props.end;
        this.state.browser.search(locus);
      }
    }
  }

  reset() {
    if (this.state.browser !== null) {
      const locus = "chr" + this.props.chrom + ":" + this.props.start + "-" + this.props.end;
      this.state.browser.search(locus);
    }
  }

  render() {
    return (
      <div>
        <div id="igv-div"></div>
          <small className="form-text text-muted">
          <button onClick={this.reset} type="button" className="btn btn-light btn-sm" style={{marginRight: "10px"}}>
            Reset view on transcript</button>
            Click on the first track's settings cog to see the three frame translation track.
      </small>
      </div>
    );
  }
}
export default IGV;
