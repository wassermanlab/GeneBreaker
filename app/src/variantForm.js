import React from 'react';

// todo: move to separate component 
function GeneralErrors(props) {
  const select_error = "No transcript is selected.";
  const sex_error = "Transcript on Y chromosome cannot be selected with XX proband."
  if (props.errors.general.transcript_select) {
    return <div className="alert alert-danger" role="alert">{select_error}</div>
  }
  else if (props.errors.general.transcript_sex) {
    return <div className="alert alert-danger" role="alert">{sex_error}</div>
  }
  return null;
}

function NavButtons(props) {
  const page = props.page;
  if (page === 1) {
    return (<div><button type="button" className="btn btn-primary float-right">Next</button></div>)
  }
  else if (page === 2) {
    return (
      <div>
        <button type="button" className="btn btn-primary float-right">Next</button>
        <button type="button" className="btn btn-primary float-right">Back</button> </div>)
  }
  else if (page === 3) {
    return (
      <div>
        <button type="button" className="btn btn-primary float-right">Next</button>
        <button type="button" className="btn btn-primary float-right">Back</button> </div>)
  }
  return null;
}

class VFrom extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      page: 1,
      errors: {
        general: {
          transcript_select: false,
          transcript_sex: false
        }
      },
      transcript_list: JSON.parse(`[]`),
      clinvar_list: [],
      clingen_list: [],
      str_list: [],
      gene_uid: "",
      genome: "hg38",
      gene_name: "",
      chr: "",
      sex: "XX",
      var1: {
        type: "",
        region: "",
        zygosity: "",
        impact: {
          clinvar: {
            id: ""
          },
          cnv: {
            start: "",
            end: "",
            copy_change: ""
          },
          indel: {
            length: "",
            start: ""
          },
          mei: {
            element: "",
            start: ""
          },
          snv: {
            snv_type: "",
            start: ""
          },
          str: {
            start: "",
            end: "",
            motif: "",
            length: ""
          }
        }
      },
      var2: {
        type: "",
        region: "",
        zygosity: "",
        impact: {
          clinvar: {
            id: ""
          },
          cnv: {
            start: "",
            end: "",
            copy_change: ""
          },
          indel: {
            length: "",
            start: ""
          },
          mei: {
            element: "",
            start: ""
          },
          snv: {
            snv_type: "",
            start: ""
          },
          str: {
            start: "",
            end: "",
            motif: "",
            length: ""
          }
        }
      }
    };

    this.handleInputChange = this.handleInputChange.bind(this);
    this.getTranscripts = this.getTranscripts.bind(this);
  }

  handleInputChange(event) {
    const target = event.target;
    const value = target.value;
    const name = target.name;

    this.setState({
      [name]: value
    },
      () => console.log(this.state));
  }

  // Go over this stuff 
  getTranscripts(event) {
    this.setState({...this.state, isFetching: true})
    fetch('http://127.0.0.1:5001/get_transcripts/' + this.state.genome + '/' + this.state.gene_name)
      .then(response => {
        response.json()})
      .then(result => this.setState({transcript_list: result, 
                                     isFetching: false}))
      .catch(e => console.log(e));
  }

  render() {
    return (
      <form>
        {/* general: genome, sex, gene_uid, chr, transcript */}
        <React.Fragment>
          {/* genome */}
          <div className="form-group">
            <label>Genome</label>
            <select className="form-control"
              name="genome"
              value={this.state.genome}
              onChange={this.handleInputChange}>
              <option value="hg38">hg38</option>
              <option value="hg19">hg19</option>
            </select>
          </div>
          {/* sex */}
          <div className="form-group">
            <label>Genome</label>
            <select className="form-control"
              name="sex"
              value={this.state.sex}
              onChange={this.handleInputChange}>
              <option value="XX">XX</option>
              <option value="XY">XY</option>
            </select>
          </div>
          {/* gene name/symbol */}
          <label>Gene symbol</label>
          <div className="input-group mb-3">
            <input type="text" className="form-control" name="gene_name" value={this.state.gene_name} onChange={this.handleInputChange} />
            <div className="input-group-append">
              <button className="btn btn-outline-secondary" type="button" onClick={this.getTranscripts}>Get Transcripts</button>
            </div>
          </div>
          {/* transcript list */}
          <div className="form-group">
            <label>Transcript</label>
            <select className="form-control"
              name="gene_uid"
              value={this.state.gene_uid}
              onChange={this.handleInputChange}
              size="5">
              <option key="0" value=""></option>
              {this.state.transcript_list.map((item, index) => (
                <option key={item.qualifiers.uid} value={item.qualifiers.uid} chrom={item.chrom}> {item.qualifiers.name} </option>
              ))}
            </select>
          </div>
          <GeneralErrors errors={this.state.errors} />
          <NavButtons page={this.state.page} />
        </React.Fragment>
        {/* var1 */}
        {/* var2 */}




      </form>
    );
  }
}

export default VFrom;
